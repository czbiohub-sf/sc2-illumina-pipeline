def helpMessage() {
  log.info"""
    Pipeline for running augur commands to build auspice visualization from main output.

    Usage:

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
      --sample_sequences            Glob pattern of sequences
      --ref                         Reference FASTA file (default: data/MN908947.3.fa)
      --ref_gb                      Reference Genbank file for augur (default: data/MN908947.3.gb)
      --blast_sequences             FASTA of sequences for BLAST alignment
      --nextstrain_sequences        FASTA of sequences to build a tree with
      --minLength                   Minimum base pair length to allow assemblies to pass QC (default: 29000)
      --maxNs                       Max number of Ns to allow assemblies to pass QC (default: 100)

    Optional arguments:
      --sample_metadata             TSV of metadata from main output
      --clades                      TSV with clades from nextstrain (default: data/clades.tsv)
      --sample_vcfs                 Glob pattern of corresponding VCF files


    Nextstrain options:
      --nextstrain_ncov             Path to nextstrain/ncov directory (default: fetches from github)
      --group_by                    string, parameter to augur filter (default: 'division year month')
      --sequences_per_group_1       Initial subsampling (default: 500)
      --sequences_per_group_2       Contextual subsampling by priority (default: 20)
      --existing_alignment          Add to existing alignment in augur align step

    Other options:
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      -resume                       Use cached results

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

ref_fasta = file(params.ref, checkIfExists: true)
ref_gb = params.ref_gb ? file(params.ref_gb, checkIfExists: true) : Channel.empty()


Channel
  .fromPath(params.sample_sequences)
  .map {file -> tuple(file.simpleName, file) }
  .into {blastconsensus_in; realign_fa; stats_fa}

Channel
  .fromPath(params.sample_sequences)
  .set {merge_fastas_ch}



process realignConsensus {
    tag { sampleName }
    publishDir "${params.outdir}/realigned-seqs", mode: 'copy'

    input:
    tuple(sampleName, path(in_fa)) from realign_fa
    path(ref_fasta)

    output:
    tuple(sampleName, path("${sampleName}.realigned.bam")) into realigned_bam
    path("${sampleName}.realigned.bam.bai")

    script:
    """
    minimap2 -ax asm5 -R '@RG\\tID:${sampleName}\\tSM:${sampleName}' \
      ${ref_fasta} ${in_fa} |
      samtools sort -O bam -o ${sampleName}.realigned.bam
    samtools index ${sampleName}.realigned.bam
    """
}

realigned_bam.into { call_variants_bam; combined_variants_bams }
combined_variants_bams = combined_variants_bams.map { it[1] }.collect()

// Allow VCF files to be optional in case we want to include assemblies that aren't produced by the core
// consensus calling pipeline
if (params.sample_vcfs) {
  Channel
    .fromPath(params.sample_vcfs)
    .map {file -> tuple(file.simpleName, file) }
    .set {variants_ch}
} 
else {
  process callVariants {
    tag { sampleName }
    publishDir "${params.outdir}/sample-variants", mode: 'copy'

    input:
    tuple(sampleName, path(in_bam)) from call_variants_bam
    path(ref_fasta)

    output:
    tuple(sampleName, path("${sampleName}.vcf")) into variants_ch
    path("${sampleName}.bcftools_stats") into bcftools_stats_ch

    script:
    """
    bcftools mpileup -f ${ref_fasta} ${in_bam} |
      bcftools call --ploidy 1 -m -P ${params.bcftoolsCallTheta} -v - \
      > ${sampleName}.vcf
    bcftools stats ${sampleName}.vcf > ${sampleName}.bcftools_stats
    """
  }
}

variants_ch.into {primer_variants_ch; assignclades_in; variants_ch}

clades = file(params.clades, checkIfExists: true)

process assignClades {
    // Use Nextstrain definitions to assign clades based on mutations

    input:
    tuple(sampleName, path(vcf)) from assignclades_in
    path(ref_gb)
    path(clades)

    output:
    tuple(sampleName, path("${sampleName}.clades")) into assignclades_out

    script:
    """
    assignclades.py \
        --reference ${ref_gb} --clades ${clades} \
        --vcf ${vcf} --sample ${sampleName}
    """
}

if (params.qpcr_primers) {
    qpcr_primers = file(params.qpcr_primers, checkIfExists: true)
} else {
    qpcr_primers = Channel.empty()
}

process searchPrimers {
    publishDir "${params.outdir}/primer-variants", mode: 'copy'

    input:
    tuple(sampleName, path(vcf)) from primer_variants_ch
    path(qpcr_primers)

    output:
    tuple(sampleName, path("${sampleName}_primers.vcf")) into primer_variants_vcf
    path("${sampleName}_primers.primer_variants_stats") into primer_stats_ch

    when:
    params.qpcr_primers

    script:
    """
    bgzip ${vcf}
    bcftools index ${vcf}.gz
    bcftools view -R ${qpcr_primers} -o ${sampleName}_primers.vcf ${vcf}.gz
    bcftools stats ${sampleName}_primers.vcf > ${sampleName}_primers.primer_variants_stats
    """
}


blast_sequences = params.blast_sequences ? file(params.blast_sequences, checkIfExists: true) : Channel.empty()

process buildBLASTDB {

    input:
    path(blast_sequences)

    output:
    path("blast_seqs.nt*") into blastdb_ch
    path("blast_sequences.fasta") into blast_clean_ch

    script:
    """
    normalize_gisaid_fasta.sh ${blast_sequences} blast_seqs_clean.fasta
    sed 's/\\.//g' blast_seqs_clean.fasta > blast_sequences.fasta

    makeblastdb -in blast_sequences.fasta -parse_seqids -title 'blastseqs' -dbtype nucl -out blast_seqs.nt
    """
}

process blastConsensus {
    tag {sampleName}
    publishDir "${params.outdir}/samples/${sampleName}", mode: 'copy'

    input:
    tuple(sampleName, path(assembly)) from blastconsensus_in
    path(dbsequences) from blast_clean_ch
    path(blastdb) from blastdb_ch
    path(ref_fasta)

    output:
    path("${sampleName}.blast.tsv")
    tuple(sampleName, path("${sampleName}_nearest_blast.fasta")) into (nearest_neighbor, collectnearest_in)

    script:
    """
    get_top_hit.py --minLength ${params.minLength} \
      --sequences ${dbsequences} \
      --sampleName ${sampleName} \
      --assembly ${assembly} \
      --default ${ref_fasta}
    """
}

process collectNearest {
    publishDir "${params.outdir}/BLAST", mode: 'copy'

    input:
    path(fastas) from collectnearest_in.map{it[1]}.collect()

    output:
    path("included_samples.fasta") into (included_fastas_ch, contextual_fastas_ch)

    script:
    """
    cat ${fastas} > all_included_samples.fasta
    seqkit rmdup all_included_samples.fasta > deduped_included_samples.fasta
    seqkit faidx -f -r deduped_included_samples.fasta '^[^MN908947.3].*' > included_samples.fasta
    """
}

stats_fa
  .join(variants_ch)
  .join(primer_variants_vcf)
  .join(assignclades_out)
  .join(nearest_neighbor)
  .set{stats_ch_in}

process computeStats {
    tag { sampleName }
    publishDir "${params.outdir}/coverage-plots", mode: 'copy',
        saveAs: { x -> x.endsWith(".png") ? x : null }

    input:
    tuple(sampleName,
          file(in_fa),
          file(vcf),
          file(primer_vcf),
          file(in_clades),
          file(neighbor_fasta)) from stats_ch_in

    output:
    path("${sampleName}.stats.json") into stats_ch

    script:
    """
    alignment_assembly_stats.py \
        --sample_name ${sampleName} \
        --assembly ${in_fa} \
        --vcf ${vcf} \
        --primervcf ${primer_vcf} \
        --neighborfasta ${neighbor_fasta} \
        --clades ${in_clades} \
        --out_prefix ${sampleName}
    """
}

process combinedVariants {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(in_bams) from combined_variants_bams
    path(ref_fasta)

    output:
    path("combined.vcf") into combined_variants_vcf

    script:
    """
    bcftools mpileup -f ${ref_fasta} ${in_bams} |
      bcftools call --ploidy 1 -m -P ${params.bcftoolsCallTheta} -v - \
      > combined.vcf
    """
}

process mergeAllAssemblies {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(in_fasta) from merge_fastas_ch.collect()

    output:
    path("combined.fa") into merged_assemblies_ch

    script:
    """
    cat ${in_fasta} > combined.fa
    """
}

process mergeAssemblyStats {
    publishDir "${params.outdir}/run_analysis-stats", mode: 'copy'

    input:
    path(in_json) from stats_ch.collect()

    output:
    path("combined.stats.tsv") into merged_stats_ch

    script:
    """
    merge_stats.py analysis ${in_json} > combined.stats.tsv
    """
}

process filterAssemblies {
    publishDir "${params.outdir}", mode: 'copy',
      saveAs: {x -> x.endsWith(".tsv") ? "run_analysis-stats/$x" : x}

    input:
    path(merged_stats) from merged_stats_ch
    path(merged_assemblies) from merged_assemblies_ch
    path(vcf) from combined_variants_vcf

    output:
    path("filtered.stats.tsv")
    path("filtered.fa") into nextstrain_ch
    path("filtered.vcf")

    script:
    """
    filter_assemblies.py \
        --vcf ${vcf} \
        --max_n ${params.maxNs} --min_len ${params.minLength} \
        --stats ${merged_stats} --fasta ${merged_assemblies} \
        --out_prefix filtered
    """
}

nextstrain_ch.into {nextstrain_ch; sample_sequences_ch}

sample_sequences_ch
  .mix(contextual_fastas_ch)
  .collect()
  .set {sample_and_contextual_ch}
// Setup nextstrain files

if (params.nextstrain_sequences && params.nextstrain_ncov) {
  ref_gb = file(params.ref_gb, checkIfExists: true)
  nextstrain_sequences = file(params.nextstrain_sequences, checkIfExists: true)

  nextstrain_ncov = params.nextstrain_ncov
  if (nextstrain_ncov[-1] != "/") {
      nextstrain_ncov = nextstrain_ncov + "/"
  }

  nextstrain_metadata_path = file(nextstrain_ncov + "data/metadata.tsv", checkIfExists: true)
  nextstrain_config = nextstrain_ncov + "config/"
  include_file = file(nextstrain_config + "include.txt", checkIfExists: true)
  exclude_file = file(nextstrain_config + "exclude.txt", checkIfExists: true)
  clades = file(nextstrain_config + "clades.tsv", checkIfExists: true)
  auspice_config = file(nextstrain_config + "auspice_config.json", checkIfExists: true)
  lat_longs = file(nextstrain_config + "lat_longs.tsv", checkIfExists: true)


  sample_metadata = params.sample_metadata ? file(params.sample_metadata, checkIfExists: true) : Channel.empty()
}

if (params.sample_metadata) {
  process combineNextstrainInputs {
      publishDir "${params.outdir}/nextstrain/data", mode: 'copy'
      stageInMode 'copy'

      input:
      path(sample_sequences) from nextstrain_ch
      path(nextstrain_sequences)
      path(nextstrain_metadata_path)
      path(included_contextual_fastas) from included_fastas_ch
      path(include_file)
      path(sample_metadata)

      output:
      path('metadata.tsv') into (nextstrain_metadata, firstfilter_metadata, extractsamples_metadata, priorities_metadata, refinetree_metadata, infertraits_metadata, tipfreq_metadata, export_metadata)
      path('deduped_sequences.fasta') into (nextstrain_in, firstfilter_in)
      path('included_sequences.txt') into (nextstrain_include, firstfilter_include)
      path("internal_samples.txt") into sample_ids
      path("external_samples.txt") into external_ids
      path("all_sequences.fasta")

      script:
      // Normalize the GISAID names using Nextstrain's bash script
      """
      cat ${included_contextual_fastas} ${nextstrain_sequences} > contextual_sequences.fasta
      normalize_gisaid_fasta.sh contextual_sequences.fasta normalized_sequences.fasta
      cat ${included_contextual_fastas} | grep '>' | awk -F '>' '{print \$2}' > included_nearest.txt
      cat included_nearest.txt ${include_file} > included_sequences.txt
      cat normalized_sequences.fasta | grep '>' | awk -F '>' '{print \$2}' > external_samples.txt
      cat ${sample_sequences} | grep '>' | awk -F '>' '{print \$2}' > internal_samples.txt
      cat internal_samples.txt >> included_sequences.txt
      seqkit rmdup normalized_sequences.fasta > sequences.fasta

      make_nextstrain_input.py --prev_metadata ${nextstrain_metadata_path} \
          --prev_sequences sequences.fasta \
          --new_sequences ${sample_sequences} \
          --new_metadata ${sample_metadata}

      seqkit rmdup all_sequences.fasta > deduped_sequences.fasta
      """
  }
}
else {
  process makeNextstrainInputs {
    publishDir "${params.outdir}/nextstrain/data", mode: 'copy'
    stageInMode 'copy'

      input:
      path(sample_sequences) from nextstrain_ch
      path(nextstrain_sequences)
      path(nextstrain_metadata_path)
      path(included_contextual_fastas) from included_fastas_ch
      path(include_file)

      output:
      path('metadata.tsv') into (nextstrain_metadata, firstfilter_metadata, extractsamples_metadata, priorities_metadata, refinetree_metadata, infertraits_metadata, tipfreq_metadata, export_metadata)
      path('deduped_sequences.fasta') into firstfilter_in
      path('included_sequences.txt') into (nextstrain_include, firstfilter_include)
      path("internal_samples.txt") into sample_ids
      path("external_samples.txt") into external_ids
      path("all_sequences.fasta")

      script:
      currdate = new java.util.Date().format('yyyy-MM-dd')
      // Normalize the GISAID names using Nextstrain's bash script
      """
      cat ${included_contextual_fastas} ${nextstrain_sequences} > contextual_sequences.fasta
      normalize_gisaid_fasta.sh contextual_sequences.fasta normalized_sequences.fasta
      cat ${included_contextual_fastas} | grep '>' | awk -F '>' '{print \$2}' > included_nearest.txt
      cat included_nearest.txt ${include_file} > included_sequences.txt
      cat normalized_sequences.fasta | grep '>' | awk -F '>' '{print \$2}' > external_samples.txt
      cat ${sample_sequences} | grep '>' | awk -F '>' '{print \$2}' > internal_samples.txt
      cat internal_samples.txt >> included_sequences.txt
      seqkit rmdup normalized_sequences.fasta > sequences.fasta

      make_nextstrain_input.py --prev_sequences sequences.fasta \
          --prev_metadata ${nextstrain_metadata_path} \
          --new_sequences ${sample_sequences} \
          --date $currdate \
          --region 'North America' \
          --country USA \
          --division 'California' \
          --location 'San Francisco' \
          --submitting_lab 'Biohub' \
          --date_submitted $currdate \
          --date $currdate 

      seqkit rmdup all_sequences.fasta > deduped_sequences.fasta
      """
  }
}

process firstFilter {
  publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

  input:
  path(sequences) from firstfilter_in
  path(metadata) from firstfilter_metadata
  path(exclude_file)
  path(include_file) from firstfilter_include

  output:
  path("filtered_sequences.fasta") into firstfiltered_ch

  script:
  String exclude_where = "date='2020' date='2020-01-XX' date='2020-02-XX' date='2020-03-XX' date='2020-04-XX' date='2020-01' date='2020-02' date='2020-03' date='2020-04'"
  """
  augur filter \
            --sequences ${sequences} \
            --metadata ${metadata} \
            --include ${include_file} \
            --exclude ${exclude_file} \
            --exclude-where ${exclude_where}\
            --min-length ${params.minLength} \
            --group-by ${params.group_by} \
            --sequences-per-group ${params.sequences_per_group_1} \
            --output filtered_sequences.fasta
  """
}

if (params.existing_alignment) {
  existing_alignment = file(params.existing_alignment, checkIfExists: true)
} else {
  // set existing_alignment to something nonempty so process will run
  existing_alignment = file(params.ref_fasta, checkIfExists: true)
}

process alignSequences {
  label "process_large"
  publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

  input:
  path(sequences) from firstfiltered_ch
  path(ref_gb)
  path(existing_alignment)
  path(sample_sequences) from sample_and_contextual_ch

  output:
  path("aligned_raw.fasta") into (firstaligned_ch, makepriorities_ch, filterstrains_in)

  script:
  if (params.existing_alignment)
  """
  cat ${sample_sequences} > sample_and_contextual.fasta
  augur align \
            --sequences sample_and_contextual.fasta \
            --reference-sequence ${ref_gb} \
            --output aligned_raw.fasta \
            --nthreads ${task.cpus} \
            --remove-reference \
            --fill-gaps \
            --existing-alignment ${existing_alignment}
  """
  else
  """
  augur align \
            --sequences ${sequences} \
            --reference-sequence ${ref_gb} \
            --output aligned_raw.fasta \
            --nthreads ${task.cpus} \
            --remove-reference \
            --fill-gaps
  """

}

process extractSampleSequences {
  publishDir "${params.outdir}/nextstrain/results", mode: 'copy'


  input:
  path(sequences) from firstaligned_ch
  path(metadata) from extractsamples_metadata
  path(include) from sample_ids
  path(exclude) from external_ids

  output:
  path("aligned_samples.fasta") into aligned_samples_ch

  script:
  """
  augur filter \
            --sequences ${sequences} \
            --metadata ${metadata} \
            --include ${include} \
            --exclude ${exclude} \
            --output aligned_samples.fasta \
  """
}

process makePriorities {
  publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

  input:
  path(sample_sequences) from aligned_samples_ch
  path(aligned_sequences) from makepriorities_ch
  path(metadata) from priorities_metadata

  output:
  path("priorities.tsv") into priorities_tsv

  script:
  """
  priorities.py --alignment ${aligned_sequences} \
                --metadata ${metadata} \
                --focal-alignment ${sample_sequences} \
                --output priorities.tsv
  """
}

process filterStrains {
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(sequences) from filterstrains_in
    path(metadata) from nextstrain_metadata
    path(include_file) from nextstrain_include
    path(exclude_file)
    path(priorities) from priorities_tsv

    output:
    path('filtered_aligned.fasta') into (aligned_ch, refinetree_alignment, ancestralsequences_alignment)

    script:
    """
    augur filter \
            --sequences ${sequences} \
            --metadata ${metadata} \
            --include ${include_file} \
            --exclude ${exclude_file} \
            --priority ${priorities} \
            --group-by ${params.group_by} \
            --sequences-per-group ${params.sequences_per_group_2} \
            --min-length ${params.minLength} \
            --output filtered_aligned.fasta \
    """
}

process buildTree {
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    cpus 4

    input:
    path(alignment) from aligned_ch

    output:
    path("tree_raw.nwk") into tree_raw_ch

    script:
    """
    augur tree \
        --method fasttree \
        --alignment ${alignment} \
        --output tree_raw.nwk \
        --nthreads ${task.cpus}
    """
}

process refineTree {
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(tree) from tree_raw_ch
    path(alignment) from refinetree_alignment
    path(metadata) from refinetree_metadata

    output:
    path('tree.nwk') into (ancestralsequences_tree, translatesequences_tree, infertraits_tree, addclades_tree, tipfreq_tree, export_tree)
    path('branch_lengths.json') into export_branch_lengths

    script:
    """
    augur refine \
        --tree ${tree} \
        --alignment ${alignment} \
        --metadata ${metadata} \
        --output-tree tree.nwk \
        --output-node-data branch_lengths.json \
        --root 'Wuhan-Hu-1/2019' \
        --timetree \
        --clock-rate 0.0008 \
        --clock-std-dev 0.0004 \
        --coalescent skyline \
        --date-inference marginal \
        --divergence-unit mutations \
        --date-confidence \
        --no-covariance \
        --clock-filter-iqd 4
    """
}

process ancestralSequences {
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(tree) from ancestralsequences_tree
    path(alignment) from ancestralsequences_alignment

    output:
    path('nt_muts.json') into (translatesequences_nodes, addclades_nuc_muts, export_nt_muts)

    script:
    """
    augur ancestral \
        --tree ${tree} \
        --alignment ${alignment} \
        --output-node-data nt_muts.json \
        --inference joint \
        --infer-ambiguous
    """
}

process translateSequences {
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(tree) from translatesequences_tree
    path(nodes) from translatesequences_nodes
    path(ref_gb)

    output:
    path('aa_muts.json') into (addclades_aa_muts, export_aa_muts)

    script:
    """
    augur translate \
        --tree ${tree} \
        --ancestral-sequences ${nodes} \
        --reference-sequence ${ref_gb} \
        --output-node-data aa_muts.json \
    """
}


process inferTraits {
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(tree) from infertraits_tree
    path(metadata) from infertraits_metadata

    output:
    path('traits.json') into export_traits

    script:
    """
    augur traits \
        --tree ${tree} \
        --metadata ${metadata} \
        --output traits.json \
        --columns country_exposure \
        --confidence \
        --sampling-bias-correction 2.5 \
    """
}

process addClades {
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(tree) from addclades_tree
    path(aa_muts) from addclades_aa_muts
    path(nuc_muts) from addclades_nuc_muts
    path(clades)

    output:
    path('clades.json') into export_clades

    script:
    """
    augur clades --tree ${tree} \
        --mutations ${nuc_muts} ${aa_muts} \
        --clades ${clades} \
        --output-node-data clades.json
    """
}

process tipFrequencies {
    publishDir "${params.outdir}/nextstrain/auspice", mode: 'copy'

    input:
    path(tree) from tipfreq_tree
    path(metadata) from tipfreq_metadata

    output:
    path('ncov_tip-frequencies.json')

    script:
    """
    augur frequencies \
        --method kde \
        --metadata ${metadata} \
        --tree ${tree} \
        --min-date 2020.0 \
        --pivot-interval 1 \
        --narrow-bandwidth 0.05 \
        --proportion-wide 0.0 \
        --output ncov_tip-frequencies.json
    """
}

process exportData {
    publishDir "${params.outdir}/nextstrain/auspice", mode: 'copy'

    input:
    path(tree) from export_tree
    path(metadata) from export_metadata
    path(branch_lengths) from export_branch_lengths
    path(nt_muts) from export_nt_muts
    path(aa_muts) from export_aa_muts
    path(traits) from export_traits
    path(auspice_config)
    path(lat_longs)
    path(clades) from export_clades

    output:
    path('ncov.json') into extractvariants_in

    script:
    """
    augur export v2 \
        --tree ${tree} \
        --metadata ${metadata} \
        --node-data ${branch_lengths} ${nt_muts} ${aa_muts} ${traits} ${clades} \
        --auspice-config ${auspice_config} \
        --lat-longs ${lat_longs} \
        --output ncov.json
    """
}

process extractNewVariants {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(ncov_json) from extractvariants_in

  output:
  path("new_mutations.csv")

  script:
  """
  extract_new_variants.py --pipeline_dir . --out_path new_mutations.csv --new_sample_string Biohub
  """
}

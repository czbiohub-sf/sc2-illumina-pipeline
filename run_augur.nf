def helpMessage() {
  log.info"""
    Pipeline for running augur commands to build auspice visualization from main output.

    Usage:

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
      --sequences                   FASTA file of consensus sequences from main output
      --include_sequences           FASTA file of closest sequences from main output
      --metadata                    TSV of metadata from main output
      --ref_gb                      Reference Genbank file for augur
      --nextstrain_sequences        FASTA of sequences to build a tree with
      --minLength                   Minimum base pair length to allow assemblies to pass QC

    Nextstrain options:
      --nextstrain_ncov             Path to nextstrain/ncov directory (default: fetches from github)
      --sequences_per_group         Initial subsampling (default: 500)

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


sample_sequences  = file(params.sequences, checkIfExists: true)
included_fastas = file(params.include_sequences, checkIfExists: true)

process makeNextstrainInput {
    publishDir "${params.outdir}/nextstrain/data", mode: 'copy'
    stageInMode 'copy'

    input:
    path(sample_sequences)
    path(nextstrain_sequences)
    path(nextstrain_metadata_path)
    path(included_fastas)
    path(include_file)

    output:
    path('metadata.tsv') into (nextstrain_metadata, firstfilter_metadata, extractsamples_metadata, priorities_metadata, refinetree_metadata, infertraits_metadata, tipfreq_metadata, export_metadata)
    path('deduped_sequences.fasta') into (nextstrain_in, firstfilter_in)
    path('included_sequences.txt') into nextstrain_include
    path("internal_samples.txt") into sample_ids
    path("external_samples.txt") into external_ids

    script:
    currdate = new java.util.Date().format('yyyy-MM-dd')
    // Normalize the GISAID names using Nextstrain's bash script
    """
    make_nextstrain_input.py -ps ${nextstrain_sequences} -pm ${nextstrain_metadata_path} -ns ${sample_sequences} --date $currdate \
    -r 'North America' -c USA -div 'California' -loc 'San Francisco County' -sublab 'Biohub' \
    -subdate $currdate

    cat ${included_fastas} | grep '>' | awk -F '>' '{print \$2}' > included_samples.txt
    
    normalize_gisaid_fasta.sh all_sequences.fasta sequences.fasta
    cat included_samples.txt ${include_file} > included_sequences.txt
    cat ${included_fastas} >> sequences.fasta
    seqkit rmdup sequences.fasta > deduped_sequences.fasta

    cat ${sample_sequences} | grep '>' | awk -F '>' '{print \$2}' > internal_samples.txt
    cat deduped_sequences.fasta | grep '>' | awk -F '>' '{print \$2}' > external_samples.txt
    """

}



process firstFilter {
  label "nextstrain"

  input:
  path(sequences) from firstfilter_in
  path(metadata) from firstfilter_metadata
  path(exclude_file)
  path(include_file)

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
            --group-by division year month \
            --sequences-per-group ${params.sequences_per_group} \
            --output filtered_sequences.fasta
  """
}

process firstAlignment{
  label "process_large"
  label "nextstrain"

  input:
  path(sequences) from firstfiltered_ch
  path(ref_gb)

  output:
  path("aligned_sequences.fasta") into (firstaligned_ch, makepriorities_ch)

  script:
  """
  augur align \
            --sequences ${sequences} \
            --reference-sequence ${ref_gb} \
            --output aligned_sequences.fasta \
            --nthreads ${task.cpus} \
            --remove-reference \
            --fill-gaps
  """

}

process extractSampleSequences {
  label 'nextstrain'

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

  input:
  path(sample_sequences) from aligned_samples_ch
  path(aligned_sequences) from makepriorities_ch
  path(metadata) from priorities_metadata

  output:
  path("priorites.tsv") into priorities_tsv

  script:
  """
  priorities.py --alignment ${aligned_sequences} \
                --metadata ${metadata} \
                --focal-alignment ${sample_sequences} \
                --output priorities.tsv
  """
}

process filterStrains {
    label 'nextstrain'
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(sequences) from nextstrain_in
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
            --group-by division year month \
            --sequences-per-group 20 \
            --min-length ${params.minLength} \
            --output filtered.fasta \
    """
}

process buildTree {
    label 'nextstrain'
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
    label 'nextstrain'
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
    label 'nextstrain'
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
    label 'nextstrain'
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
    label 'nextstrain'
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
    label 'nextstrain'
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
    label 'nextstrain'
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
    label 'nextstrain'
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
    path('ncov.json')

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

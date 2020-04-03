def helpMessage() {
	log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --reads '*_R{1,2}_001.fastq.gz' --ref reference.fasta --primers primers.bed

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
      --reads                       Path to reads, must be in quotes
      --primers                     Path to BED file of primers
      --ref                         Path to FASTA reference sequence
      --gisaid_sequences            GISAID FASTA


    Options:
      --single_end [bool]           Specifies that the input is single-end reads
      --skip_trim_adapters [bool]   Skip trimming of illumina adapters. (NOTE: this does NOT skip the step for trimming spiked primers)
      --maxNs                       Max number of Ns to allow assemblies to pass QC
      --minLength                   Minimum base pair length to allow assemblies to pass QC
      --no_reads_quast              Run QUAST without aligning reads
      --gisaid_metadata             Metadata for GISAID sequences (default: fetches from github.com/nextstrain/ncov)
      --include_strains             File with included strains after augur filter (default: fetches from github.com/nextstrain/ncov)
      --exclude_strains             File with excluded strains for augur filter (default: fetches from github.com/nextstrain/ncov)
      --weights                     File with weights for augur traits (default: fetches from github.com/nextstrain/ncov)
      --clades                      File with clade for augur clades (default: fetches from github.com/nextstrain/ncov)
      --auspice_config              Config file for auspice (default: fetches from github.com/nextstrain/ncov)
      --lat_longs                   File with latitudes and longitudes for locations (default: fetches from github.com/nextstrain/ncov)

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

multiqc_config = file(params.multiqc_config, checkIfExists: true)

if (params.readPaths) {
    if (params.single_end){
        Channel
        .fromList(params.readPaths)
        .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true)] ] }
        .into {reads_ch; quast_reads}
    } else {
        Channel
        .fromList(params.readPaths)
        .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
        .into {reads_ch; quast_reads}
    }
} else {
    Channel
        .fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
        .into {reads_ch; quast_reads}
}

untrimmed_ch = params.skip_trim_adapters ? Channel.empty() : reads_ch

ref_fasta = file(params.ref, checkIfExists: true)
ref_gb = file(params.ref_gb, checkIfExists: true)
primer_bed = file(params.primers, checkIfExists: true)


process trimReads {
    tag { sampleName }

    cpus 2

    input:
    tuple(sampleName, file(reads)) from untrimmed_ch

    output:
    tuple(sampleName, file("*_val_*.fq.gz")) into trimmed_ch
    path("*") into trimmed_reports

    when:
    !params.skip_trim_adapters

    script:
    """
    trim_galore --fastqc --paired ${reads}
    """
}

unaligned_ch = params.skip_trim_adapters ? reads_ch : trimmed_ch

process alignReads {
    tag { sampleName }

    cpus 4

    input:
    tuple(sampleName, file(reads)) from unaligned_ch

    output:
    tuple(sampleName, file("${sampleName}.bam")) into aligned_reads

    script:
    """
    minimap2 -ax sr -R '@RG\\tID:${sampleName}\\tSM:${sampleName}' ${ref_fasta} ${reads} |
      samtools sort -@ 2 -O bam -o ${sampleName}.bam
    """
}

process trimPrimers {
	tag { sampleName }
	publishDir "${params.outdir}/samples/${sampleName}", mode: 'copy'

    input:
    tuple(sampleName, file(alignment)) from aligned_reads

    output:
    tuple(sampleName, file("${sampleName}.primertrimmed.bam")) into trimmed_bam_ch
    file("${sampleName}.primertrimmed.bam") into combined_vcf_bam

    script:
    """
    samtools view -F4 -q ${params.samQualThreshold} -o ivar.bam ${alignment}
    samtools index ivar.bam
    ivar trim -e -i ivar.bam -b ${primer_bed} -p ivar.out
    samtools sort -O bam -o ${sampleName}.primertrimmed.bam ivar.out.bam
    """
}

trimmed_bam_ch.into { quast_bam; consensus_bam; samtools_stats_bam;
                     sample_vcf_bam }

process samtoolsStats {
    tag { sampleName }
    input:
    tuple(sampleName, file(in_bam)) from samtools_stats_bam

    output:
    file("${sampleName}.samtools_stats") into samtools_stats_out


    script:
    """
    samtools stats ${in_bam} > ${sampleName}.samtools_stats
    """
}

process makeConsensus {
	tag { sampleName }
	publishDir "${params.outdir}/samples/${sampleName}", mode: 'copy'

	input:
	tuple(sampleName, path(bam)) from consensus_bam

	output:
	tuple(sampleName, path("${sampleName}.fa")) into quast_ch
        path("${sampleName}.fa") into merge_fastas_ch
        path("${sampleName}.stats.json") into stats_ch
        path("${sampleName}.depths.png")

	script:
	"""
        samtools index ${bam}
	samtools mpileup -A -d ${params.mpileupDepth} -Q0 ${bam} |
	  ivar consensus -q ${params.ivarQualThreshold} -t ${params.ivarFreqThreshold} -m ${params.minDepth} -n N -p ${sampleName}.primertrimmed.consensus
        clean_summarize_assembly.py \
             ${sampleName} \
             ${sampleName}.primertrimmed.bam \
             ${sampleName}.primertrimmed.consensus.fa \
             ${sampleName}
	"""
}

process sampleVariants {
    tag { sampleName }
    publishDir "${params.outdir}/samples/${sampleName}", mode: 'copy'

    input:
    tuple(sampleName, path(in_bam)) from sample_vcf_bam

    output:
    path("${sampleName}.vcf") into sample_vcf_out

    script:
    """
    bcftools mpileup -f ${ref_fasta} \
               ${in_bam} |
          bcftools call --ploidy 1 -m -P ${params.bcftoolsCallTheta} -v - |
          bcftools view -e 'DP<${params.minDepth}' > ${sampleName}.vcf
    """
}

process combinedVariants {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(in_bam) from combined_vcf_bam.collect()
    path(vcfs) from sample_vcf_out.collect()

    output:
    path("combined.vcf") into ch_vcf

    // use theta=.0003, assuming about 10 mutations between 2 random samples
    script:
    """
    printf "%s\\n" ${vcfs} | xargs -I % bgzip %
    printf "%s\\n" ${vcfs} | xargs -I % tabix %.gz
    printf "%s\\n" ${in_bam} | xargs -I % samtools index %
    bcftools merge \$(printf "%s.gz\n" ${vcfs}) | bcftools query -f '%CHROM\\t%POS\\t%END\\n' > variant_positions.txt
    bcftools mpileup -f ${ref_fasta} -R variant_positions.txt \
               ${in_bam} |
          bcftools call --ploidy 1 -m -P 0.0003 -v - > combined.vcf
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
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(in_json) from stats_ch.collect()

    output:
    path("combined.stats.tsv") into merged_stats_ch

    script:
    """merge_stats.py ${in_json} > combined.stats.tsv"""
}

process filterAssemblies {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(merged_stats) from merged_stats_ch
    path(merged_assemblies) from merged_assemblies_ch
    path(vcf) from ch_vcf

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

process quast {
    tag { sampleName }
    publishDir "${params.outdir}/QUAST", mode: 'copy'

    input:
    tuple(sampleName, path(assembly)) from quast_ch
    tuple(sampleName, path(bam)) from quast_bam
    tuple(sample, path(reads)) from quast_reads

    output:
    // Avoid name clash with other samples for MultiQC
    path("${sampleName}/*")
    path("${sampleName}/report.tsv") into multiqc_quast

    script:
    if (params.no_reads_quast)
    """
    quast --min-contig 0 -o ${sampleName} -r ${ref_fasta} -t ${task.cpus} --ref-bam ${bam} $assembly
    """
    else
    """
    quast --min-contig 0 -o ${sampleName} -r ${ref_fasta} -t ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} --ref-bam ${bam} $assembly
    """
}

// Set up GISAID files
gisaid_sequences = file(params.gisaid_sequences, checkIfExists: true)
gisaid_metadata = file(params.gisaid_metadata, checkIfExists: true)

process makeNextstrainInput {
    publishDir "${params.outdir}/nextstrain/data", mode: 'copy'
    stageInMode 'copy'

    input:
    path(sample_sequences) from nextstrain_ch
    path(gisaid_sequences)
    path(gisaid_metadata)

    output:
    path('metadata.tsv') into (nextstrain_metadata, refinetree_metadata, infertraits_metadata, tipfreq_metadata, export_metadata)
    path('sequences.fasta') into nextstrain_sequences

    script:
    currdate = new java.util.Date().format('yyyy-MM-dd')
    // Normalize the GISAID names using Nextstrain's bash script
    """
    make_nextstrain_input.py -ps ${gisaid_sequences} -pm ${gisaid_metadata} -ns ${sample_sequences} --date $currdate \
    -r 'North America' -c USA -div 'California' -loc 'San Francisco County' -origlab 'Biohub' -sublab 'Biohub' \
    -subdate $currdate

    normalize_gisaid_fasta.sh all_sequences.fasta sequences.fasta
    """

}

include_file = file(params.include_strains, checkIfExists: true)
exclude_file = file(params.exclude_strains, checkIfExists: true)

process filterStrains {
    label 'nextstrain'
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(sequences) from nextstrain_sequences
    path(metadata) from nextstrain_metadata
    path(include_file)
    path(exclude_file)

    output:
    path('filtered.fasta') into filtered_sequences_ch

    script:
    String exclude_where = "date='2020' date='2020-01-XX' date='2020-02-XX' date='2020-03-XX' date='2020-04-XX' date='2020-01' date='2020-02' date='2020-03' date='2020-04'"
    """
    augur filter \
            --sequences ${sequences} \
            --metadata ${metadata} \
            --include ${include_file} \
            --exclude ${exclude_file} \
            --exclude-where ${exclude_where} \
            --min-length 25000 \
            --group-by 'division year month' \
            --sequences-per-group 300 \
            --output filtered.fasta
    """

}

process alignSequences {
    label 'nextstrain'
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    cpus 4

    input:
    path(filtered_sequences) from filtered_sequences_ch
    path(ref_fasta)

    output:
    path('aligned.fasta') into (aligned_ch, refinetree_alignment, ancestralsequences_alignment)

    script:
    """
    augur align \
        --sequences ${filtered_sequences} \
        --reference-sequence ${ref_fasta} \
        --output aligned.fasta \
        --nthreads ${task.cpus} \
        --remove-reference \
        --fill-gaps
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

weights = file(params.weights, checkIfExists: true)

process inferTraits {
    label 'nextstrain'
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(tree) from infertraits_tree
    path(metadata) from infertraits_metadata
    path(weights)

    output:
    path('traits.json') into export_traits

    script:
    """
    augur traits \
        --tree ${tree} \
        --metadata ${metadata} \
        --weights ${weights} \
        --output traits.json \
        --columns country_exposure \
        --confidence \
        --sampling-bias-correction 2.5 \
    """
}

clades = file(params.clades, checkIfExists: true)

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

auspice_config = file(params.auspice_config, checkIfExists: true)
lat_longs = file(params.lat_longs, checkIfExists: true)

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

process multiqc {
	publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path(trim_galore_results) from trimmed_reports.collect().ifEmpty([])
    path("quast_results/*/*") from multiqc_quast.collect()
    path(samtools_stats) from samtools_stats_out.collect()
    path(multiqc_config)

    output:
    path("*multiqc_report.html")
    path("*_data")
    path("multiqc_plots")

	script:
	"""
	multiqc -f --config ${multiqc_config} ${trim_galore_results}  ${samtools_stats} quast_results/
	"""
}

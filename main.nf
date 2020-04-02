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


    Options:
      --single_end [bool]           Specifies that the input is single-end reads
      --skip_trim_adapters [bool]   Skip trimming of illumina adapters. (NOTE: this does NOT skip the step for trimming spiked primers)
      --maxNs                       Max number of Ns to allow assemblies to pass QC
      --minLength                   Minimum base pair length to allow assemblies to pass QC
      --no_reads_quast              Run QUAST without aligning reads

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
    file("${sampleName}.primertrimmed.bam") into vcf_in_ch

    script:
    """
    samtools view -F4 -q ${params.samQualThreshold} -o ivar.bam ${alignment}
    samtools index ivar.bam
    ivar trim -e -i ivar.bam -b ${primer_bed} -p ivar.out
    samtools sort -O bam -o ${sampleName}.primertrimmed.bam ivar.out.bam
    """
}

trimmed_bam_ch.into { quast_bam; consensus_bam; samtools_stats_in }

process samtoolsStats {
    tag { sampleName }
    input:
    tuple(sampleName, file(in_bam)) from samtools_stats_in

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
	  ivar consensus -q ${params.ivarQualThreshold} -t ${params.ivarFreqThreshold} -m ${params.ivarMinDepth} -n N -p ${sampleName}.primertrimmed.consensus
        clean_summarize_assembly.py \
             ${sampleName} \
             ${sampleName}.primertrimmed.bam \
             ${sampleName}.primertrimmed.consensus.fa \
             ${sampleName}
	"""
}

process callVariants {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(in_bam) from vcf_in_ch.collect()

    output:
    path("combined.vcf") into ch_vcf

    // use theta=.0003, assuming about 10 mutations between 2 random samples
    script:
    """
    bcftools mpileup -f ${ref_fasta} \
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
    path("filtered.fa")
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

process multiqc {
	publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path(trim_galore_results) from trimmed_reports.collect().ifEmpty([])
    path("quast_results/*/*") from multiqc_quast.collect()
    path(samtools_stats) from samtools_stats_out.collect()

    output:
    path("*multiqc_report.html")
    path("*_data")
    path("multiqc_plots")

	script:
	"""
	multiqc -f --config ${multiqc_config} ${trim_galore_results}  ${samtools_stats} quast_results/
	"""
}

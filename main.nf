def helpMessage() {
	log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --reads '*_R{1,2}_001.fastq.gz' --ref reference.fasta --primers primers.bed

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
      --reads                      Path to reads, must be in quotes
      --primers                     Path to BED file of primers
      --ref                         Path to FASTA reference sequence


    Options:
      --single_end [bool]             Specifies that the input is single-end reads
      --skip_trim_adaptors [bool]     Skip trimming of illumina adaptors. (NOTE: this does NOT skip the step for trimming spiked primers)

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

ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)

if (params.readPaths) {
    Channel
        .from(params.readPaths)
        .map { row -> [ row[0], row[1].each{file(it)} ] }
        .set {reads_ch}
} else {
    Channel
        .fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
        .set {reads_ch}
}

untrimmed_ch = params.skip_trim_adaptors ? Channel.empty() : reads_ch

ch_fasta = file(params.ref, checkIfExists: true)
ch_bed = file(params.primers, checkIfExists: true)

process trimReads {
    tag { sampleName }

    cpus 2

    input:
    tuple(sampleName, file(reads)) from untrimmed_ch

    output:
    tuple(sampleName, file("*_val_*.fq.gz")) into trimmed_ch
    path("*") into trimmed_reports

    script:
    """
    trim_galore --fastqc --paired ${reads}
    """
}

unaligned_ch = params.skip_trim_adaptors ? reads_ch : trimmed_ch

process alignReads {
    tag { sampleName }
    publishDir "${params.outdir}/${sampleName}"

    cpus 4

    input:
    tuple(sampleName, file(reads)) from unaligned_ch
    path(ref_fasta) from ch_fasta

    output:
    tuple(sampleName, file("${sampleName}.bam")) into aligned_reads

    script:
    """
    minimap2 -ax sr ${ref_fasta} ${reads} |
      samtools sort -@ 2 -O bam -o ${sampleName}.bam
    """
}

process trimPrimers {
	tag { sampleName }
	publishDir "${params.outdir}/${sampleName}"

    input:
    tuple(sampleName, file(alignment)) from aligned_reads
    path(primer_bed) from ch_bed

    output:
    tuple(sampleName, file("${sampleName}.primertrimmed.bam")) into trimmed_bam_ch

    script:
    """
    samtools view -F4 -o ivar.bam ${alignment}
    samtools index ivar.bam
    ivar trim -e -i ivar.bam -b ${primer_bed} -p ivar.out
    samtools sort -O bam -o ${sampleName}.primertrimmed.bam ivar.out.bam
    """
}

process makeConsensus {
	tag { sampleName }
	publishDir "${params.outdir}/${sampleName}"

	input:
	tuple(sampleName, path(bam)) from trimmed_bam_ch

	output:
	tuple(sampleName, path("${sampleName}.primertrimmed.consensus.fa"))

	script:
	"""
	samtools mpileup -A -d ${params.mpileupDepth} -Q0 ${bam} | 
	  ivar consensus -q ${params.ivarQualThreshold} -t ${params.ivarFreqThreshold} -m ${params.ivarMinDepth} -n N -p ${sampleName}.primertrimmed.consensus
	"""
}

process multiqc {
	publishDir "${params.outdir}/MultiQC"

	input:
	path(multiqc_config) from ch_multiqc_config
	path(trim_galore_results) from trimmed_reports.collect()

	output:
	path("*multiqc_report.html")
	path("*_data")
	path("multiqc_plots")

	script:
	"""
	multiqc -f --config ${multiqc_config} ${trim_galore_results}
	"""
}

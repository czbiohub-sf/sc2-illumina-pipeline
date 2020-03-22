Channel
    .fromFilePairs(params.fastqs)
    .set { reads_ch }

ch_fasta = file(params.ref, checkIfExists: true)
ch_bed = file(params.primers, checkIfExists: true)

process trimReads {
    tag { sampleName }

    cpus 2

    input:
    tuple(sampleName, file(reads)) from reads_ch

    output:
    tuple(sampleName, file("*_val_*.fq.gz")) into trimmed_ch

    script:
    """
    trim_galore --paired ${reads}
    """
}

process alignReads {
    tag { sampleName }
    publishDir "${params.outdir}/${sampleName}"

    cpus 4

    input:
    tuple(sampleName, file(reads)) from trimmed_ch
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

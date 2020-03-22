Channel
    .fromFilePairs(params.fastqs)
    .into { read_files_fastqc; read_files_minimap2 }

ch_fasta = file(params.ref, checkIfExists: true) 

process fastqc {
    conda "bioconda::fastqc"
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }

    cpus 2

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}

process minimap2 {
    conda "bioconda::minimap2 bioconda::samtools"
    publishDir "${params.outdir}", mode: 'copy'

    cpus 4

    input:
    tuple(val(name), file(reads)) from read_files_minimap2
    path(ref_fasta) from ch_fasta

    output:
    file "*.bam*" into minimap2_results

    script:
    """
    minimap2 -ax sr ${ref_fasta} ${reads} |
      samtools sort -@ 2 -O bam -o ${name}.bam
    samtools index ${name}.bam
    """
}

process samtools_flagstat {
    conda "bioconda::minimap2 bioconda::samtools"
    publishDir "${params.outdir}", mode: 'copy'

    cpus 1

    input:
    tuple(path(bam), path(bai)) from minimap2_results

    output:
    file "*.flagstat" into samtools_flagstat_results

    """
    samtools flagstat ${bam} > ${bam}.flagstat
    """
}

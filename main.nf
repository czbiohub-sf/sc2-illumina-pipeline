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
      --gisaid_metadata             Path to GISAID metadata from Nextstrain
      --gisaid_sequences            Path to GISAID sequences 
      --collection_date             Collection date (e.g. '2020-03-08')


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

ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)

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

ch_fasta = file(params.ref, checkIfExists: true)
ch_bed = file(params.primers, checkIfExists: true)

// Set up files for QUAST
quast_ref = file(params.ref, checkIfExists: true)

// GISAID files
gisaid_sequences = file(params.gisaid_sequences, checkIfExists: true)
gisaid_metadata = file(params.gisaid_metadata, checkIfExists: true)

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
    publishDir "${params.outdir}/samples/${sampleName}", mode: 'copy'

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
	publishDir "${params.outdir}/samples/${sampleName}", mode: 'copy'

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
	publishDir "${params.outdir}/samples/${sampleName}", mode: 'copy'

	input:
	tuple(sampleName, path(bam)) from trimmed_bam_ch

	output:
	tuple(sampleName, path("${sampleName}.primertrimmed.consensus.fa")) into quast_ch

	script:
	"""
	samtools mpileup -A -d ${params.mpileupDepth} -Q0 ${bam} | 
	  ivar consensus -q ${params.ivarQualThreshold} -t ${params.ivarFreqThreshold} -m ${params.ivarMinDepth} -n N -p ${sampleName}.primertrimmed.consensus
	"""
}

process quast {
	tag { sampleName }
	publishDir "${params.outdir}/QUAST", mode: 'copy'

	input:
	tuple(sampleName, path(assembly)) from quast_ch
	path(quast_ref)
	tuple(sample, path(reads)) from quast_reads

	output:
    // Avoid name clash with other samples for MultiQC
    path("${sampleName}/*")
    path("${sampleName}/report.tsv") into multiqc_quast
    tuple(sampleName, path(assembly), path("${sampleName}/report.tsv")) into sort_assemblies_ch

    script:
    if (params.no_reads_quast)
    """
    quast --min-contig 0 -o ${sampleName} -r ${quast_ref} -t ${task.cpus} $assembly
    """
    else
    """
    quast --min-contig 0 -o ${sampleName} -r ${quast_ref} -t ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} $assembly
    """
}

process sortAssemblies {
    tag {sampleName}

    input:
    tuple(sampleName, path(assembly), path(report_tsv)) from sort_assemblies_ch

    output:
    path("passed_QC/*") into passed_qc_asm

    script:
    // create placeholder file so nextflow always has an output
    """
    mkdir -p passed_QC
    touch passed_QC/${sampleName}_foo.txt
    mkdir -p failed_QC
    parse_quast_unaligned_report.py --report ${report_tsv} --assembly ${assembly} -n ${params.maxNs} -l ${params.minLength}
    """
}

process combineFiles {
    publishDir "${params.outdir}/submission_files", mode: 'copy'
    
    input:
    path(asm_files) from passed_qc_asm.collect()

    output:
    path("combined_sequences.fasta") into nextstrain_ch
    path("all_sequences/*")

    script:
    """
    mkdir -p all_sequences
    touch combined_sequences.fasta
    for f in ${asm_files}
    do
    if [[ \$f != *foo.txt* ]]
    then
    cat \$f >> combined_sequences.fasta
    cp \$f all_sequences/
    fi
    done
    """
}

process makeNextstrainInput {
    publishDir "${params.outdir}/nextstrain/data", mode: 'copy'
    stageInMode 'copy'

    input:
    path(sample_sequences) from nextstrain_ch
    path(gisaid_sequences)
    path(gisaid_metadata)

    output:
    path('metadata.tsv') into nextstrain_metadata
    path('norm_sequences.fasta') into nextstrain_sequences

    script:
    currdate = new java.util.Date().format('yyyy-MM-dd')
    // Need to normalize the GISAID fasta here to avoid awk error in runNextstrain
    """
    make_nextstrain_input.py -ps ${gisaid_sequences} -pm ${gisaid_metadata} -ns ${sample_sequences} --date ${params.collection_date} \
    -r 'North America' -c USA -div 'California' -loc 'San Francisco County' -origlab 'Biohub' -sublab 'Biohub' \
    -subdate $currdate

    normalize_gisaid_fasta.sh sequences.fasta norm_sequences.fasta
    """

}

// process fetchNextstrain {
//     publishDir "${params.outdir}/nextstrain"
//     echo true

//     input:
//     path(metadata) from nextstrain_metadata
//     path(sequences) from nextstrain_sequences

//     output:


//     script:
//     """
//     wget https://github.com/nextstrain/ncov/archive/master.zip
//     unzip master.zip
//     cp ${sequences} ncov-master/data/sequences.fasta
//     cp ${metadata} ncov-master/data/
//     cd ncov-master
//     snakemake -p
//     """
// }

process multiqc {
	publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path(multiqc_config) from ch_multiqc_config
    path(trim_galore_results) from trimmed_reports.collect().ifEmpty([])
    path("quast_results/*/*") from multiqc_quast.collect()

    output:
    path("*multiqc_report.html")
    path("*_data")
    path("multiqc_plots")

	script:
	"""
	multiqc -f --config ${multiqc_config} ${trim_galore_results}  quast_results/
	"""
}

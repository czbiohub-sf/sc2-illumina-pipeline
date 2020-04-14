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
      --ref_gb                      Reference Genbank file for augur
      --blast_sequences             FASTA of sequences to BLAST for nearest neighbors
      --nextstrain_sequences        FASTA of sequences to build a tree with


    Consensus calling options:
      --kraken2_db                  Path to kraken db (default: "")
      --exclude_samples             comma-separated string of samples to exclude from analysis
      --single_end [bool]           Specifies that the input is single-end reads
      --skip_trim_adapters [bool]   Skip trimming of illumina adapters. (NOTE: this does NOT skip the step for trimming spiked primers)
      --maxNs                       Max number of Ns to allow assemblies to pass QC
      --minLength                   Minimum base pair length to allow assemblies to pass QC
      --no_reads_quast              Run QUAST without aligning reads
      --clades                      TSV file with columns clade, gene, site, alt (augur clades format)
      --qpcr_primers                BED file with positions of qPCR primers to check for variants


    Nextstrain options:
      --nextstrain_ncov             Path to nextstrain/ncov directory (default: fetches from github)
      --subsample_gisaid            Subsample GISAID for alignment


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

ref_fasta = file(params.ref, checkIfExists: true)
ref_gb = file(params.ref_gb, checkIfExists: true)
primer_bed = file(params.primers, checkIfExists: true)

if (params.readPaths) {
    if (params.single_end){
        Channel
        .fromList(params.readPaths)
        .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true)] ] }
        .set { reads_ch }
    } else {
        Channel
        .fromList(params.readPaths)
        .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
        .set { reads_ch }
    }
} else {
    Channel
        .fromFilePairs(params.reads, size: params.single_end ? 1 : 2)
        .set { reads_ch }
}

// remove excluded samples from reads_ch
exclude_samples = params.exclude_samples.split(",")
reads_ch = reads_ch.filter { !exclude_samples.contains(it[0]) }

reads_ch.into { unaligned_reads; stats_reads; }
reads_ch = unaligned_reads

if (params.kraken2_db == "") {
    // skip kraken
    kraken2_reads_in = Channel.empty()
    kraken2_db = Channel.empty()
} else {
    // send reads to kraken, and empty the reads channel
    kraken2_reads_in = reads_ch
    reads_ch = Channel.empty()
    kraken2_db = file(params.kraken2_db, checkIfExists: true)
}

process filterReads {
    tag { sampleName }
    label 'process_large'

    input:
    path(db) from kraken2_db
    path(ref_fasta)
    tuple(sampleName, file(reads)) from kraken2_reads_in

    output:
    tuple(sampleName, file("${sampleName}_covid_*.fq.gz")) into kraken2_reads_out

    script:
    """
    minimap2 -ax sr ${ref_fasta} ${reads} |
      samtools sort -@ ${task.cpus-1} -n -O bam -o mapped.bam
    samtools fastq -@ ${task.cpus-1} -G 12 -1 paired1.fq.gz -2 paired2.fq.gz \
       -0 /dev/null -s /dev/null -n -c 6 \
       mapped.bam
    rm mapped.bam

    LINES=\$(zcat paired1.fq.gz | wc -l)
    if [ "\$LINES" -gt 0 ];
    then
        kraken2 --db ${db} \
          --report ${sampleName}.kraken2_report \
          --classified-out "${sampleName}_classified#.fq" \
          --output - \
          --memory-mapping --gzip-compressed --paired \
          paired1.fq.gz paired2.fq.gz

        rm paired1.fq.gz paired2.fq.gz

        grep --no-group-separator -A3 "kraken:taxid|2697049" \
             ${sampleName}_classified_1.fq \
             > ${sampleName}_covid_1.fq || [[ \$? == 1 ]]

        grep --no-group-separator -A3 "kraken:taxid|2697049" \
             ${sampleName}_classified_2.fq \
             > ${sampleName}_covid_2.fq || [[ \$? == 1 ]]

        gzip ${sampleName}_covid_1.fq
        gzip ${sampleName}_covid_2.fq

        rm ${sampleName}_classified_*.fq
    else
        mv paired1.fq.gz ${sampleName}_covid_1.fq.gz
        mv paired2.fq.gz ${sampleName}_covid_2.fq.gz
    fi
    """
}

//send kraken output back to the reads channel
reads_ch = reads_ch.concat(kraken2_reads_out)

if (params.skip_trim_adapters) {
    // skip trimming
    trimgalore_reads_in = Channel.empty()
} else {
    // send reads to trim_galore, and empty the reads channel
    trimgalore_reads_in = reads_ch
    reads_ch = Channel.empty()
}

process trimReads {
    tag { sampleName }
    publishDir "${params.outdir}/trimmed-reads", mode: 'copy',
        saveAs: { x -> x.endsWith(".fq.gz") ? x : null }

    cpus 2

    input:
    tuple(sampleName, file(reads)) from trimgalore_reads_in

    output:
    tuple(sampleName, file("*_val_*.fq.gz")) into trimgalore_reads_out
    path("*") into trimmed_reports

    script:
    """
    LINES=\$(zcat ${reads[0]} | wc -l)
    if [ "\$LINES" -gt 0 ];
    then
        trim_galore --fastqc --paired ${reads}
        TRIMMED=\$(zcat ${sampleName}_covid_1_val_1.fq.gz | wc -l)
        if [ "\$TRIMMED" == 0 ];
        then
            rm -r *fastqc.zip
        fi
    else
        cp ${reads[0]} ${sampleName}_covid_1_val_1.fq.gz
        cp ${reads[1]} ${sampleName}_covid_2_val_2.fq.gz
    fi
    """
}

// send trim_galore output back to the reads channel
reads_ch = reads_ch.concat(trimgalore_reads_out)

// send reads to minimap2 and quast
reads_ch.into { minimap2_reads_in; quast_reads }

process alignReads {
    tag { sampleName }

    cpus 4

    input:
    tuple(sampleName, file(reads)) from minimap2_reads_in
    path(ref_fasta)

    output:
    tuple(sampleName, file("${sampleName}.bam")) into bam2trimPrimers

    script:
    """
    minimap2 -ax sr -R '@RG\\tID:${sampleName}\\tSM:${sampleName}' ${ref_fasta} ${reads} |
      samtools sort -@ 2 -O bam -o ${sampleName}.bam
    """
}

process trimPrimers {
    tag { sampleName }
    publishDir "${params.outdir}/aligned-reads", mode: 'copy'

    input:
    tuple(sampleName, file(alignment)) from bam2trimPrimers
    path(primer_bed)

    output:
    tuple(sampleName, file("${sampleName}.primertrimmed.bam")) into trimmed_bam_ch;
    file("${sampleName}.primertrimmed.bam.bai")

    script:
    """
    samtools view -F4 -q ${params.samQualThreshold} -o ivar.bam ${alignment}
    samtools index ivar.bam
    ivar trim -e -i ivar.bam -b ${primer_bed} -p ivar.out
    samtools sort -O bam -o ${sampleName}.primertrimmed.bam ivar.out.bam
    samtools index ${sampleName}.primertrimmed.bam
    """
}

trimmed_bam_ch.into { quast_bam; consensus_bam; stats_bam; intrahost_bam }

process makeConsensus {
	tag { sampleName }
	publishDir "${params.outdir}/consensus-seqs", mode: 'copy'

	input:
	tuple(sampleName, path(bam)) from consensus_bam

	output:
	tuple(sampleName, path("${sampleName}.consensus.fa")) into (consensus_fa, quast_ch, compareconsensus_ch, consensustoneighbor_ch)

	script:
	"""
        samtools index ${bam}
	samtools mpileup -A -d ${params.mpileupDepth} -Q0 ${bam} |
	  ivar consensus -q ${params.ivarQualThreshold} -t ${params.ivarFreqThreshold} -m ${params.minDepth} -n N -p ${sampleName}.primertrimmed.consensus
        echo '>${sampleName}' > ${sampleName}.consensus.fa
        seqtk seq -l 50 ${sampleName}.primertrimmed.consensus.fa | tail -n +2 >> ${sampleName}.consensus.fa
	"""
}

blast_sequences = file(params.blast_sequences, checkIfExists: true)

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
    tuple(sampleName, path(assembly)) from compareconsensus_ch
    path(dbsequences) from blast_clean_ch
    path(blastdb) from blastdb_ch
    path(ref_fasta)

    output:
    path("${sampleName}.blast.tsv")
    tuple(sampleName, path("${sampleName}_nearest_blast.fasta")) into (nearest_neighbor, collectnearest_ch)

    script:
    """
    get_top_hit.py --minLength ${params.minLength} --sequences ${dbsequences} --sampleName ${sampleName} --assembly ${assembly} --default ${ref_fasta}
    """
}

process collectNearest {

    input:
    path(fastas) from collectnearest_ch.map{it[1]}.collect()

    output:
    path("included_samples.txt") into included_samples_ch
    path("included_samples.fasta") into included_fastas_ch

    script:
    if (params.subsample_gisaid)
    """
    cat ${fastas} > included_samples.fasta
    cat ${fastas} | grep '>' | awk -F '>' '{print \$2}' > included_samples.txt
    """
    else
    """
    touch included_samples.fasta
    touch included_samples.txt
    """
}

nearest_neighbor
    .join(consensustoneighbor_ch)
    .set{nearest_blast_ch}

process nearestVariants {
    tag {sampleName}
    publishDir "${params.outdir}/sample-variants", mode: 'copy'

    input:
    tuple(sampleName, path(nearest_blast), path(assembly)) from nearest_blast_ch

    output:
    path("${sampleName}.nearest_realigned.bcftools_stats") into nearest_realigned_stats
    tuple(sampleName, path("${sampleName}.nearest_realigned.vcf")) into nearest_realigned_vcf


    script:
    """
    minimap2 -ax asm5 -R '@RG\\tID:${sampleName}\\tSM:${sampleName}' \
      ${nearest_blast} ${assembly} |
      samtools sort -O bam -o ${sampleName}.nearest_realigned.bam
    samtools index ${sampleName}.nearest_realigned.bam
    bcftools mpileup -f ${nearest_blast}  ${sampleName}.nearest_realigned.bam |
      bcftools call --ploidy 1 -m -P ${params.bcftoolsCallTheta} -v - \
      > ${sampleName}.nearest_realigned.vcf
    bcftools stats ${sampleName}.nearest_realigned.vcf > ${sampleName}.nearest_realigned.bcftools_stats
    """
}

process quast {
   tag { sampleName }
   publishDir "${params.outdir}/QUAST", mode: 'copy'

   input:
   tuple(sampleName, path(assembly)) from quast_ch
   tuple(sampleName, path(bam)) from quast_bam
   tuple(sample, path(reads)) from quast_reads
   path(ref_fasta)

   output:
   // Avoid name clash with other samples for MultiQC
   path("${sampleName}/*") into multiqc_quast

   script:
   if (params.no_reads_quast)
   """
   run_quast.py --noreads --assembly ${assembly} --sample ${sampleName} --ref ${ref_fasta} \
    --threads ${task.cpus} --bam ${bam}
   """
   else
   """
   run_quast.py --noreads --assembly ${assembly} --sample ${sampleName} --ref ${ref_fasta} \
    --threads ${task.cpus} --bam ${bam} --R1 ${reads[0]} --R2 ${reads[1]}
   """
}

consensus_fa.into { quast_ch; stats_fa; merge_fastas_ch; realign_fa }
merge_fastas_ch = merge_fastas_ch.map { it[1] }

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

process callVariants {
    tag { sampleName }
    publishDir "${params.outdir}/sample-variants", mode: 'copy'

    input:
    tuple(sampleName, path(in_bam)) from call_variants_bam
    path(ref_fasta)

    output:
    tuple(sampleName, path("${sampleName}.vcf")) into (primer_variants_ch, assignclades_in)
    path("${sampleName}.bcftools_stats") into bcftools_stats_ch

    script:
    """
    bcftools mpileup -f ${ref_fasta} ${in_bam} |
      bcftools call --ploidy 1 -m -P ${params.bcftoolsCallTheta} -v - \
      > ${sampleName}.vcf
    bcftools stats ${sampleName}.vcf > ${sampleName}.bcftools_stats
    """
}

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

stats_reads
    .join(stats_bam)
    .join(stats_fa)
    .join(primer_variants_vcf)
    .join(assignclades_out)
    .join(nearest_realigned_vcf)
    .set { stats_ch_in }

process computeStats {
    tag { sampleName }
    publishDir "${params.outdir}/coverage-plots", mode: 'copy',
        saveAs: { x -> x.endsWith(".png") ? x : null }

    input:
    tuple(sampleName,
          file(reads),
          file(trimmed_filtered_bam),
          file(in_fa),
          file(primer_vcf),
          file(in_clades),
          file(neighbor_vcf)) from stats_ch_in

    output:
    file("${sampleName}.samtools_stats") into samtools_stats_out
    path("${sampleName}.stats.json") into stats_ch
    path("${sampleName}.depths.png")

    script:
    """
    samtools index ${trimmed_filtered_bam}
    samtools stats ${trimmed_filtered_bam} > ${sampleName}.samtools_stats
    alignment_assembly_stats.py \
        --sample_name ${sampleName} \
        --cleaned_bam ${trimmed_filtered_bam} \
        --samtools_stats ${sampleName}.samtools_stats \
        --assembly ${in_fa} \
        --primervcf ${primer_vcf} \
        --neighborvcf ${neighbor_vcf} \
        --clades ${in_clades} \
        --out_prefix ${sampleName} \
        --reads ${reads}
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

if (!params.intrahost_variants) {
    intrahost_bam = Channel.empty()
} else {
    intrahost_bam = intrahost_bam.map { it[1] }.collect()
}

process intrahostVariants {
    publishDir "${params.outdir}/",
        mode: 'copy'

    cpus params.intrahost_variants_cpu

    input:
    path(bam) from intrahost_bam
    path(ref_fasta)

    output:
    path("intrahost-variants." +
         "ploidy${params.intrahost_ploidy}-" +
         "minfrac${params.intrahost_min_frac}.vcf")

    script:
    """
    ls ${bam} | xargs -I % samtools index %
    samtools faidx ${ref_fasta}
    freebayes-parallel <(fasta_generate_regions.py ${ref_fasta}.fai 1000) ${task.cpus} \
        --ploidy ${params.intrahost_ploidy} \
        --min-alternate-fraction ${params.intrahost_min_frac} \
        -f ${ref_fasta} ${bam} |
        bcftools view -g het \
        > intrahost-variants.ploidy${params.intrahost_ploidy}-minfrac${params.intrahost_min_frac}.vcf
    """
}

// Set up GISAID files

if (params.nextstrain_ncov != "" && params.nextstrain_sequences != "") {
    nextstrain_sequences_ch = Channel.from(file(params.nextstrain_sequences, checkIfExists: true))

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
} else {
    nextstrain_sequences_ch = Channel.empty()
    nextstrain_metadata_path = Channel.empty()
    nextstrain_config = Channel.empty()
    include_file = Channel.empty()
    exclude_file = Channel.empty()
    clades = Channel.empty()
    auspice_config = Channel.empty()
    lat_longs = Channel.empty()
}

process makeNextstrainInput {
    publishDir "${params.outdir}/nextstrain/data", mode: 'copy'
    stageInMode 'copy'

    input:
    path(sample_sequences) from nextstrain_ch
    path(nextstrain_sequences) from nextstrain_sequences_ch
    path(nextstrain_metadata_path)

    output:
    path('metadata.tsv') into (nextstrain_metadata, refinetree_metadata, infertraits_metadata, tipfreq_metadata, export_metadata)
    path('sequences.fasta') into nextstrain_sequences

    script:
    currdate = new java.util.Date().format('yyyy-MM-dd')
    // Normalize the GISAID names using Nextstrain's bash script
    """
    make_nextstrain_input.py -ps ${nextstrain_sequences} -pm ${nextstrain_metadata_path} -ns ${sample_sequences} --date $currdate \
    -r 'North America' -c USA -div 'California' -loc 'San Francisco County' -origlab 'Biohub' -sublab 'Biohub' \
    -subdate $currdate

    normalize_gisaid_fasta.sh all_sequences.fasta sequences.fasta
    """

}

process filterStrains {
    label 'nextstrain'
    publishDir "${params.outdir}/nextstrain/results", mode: 'copy'

    input:
    path(sequences) from nextstrain_sequences
    path(metadata) from nextstrain_metadata
    path(include_file)
    path(exclude_file)
    path(included_samples) from included_samples_ch
    path(included_fastas) from included_fastas_ch

    output:
    path('filtered.fasta') into filtered_sequences_ch

    script:
    String exclude_where = "date='2020' date='2020-01-XX' date='2020-02-XX' date='2020-03-XX' date='2020-04-XX' date='2020-01' date='2020-02' date='2020-03' date='2020-04'"
    if (params.subsample_gisaid)
    """
    cat ${included_samples} >> ${include_file}
    cat ${included_fastas} >> ${sequences}
    seqkit rmdup ${sequences} > deduped_sequences.fasta
    augur filter \
            --sequences deduped_sequences.fasta \
            --metadata ${metadata} \
            --include ${include_file} \
            --exclude ${exclude_file} \
            --exclude-where ${exclude_where} \
            --min-length ${params.minLength} \
            --group-by division year month \
            --sequences-per-group 1 \
            --output filtered.fasta
    """

    else
    """
    augur filter \
            --sequences ${sequences} \
            --metadata ${metadata} \
            --include ${include_file} \
            --exclude ${exclude_file} \
            --exclude-where ${exclude_where} \
            --min-length ${params.minLength} \
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

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path(trim_galore_results) from trimmed_reports.collect().ifEmpty([])
    path("quast_results/*/*") from multiqc_quast.collect()
    path(samtools_stats) from samtools_stats_out.collect()
    path(multiqc_config)
    path(bcftools_stats) from bcftools_stats_ch.collect().ifEmpty([])
    path(primer_stats) from primer_stats_ch.collect().ifEmpty([])
    path(nearest_realigned) from nearest_realigned_stats.collect().ifEmpty([])

    output:
    path("*multiqc_report.html")
    path("*_data")
    path("multiqc_plots")

    // TODO: add trim_galore results (currently breaking for empty fastqs?)
    script:
    """
    multiqc -f -ip --config ${multiqc_config} \
        ${samtools_stats} quast_results/ ${bcftools_stats} \
        ${primer_stats} ${nearest_realigned} ${trim_galore_results}
    """
}

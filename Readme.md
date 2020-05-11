# SARS-CoV-2 Genomic Epidemiology Pipeline

This pipeline has two major components:

1. Assembly and QC of SARs-CoV-2 genomes from raw sequencing reads (fastq files).
2. Constructing a phylogenetic tree containing high-quality genomes from the pipeline and contextual reference genomes from GISAID.

For assembly, we use a simple approach based on alignment to reference.

For phylogenetics, we use the [augur](https://nextstrain.org/docs/bioinformatics/introduction-to-augur) pipeline developed by the Nextstrain team.

# Typical usage

For generating consensus genomes from reads:

```{sh}
nextflow run call_consensus.nf -profile docker \
    --reads 's3://path/to/reads/*_R{1,2}_001.fastq.gz*' \
    --outdir 's3://path/to/outdir'
```

## Outputs

Below is a non-exhaustive list of the outputs from the `call_consensus.nf` pipeline.

```
├── combined.fa
├── filtered.fa
├── combined.vcf
├── call_consensus-stats
│   ├── combined.stats.tsv
│   └── filtered.stats.tsv
├── MultiQC
│   └── multiqc_report.html
├── aligned-reads
│   ├── sample1.primertrimmed.bam
│   ├── sample1.primertrimmed.bam.bai
│   ├── sample2.primertrimmed.bam
│   └── ...
├── coverage-plots
│   ├── sample1.depths.png
│   ├── sample2.depths.bam
│   └── ...
├── trimmed-reads
│   ├── sample1_1_val_1.fq.gz
│   ├── sample1_2_val_2.fq.gz
│   ├── sample2_1_val_1.fq.gz
│   ├── sample2_2_val_2.fq.gz
│   └── ...
```

- `combined.fa`: Fasta of consensus sequences for all samples.
- `filtered.fa`: Fasta of consensus sequences for samples passing
  filters (default: require 25000 nonmissing genome).
- `combined.vcf`: VCF of variants for all samples.
- `call_consensus-stats`: Folder containing basic QC stats for consensus genomes.
- `MultiQC/multiqc_report.html`: QC report from MultiQC (aggregates stats from quast,
  samtools, trim-galore, bcftools).
- `aligned-reads`:

# Testing

Simple test to make sure things aren't broken:

```{sh}
nextflow run call_consensus.nf --profile docker,test
```

# Benchmarking

Simple benchmark (for mapping, not speed). Run after algorithm changes to see how accuracy might be affected. Result in `benchmark/call_consensus-stats/combined.stats.tsv`

```{sh}
nextflow run call_consensus.nf --profile docker,benchmark
```

# Analysis pipeline

`run_analysis.nf` runs the nextstrain augur/auspice pipeline, as well as several other ad hoc analyses. (To be refactored soon).

# Pipeline Overview

TODO

# Sequencing Protocols

The current pipeline is designed to work with the MSSPE spiked primer enrichment [protocol](https://www.protocols.io/private/32717E8D59E211EABDB40242AC110003?step=4). Future versions will adapt to amplicon-based short- and long- read sequencing.

# Acknowledgments

Initial version of this pipeline was based on
https://github.com/connor-lab/ncov2019-artic-nf

Augur/Auspice pipeline based on
https://github.com/nextstrain/ncov

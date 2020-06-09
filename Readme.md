# SARS-CoV-2 Consensus Genome Pipeline

This pipeline generates consensus SARS-CoV-2 genomes from fastq
files. We are using it on the following types of sequencing data:
1. Metagenomic sequencing enriched for SARS-CoV-2 reads
   ([protocols.io](https://www.protocols.io/private/32717E8D59E211EABDB40242AC110003?step=4)).
2. Amplicon-based short-read sequencing (using ARTIC v3 protocol).

# Typical usage

For generating consensus genomes from reads:

```{sh}
nextflow run main.nf -profile docker \
    --reads '[s3://]path/to/reads/*_R{1,2}_001.fastq.gz*' \
    --kraken2_db '[s3://]path/to/kraken2db' \
    --outdir '[s3://]path/to/outdir'
```

The kraken2db can be downloaded from https://genexa.ch/sars2-bioinformatics-resources/.

## Outputs

Below is a non-exhaustive list of the outputs from the `main.nf` pipeline.

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
- `aligned-reads`: Bamfile of trimmed reads mapping to reference.

# Testing

Simple test to make sure things aren't broken:

```{sh}
nextflow run main.nf -profile docker,test
```

# Benchmarking

Simple benchmark (for mapping, not speed). Run after algorithm changes to see how accuracy might be affected. Result in `benchmark/call_consensus-stats/combined.stats.tsv`

```{sh}
nextflow run main.nf --profile docker,benchmark
```

# Analysis pipeline

Analysis pipeline has moved here: https://github.com/czbiohub/sc2-ngs-analysis

# Pipeline Overview

TODO

# Acknowledgments

Initial version of this pipeline was based on
https://github.com/connor-lab/ncov2019-artic-nf

# SARS-CoV-2 Consensus Genome Pipeline

This pipeline generates consensus SARS-CoV-2 genomes from fastq
files. We are using it on the following types of sequencing data:
1. Metagenomic sequencing enriched for SARS-CoV-2 reads
   ([protocols.io](https://www.protocols.io/private/32717E8D59E211EABDB40242AC110003?step=4)).
2. Amplicon-based short-read sequencing (using ARTIC v3 protocol).

<!-- MarkdownTOC -->

- [Typical usage](#typical-usage)
- [Testing](#testing)
- [Benchmarking](#benchmarking)
- [Documentation](#documentation)
- [Acknowledgments](#acknowledgments)

<!-- /MarkdownTOC -->

# Typical usage

For generating consensus genomes from reads:

```{sh}
nextflow run czbiohub/sc2-msspe-bioinfo -profile docker \
    --reads '[s3://]path/to/reads/*_R{1,2}_001.fastq.gz*' \
    --kraken2_db '[s3://]path/to/kraken2db' \
    --outdir '[s3://]path/to/outdir'
```

The kraken2db can be downloaded from https://genexa.ch/sars2-bioinformatics-resources/.

# Testing

Simple test to make sure things aren't broken:

```{sh}
nextflow run czbiohub/sc2-msspe-bioinfo -profile docker,test
```

# Benchmarking

Simple benchmark (for mapping, not speed). Run after algorithm changes to see how accuracy might be affected. Result in `benchmark/call_consensus-stats/combined.stats.tsv`

```{sh}
nextflow run czbiohub/sc2-msspe-bioinfo --profile docker,benchmark
```

# Documentation

The czbiohub/sc2-msspe-bioinfo pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/running.md)
3. [Pipeline overview](docs/overview.md)
4. [Output](docs/output.md)


# Acknowledgments

Initial version of this pipeline was based on
https://github.com/connor-lab/ncov2019-artic-nf

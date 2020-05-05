# SARS-CoV-2 Genomic Epidemiology Pipeline

This pipeline has two major components:

1. Assembly and QC of SARs-CoV-2 genomes from raw sequencing reads (fastq files).
2. Constructing a phylogenetic tree containing high-quality genomes from the pipeline and contextual reference genomes from GISAID.

For assembly, we use a simple approach based on alignment to reference.

For phylogenetics, we use the [augur](https://nextstrain.org/docs/bioinformatics/introduction-to-augur) pipeline developed by the Nextstrain team.

# Sequencing Protocols

The current pipeline is designed to work with the MSSPE spiked primer enrichment [protocol](https://www.protocols.io/private/32717E8D59E211EABDB40242AC110003?step=4). Future versions will adapt to amplicon-based short- and long- read sequencing.

# Acknowledgments

Initial version of this pipeline was based on
https://github.com/connor-lab/ncov2019-artic-nf

Augur/Auspice pipeline based on
https://github.com/nextstrain/ncov

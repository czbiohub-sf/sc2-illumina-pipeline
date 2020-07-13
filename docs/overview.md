# Pipeline overview
<!-- MarkdownTOC -->

- [Processes](#processes)
	- [`filterRefReads`](#filterrefreads)
	- [`quantifyERCCs`](#quantifyerccs)
	- [`filterReads`](#filterreads)
	- [`trimReads`](#trimreads)
	- [`alignReads`](#alignreads)
	- [`trimPrimers`](#trimprimers)
	- [`makeConsensus`](#makeconsensus)
	- [`quast`](#quast)
	- [`callVariants`](#callvariants)
	- [`computeStats`](#computestats)
	- [`combinedVariants`](#combinedvariants)
	- [`mergeAllAssemblies`](#mergeallassemblies)
	- [`mergeAssemblyStats`](#mergeassemblystats)
	- [`filterAssemblies`](#filterassemblies)
	- [`multiqc`](#multiqc)

<!-- /MarkdownTOC -->
## Processes

### `filterRefReads`

Filters host reads by mapping with `minimap2` to the file given with `--ref_host` and keeping unmapped reads. This process can be skipped with the flag `--skip_filter_ref`. The filtered reads will be output into the folder `filtered-reads`. These reads are send to [`filterReads`](#filterreads).

### `quantifyERCCs`

Quantify ERCC reads by mapping with `minimap2` to the file given with `--ercc_fasta`. The individual output files are saved to the folder `ercc-stats` and the information is summarized in `combined.stats.tsv`.

### `filterReads`

Filter for only SARS-CoV-2 reads. This is accomplished by first mapping with `minimap2` to the file given with `--ref` and only keeping mapped reads. These reads are then run through `kraken2` against the database provided by `--kraken2_db`. Only reads that are assigned uniquely to SARS-CoV-2 are kept. These reads are sent to [`trimReads`](#trimreads).

### `trimReads`

Trim adapter sequences and perform quality-trimming using [Trim Galore](https://github.com/FelixKrueger/TrimGalore). Skip this step with the flag `--skip_trim_adapters`. The trimmed reads will be published to the folder `trimmed-reads`. These will be sent to [`alignReads`](#alignreads) and [`callVariants`](#callvariants).

### `alignReads`

Align the reads to the reference genome provided by `--ref` with `minimap2` and sort by position with `samtools sort`. The BAM files are sent to [`trimPrimers`](#trimprimers)

### `trimPrimers`

Trim the primers provided in `--primers` from the BAM using `ivar trim`. The primers will be soft-clipped from the ends of the reads in the BAM file. This step is necessary to properly call consensus alleles in these regions as the primer sequences can be artificially amplified. The trimmed BAMs will be published to the folder `aligned-reads`.

### `makeConsensus`

Call a consensus genome using `ivar consensus`. See [here](running.md#primer-trimming) for parameters. These genomes will be published to the folder `consensus-seqs`.

### `quast`

Run `QUAST` on the consensus genomes. The results will be published to the folder `QUAST`.

### `callVariants`

Call variant sites (SNPs) use `samtools mpileup` and `bcftools call`. This will also produce files from `bcftools stats`. The VCF files and stats files are published to `sample-variants`.

### `computeStats`

Compute alignment and assembly statistics using `samtools stats` and the python script [`alignment_assembly_stats.py`](../bin/alignment_assembly_stats.py). Coverage plots are created during this process and published to `coverage-plots`.

### `combinedVariants`

Creates a joint VCF using `samtools` and `bcftools`. `combined.vcf` will be published to the output directory. This is skipped by default.

### `mergeAllAssemblies`

Combine all the consensus genomes into one file `combined.fa`, which will be published to the output directory.

### `mergeAssemblyStats`

Parse previous files to create a summary statistics file `combined.stats.tsv`. This will be published to `call_consensus-stats/combined.stats.tsv`.

### `filterAssemblies`

Filter the assemblies based on the criteria given with `--maxNs` and `--minLength`. This outputs `filtered.fa`, published to the output directory, and `filtered.stats.tsv`, published to `call_consensus-stats`.

### `multiqc`

Run MultiQC to collect the reports from `Trim Galore` (`FastQC`), `QUAST`, `samtools` and `bcftools`. The output will be published to the folder `MultiQC`.
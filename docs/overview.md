# Pipeline overview
<!-- MarkdownTOC -->

- [Processes](#processes)
	- [`filterRefReads`](#filterrefreads)
	- [`quantifyERCCs`](#quantifyerccs)
	- [`filterReads`](#filterreads)
	- [`trimReads`](#trimreads)

<!-- /MarkdownTOC -->
## Processes

### `filterRefReads`

Filters host reads by mapping with `minimap2` to the file given with `--ref_host` and keeping unmapped reads. This process can be skipped with the flag `--skip_filter_ref`. The filtered reads will be output into the folder `filtered-reads`.

### `quantifyERCCs`

Quantify ERCC reads by mapping with `minimap2` to the file given with `--ercc_fasta`. The individual output files are saved to the folder `ercc-stats` and the information is summarized in `combined.stats.tsv`.

### `filterReads`

Filter for only SARS-CoV-2 reads. This is accomplished by first mapping with `minimap2` to the file given with `--ref` and only keeping mapped reads. These reads are then run through `kraken2` against the database provided by `--kraken2_db`. Only reads that are assigned uniquely to SARS-CoV-2 are kept.

### `trimReads`
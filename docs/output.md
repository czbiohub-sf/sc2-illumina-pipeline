# Output
<!-- MarkdownTOC -->

- [Summary](#summary)
- [`combined.fa`](#combinedfa)
- [`filtered.fa`](#filteredfa)
- [`MultiQC`](#multiqc)
- [`QUAST`](#quast)
- [`aligned-reads`](#aligned-reads)
- [`call_consensus-stats`](#call_consensus-stats)
	- [`combined.stats.tsv`](#combinedstatstsv)
		- [`sample_name`](#sample_name)
		- [`depth_avg`](#depth_avg)
		- [`mapped_reads`](#mapped_reads)

<!-- /MarkdownTOC -->
## Summary

In the output directory specified by `--outdir`, the following folders and files will be created:

```
outdir
├── combined.fa
├── filtered.fa
├── MultiQC
│   ├── multiqc_data
│   └── multiqc_plots
│       ├── pdf
│       ├── png
│       └── svg
├── QUAST
│   ├── sample1
│   │   ├── aligned_stats
│   │   ├── basic_stats
│   │   ├── contigs_reports
│   │   │   └── minimap_output
│   │   ├── genome_stats
│   │   └── icarus_viewers
│   └── ...
├── aligned-reads
│   ├── sample1.primertrimmed.bam
│   ├── sample1.primertrimmed.bam.bai
│   └── ...
├── call_consensus-stats
│   ├── combined.stats.tsv
│   └── filtered.stats.tsv
├── consensus-seqs
│   ├── sample1.consensus.fa
│   └── ...
├── coverage-plots
│   ├── sample1.depths.png
│   └── ...
├── ercc-stats
│   ├── sample1.ercc_stats
│   └── ...
├── sample-variants
│   ├── sample1.bcftools_stats
│   ├── sample1.vcf.gz
│   ├── sample1.vcf.gz.tbi
│   └── ...
└── trimmed-reads
    ├── sample1_covid_1_val_1.fq.gz
    ├── sample1_covid_2_val_2.fq.gz
    └── ...
```

## `combined.fa`

A multi-FASTA file containing all consensus sequences.

## `filtered.fa`

A multi-FASTA file containing only sequences that pass filtering by `--maxNs` and `--minLength`.

## `MultiQC`

All MultiQC output files.

## `QUAST`

All QUAST output files.

## `aligned-reads`

Primer-trimmed BAM files from alignment of reads to the reference provided by `--ref`. These files are the output of `ivar trim`.

## `call_consensus-stats`

### `combined.stats.tsv`

Summary statistics for each sample run through the pipeline. Each column is described below.

#### `sample_name`

Sample name derived from excluding the read suffixes.

#### `depth_avg`

Average depth across the consensus genome, calculated using `samtools depth`.

#### `mapped_reads`
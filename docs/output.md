# Output
<!-- MarkdownTOC -->

- [Summary](#summary)
- [`combined.fa`](#combinedfa)
- [`filtered.fa`](#filteredfa)
- [`MultiQC/`](#multiqc)
	- [`multiqc_report.html`](#multiqc_reporthtml)
- [`QUAST/`](#quast)
- [`aligned-reads/`](#aligned-reads)
- [`call_consensus-stats/`](#call_consensus-stats)
	- [`combined.stats.tsv`](#combinedstatstsv)
	- [`filtered.stats.tsv`](#filteredstatstsv)
- [`consensus-seqs/`](#consensus-seqs)
- [`coverage-plots/`](#coverage-plots)
- [`ercc-stats/`](#ercc-stats)
- [`filtered-sars2-reads/`](#filtered-sars2-reads)
- [`sample-variants/`](#sample-variants)
- [`trimmed-reads/`](#trimmed-reads)

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
├── filtered-sars2-reads			# optional
│   ├── sample1_covid_1_val_1.fq.gz
│   ├── sample1_covid_2_val_2.fq.gz
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

## `MultiQC/`

All MultiQC output files.

### `multiqc_report.html`

QC report from MultiQC (aggregates stats from quast, samtools, trim-galore, bcftools).

## `QUAST/`

All QUAST output files.

## `aligned-reads/`

Primer-trimmed BAM files from alignment of reads to the reference provided by `--ref`. These files are the output of `ivar trim`.

## `call_consensus-stats/`

### `combined.stats.tsv`

Summary QC statistics for each sample run through the pipeline. Each column is described below. This is the primary QC file.

Column header| Description
-----|-----
`sample_name`| Sample name derived from excluding the read suffixes.
`depth_avg`| Average depth across the consensus genome, calculated using `samtools depth`.
`mapped_reads`| Total number of reads, paired or single, that are mapped.
`total_reads`| Total number of reads in the input file(s) for the sample.
`n_actg`| Total number of nonambiguous bases (ACTG) in the consensus genome.
`n_missing`| Total number of bases called N in the consensus genome. The N character represents any base; these characters are called when the depth at the position is less than `--minDepth`.
`n_gap`| Total number of gaps in the consensus genome.
`n_ambiguous`| Total number of [IUPAC ambiguous (degenerate) bases](https://www.bioinformatics.org/sms/iupac.html). These are called at positions where the most frequent allele does not meet the threshold set by `--ivarFreqThreshold`.
`ref_snps`| Number of SNPs (single nucleotide polymorphisms) in the assembly relative to the reference genome.
`ref_mnps`| Number of MNPs (multiple nucleotide polymorphisms) in the assembly relative to the reference genome.
`ref_indels`| Number of indels with respect to the reference genome.
`depth_q.01`| This is the q-th quantile where `q=.01`. This is the value at which 1% of the depth values are below. If this value is higher, then the depth over the genome is high.
`depth_q.05`| This is the q-th quantile where `q=.05`. This the value at which 5% of the depth values per base are below.
`depth_q.1`| This is the q-th quantile where `q=.1`. This is the value at which 10% of the depth values per base are below.
`depth_q.25`| This is the q-th quantile where `q=.25`. This is the value at which 25% of the depth values per base are below.
`depth_q.5`| This is the q-th quantile where `q=.5`. This is the value at which 50% of the depth values per base are below. This is equivalent to the median depth.
`depth_q.75`| This is the q-th quantile where `q=.75`. This is the value at which 75% of the depth values per base are below.
`depth_frac_above_10x`| This is the proportion of bases with respect to the reference that have at least 10x coverage.
`depth_frac_above_25x`| This is the proportion of bases with respect to the reference that have at least 25x coverage.
`depth_frac_above_50x`| This is the proportion of bases with respect to the reference that have at least 50x coverage.
`depth_frac_above_100x`| This is the proportion of bases with respect to the reference that have at least 100x coverage.
`mapped_paired`| This is the number of reads where the pair maps to the reference.
`paired_inward`| This is the number of reads that are paired and map with an inward orientation (-> <-).
`paired_outward`| This is the number of reads that are paired and map with an outward orientation (<- ->).
`paired_other_orientation`| This is the number of reads that are paired and do not map in either outward or inward orientation (e.g. <- <-).
`ercc_mapped_reads`| This is the number of paired or single reads that map to the file provided by `--ercc_fasta`.
`ercc_mapped_paired`| This is the number of paired reads that map to the file provided by `--ercc_fasta`.

### `filtered.stats.tsv`

This file contains the same columns as in `combined.stats.tsv`, but is filtered only for samples that pass the filters imposed by `--maxNs` and `--minLength`.

## `consensus-seqs/`

All individual sample consensus sequences are in this folder.

## `coverage-plots/`

Line plots of coverage for each sample.

## `ercc-stats/`

The outputs of `samtools stats` on the alignments to the file `--ercc_fasta`.

## `filtered-sars2-reads/`

This folder is only created when the flag `--save_sars2_filtered_reads` is set. It contains untrimmed reads that are host-filtered and map to the reference genome.

## `sample-variants/`

VCFs for each sample with respect to the reference.

## `trimmed-reads/`

Sample reads that have undergone all filtering steps and adapter-trimming.
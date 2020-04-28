#!/usr/bin/env python3

import argparse
import collections
import json
import re
import subprocess
import pysam
from Bio import SeqIO
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("--sample_name")
parser.add_argument("--cleaned_bam")
parser.add_argument("--assembly")
parser.add_argument("--samtools_stats")
parser.add_argument("--ercc_stats")
parser.add_argument("--vcf", help="reference SNPs")
parser.add_argument("--primervcf")
parser.add_argument("--neighborvcf")
parser.add_argument("--neighborfasta")
parser.add_argument("--clades")
parser.add_argument("--out_prefix")
parser.add_argument("--reads", nargs="+")
args = parser.parse_args()

stats = {"sample_name": args.sample_name}

# keep extraction of neighbor name from VCF for compatibility with large pipeline
if args.neighborvcf:
    neighbor_vcf = pysam.VariantFile(args.neighborvcf)
    nearest_neighbor = list(neighbor_vcf.header.contigs)[0]
    stats["nearest_sequence"] = nearest_neighbor

# separate extraction of neighbor name from FASTA if no VCF
elif args.neighborfasta:
    neighbor_fasta = SeqIO.read(args.neighborfasta, 'fasta')
    nearest_neighbor = neighbor_fasta.name
    stats["nearest_sequence"] = nearest_neighbor


if args.cleaned_bam:
    samfile = pysam.AlignmentFile(args.cleaned_bam, "rb")
    ref_len, = samfile.lengths
    depths = [0] * ref_len
    for column in samfile.pileup():
        depths[column.reference_pos] = column.nsegments
    depths = np.array(depths)

    stats["depth_avg"] = depths.mean()
    stats["depth_q.25"] = np.quantile(depths, .25)
    stats["depth_q.5"] = np.quantile(depths, .5)
    stats["depth_q.75"] = np.quantile(depths, .75)
    stats["depth_frac_above_10x"] = (depths >= 10).mean()
    stats["depth_frac_above_25x"] = (depths >= 30).mean()
    stats["depth_frac_above_50x"] = (depths >= 30).mean()
    stats["depth_frac_above_100x"] = (depths >= 100).mean()

    ax = sns.lineplot(np.arange(1, ref_len+1), depths)
    ax.set_title(args.sample_name)
    ax.set(xlabel="position", ylabel="depth")
    plt.yscale("symlog")
    plt.savefig(args.out_prefix + ".depths.png")

seq, = SeqIO.parse(args.assembly, "fasta")
stats["allele_counts"] = dict(collections.Counter(str(seq.seq)))

if args.reads:
    fq_lines = subprocess.run(" ".join(["zcat"] + list(args.reads)) + " | wc -l",
                              shell=True, stdout=subprocess.PIPE).stdout
    stats["total_reads"] = int(int(fq_lines) / 4)

if args.samtools_stats:
    with open(args.samtools_stats) as f:
        sam_stats_re = re.compile(r"SN\s+([^\s].*):\s+(\d+)")
        for line in f:
            matched = sam_stats_re.match(line)
            if matched:
                if matched.group(1) == "reads mapped":
                    stats["mapped_reads"] = int(matched.group(2))
                elif matched.group(1) == "reads mapped and paired":
                    stats["mapped_paired"] = int(matched.group(2))
                elif matched.group(1) == "inward oriented pairs":
                    stats["paired_inward"] = int(matched.group(2)) * 2
                elif matched.group(1) == "outward oriented pairs":
                    stats["paired_outward"] = int(matched.group(2)) * 2
                elif matched.group(1) == "pairs with other orientation":
                    stats["paired_other_orientation"] = int(matched.group(2)) * 2
                # TODO: number of discordant read pairs

if args.ercc_stats:
    with open(args.ercc_stats) as f:
        ercc_stats_re = re.compile(r"SN\s+([^\s].*):\s+(\d+)")
        for line in f:
            matched = ercc_stats_re.match(line)
            if matched:
                if matched.group(1) == "reads mapped":
                    stats["ercc_mapped_reads"] = int(matched.group(2))
                elif matched.group(1) == "reads mapped and paired":
                    stats["ercc_mapped_paired"] = int(matched.group(2))


def countVCF(vcf_file, snpcol, mnpcol, indelcol, statsdict):
    vcf = pysam.VariantFile(vcf_file)
    statsdict[snpcol] = 0
    statsdict[mnpcol] = 0
    statsdict[indelcol] = 0
    for rec in vcf.fetch():
        allele_lens = set([len(a) for a in [rec.ref] + list(rec.alts)])
        if len(allele_lens) > 1:
            statsdict[indelcol] += 1
        else:
            l, = allele_lens
            if l == 1:
                statsdict[snpcol] += 1
            else:
                statsdict[mnpcol] += 1
    return statsdict
if args.vcf:
    stats = {**stats, **countVCF(args.vcf, 'ref_snps', 'ref_mnps', 'ref_indels', stats)}
if args.primervcf:
    stats = {**stats, **countVCF(args.primervcf, 'primer_snps', 'primer_mnps', 'primer_indels', stats)}
if args.neighborvcf:
    stats = {**stats, **countVCF(args.neighborvcf, 'new_snps', 'new_mnps', 'new_indels', stats)}
if args.clades:
    stats["clade"] = []
    with open(args.clades) as f:
        for line in f:
            stats["clade"].append(line.strip())

with open(args.out_prefix + ".stats.json", "w") as f:
    json.dump(stats, f, indent=2)

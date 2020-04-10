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
parser.add_argument("--vcf")
parser.add_argument("--primervcf")
parser.add_argument("--clades")
parser.add_argument("--out_prefix")
parser.add_argument("--reads", nargs="+")
args = parser.parse_args()

stats = {"sample_name": args.sample_name}

samfile = pysam.AlignmentFile(args.cleaned_bam, "rb")
ref_len, = samfile.lengths
depths = [0] * ref_len
for column in samfile.pileup():
    depths[column.reference_pos] = column.nsegments
depths = np.array(depths)
sns.lineplot(np.arange(1, ref_len+1), depths)
plt.yscale("symlog")
plt.savefig(args.out_prefix + ".depths.png")
stats["avg_depth"] = depths.mean()

seq, = SeqIO.parse(args.assembly, "fasta")
stats["allele_counts"] = dict(collections.Counter(str(seq.seq)))

fq_lines = subprocess.run(" ".join(["zcat"] + list(args.reads)) + " | wc -l",
                          shell=True, stdout=subprocess.PIPE).stdout
stats["total_reads"] = int(int(fq_lines) / 4)

with open(args.samtools_stats) as f:
    sam_stats_re = re.compile(r"SN\s+([^\s].*):\s+(\d+)")
    for line in f:
        matched = sam_stats_re.match(line)
        if matched:
            if matched.group(1) == "reads mapped":
                stats["mapped"] = int(matched.group(2))
            elif matched.group(1) == "reads mapped and paired":
                stats["mapped_paired"] = int(matched.group(2))
            elif matched.group(1) == "inward oriented pairs":
                stats["pairs_inward"] = int(matched.group(2)) * 2
            elif matched.group(1) == "outward oriented pairs":
                stats["pairs_outward"] = int(matched.group(2)) * 2
            elif matched.group(1) == "pairs with other orientation":
                stats["pairs_other_orientation"] = int(matched.group(2)) * 2
            # TODO: number of discordant read pairs

vcf = pysam.VariantFile(args.vcf)
stats["snps"] = 0
stats["mnps"] = 0
stats["indels"] = 0
for rec in vcf.fetch():
    allele_lens = set([len(a) for a in [rec.ref] + list(rec.alts)])
    if len(allele_lens) > 1:
        stats["indels"] += 1
    else:
        l, = allele_lens
        if l == 1:
            stats["snps"] += 1
        else:
            stats["mnps"] += 1

primervcf = pysam.VariantFile(args.primervcf)
stats["primer_snps"] = 0
stats["primer_mnps"] = 0
stats["primer_indels"] = 0
for rec in primervcf.fetch():
    allele_lens = set([len(a) for a in [rec.ref] + list(rec.alts)])
    if len(allele_lens) > 1:
        stats["primer_indels"] += 1
    else:
        l, = allele_lens
        if l == 1:
            stats["primer_snps"] += 1
        else:
            stats["primer_mnps"] += 1

stats["clades"] = []
with open(args.clades) as f:
    for line in f:
        stats["clades"].append(line.strip())

with open(args.out_prefix + ".stats.json", "w") as f:
    json.dump(stats, f, indent=2)

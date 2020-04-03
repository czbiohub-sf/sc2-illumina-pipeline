#!/usr/bin/env python3

import argparse
import collections
import json
import re
import pysam
from Bio import SeqIO
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("sample_name")
parser.add_argument("raw_bam")
parser.add_argument("trimmed_filtered_bam")
parser.add_argument("in_assembly")
parser.add_argument("in_sam_stats")
parser.add_argument("out_prefix")
args = parser.parse_args()

raw_samfile = pysam.AlignmentFile(args.raw_bam, "rb")
samfile = pysam.AlignmentFile(args.trimmed_filtered_bam, "rb")

ref_len, = samfile.lengths
depths = [0] * ref_len
for column in samfile.pileup():
    depths[column.reference_pos] = column.nsegments
depths = np.array(depths)

sns.lineplot(np.arange(1, ref_len+1), depths)
plt.yscale("symlog")
plt.savefig(args.out_prefix + ".depths.png")

seq, = SeqIO.parse(args.in_assembly, "fasta")

allele_counts = dict(collections.Counter(str(seq.seq)))

# TODO: number of discordant read pairs
stats = {
    "sample_name": args.sample_name,
    "avg_depth": depths.mean(),
    "allele_counts": allele_counts,
    "total_reads": raw_samfile.mapped + raw_samfile.unmapped,
}

with open(args.in_sam_stats) as f:
    sam_stats_re = re.compile(r"SN\s+([^\s].*):\s+(\d+)")
    for line in f:
        matched = sam_stats_re.match(line)
        if matched:
            if matched.group(1) == "reads mapped":
                stats["mapped"] = int(matched.group(2))
            elif matched.group(1) == "reads mapped and paired":
                stats["mapped_paired"] = int(matched.group(2))

with open(args.out_prefix + ".stats.json", "w") as f:
    json.dump(stats, f, indent=2)

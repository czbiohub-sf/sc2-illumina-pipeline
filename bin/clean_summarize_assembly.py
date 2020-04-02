#!/usr/bin/env python3

import argparse
import collections
import json
import pysam
from Bio import SeqIO
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("sample_name")
parser.add_argument("in_bam")
parser.add_argument("in_assembly")
parser.add_argument("out_prefix")
args = parser.parse_args()

samfile = pysam.AlignmentFile(args.in_bam, "rb")

ref_len, = samfile.lengths
depths = [0] * ref_len
for column in samfile.pileup():
    depths[column.reference_pos] = column.nsegments
depths = np.array(depths)

sns.lineplot(np.arange(1, ref_len+1), depths)
plt.yscale("symlog")
plt.savefig(args.out_prefix + ".depths.png")

seq, = SeqIO.parse(args.in_assembly, "fasta")
seq.id = args.sample_name
seq.name = args.sample_name
seq.description = args.sample_name

SeqIO.write([seq], args.out_prefix + ".fa", "fasta")

allele_counts = dict(collections.Counter(str(seq.seq)))

# TODO: number of discordant read pairs
stats = {
    "sample_name": args.sample_name,
    "avg_depth": depths.mean(),
    "allele_counts": allele_counts
}

with open(args.out_prefix + ".stats.json", "w") as f:
    json.dump(stats, f, indent=2)

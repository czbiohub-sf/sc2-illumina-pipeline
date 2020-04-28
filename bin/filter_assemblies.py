#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--max_n', type=int)
parser.add_argument('--min_len', type=int)
parser.add_argument('--stats')
parser.add_argument('--fasta')
parser.add_argument('--vcf')
parser.add_argument('--out_prefix')
args = parser.parse_args()

stats_df = pd.read_csv(args.stats, sep="\t")

filtered_rows = []
for _, row in stats_df.iterrows():
    if row["n_missing"] <= args.max_n and row["n_actg"] >= args.min_len:
        filtered_rows.append(row)
if filtered_rows:
    filtered_rows = pd.DataFrame(filtered_rows)
else:
    filtered_rows = pd.DataFrame(columns=stats_df.columns)
filtered_rows.to_csv(args.out_prefix + ".stats.tsv", sep="\t", index=False)

samples_to_keep = set(filtered_rows["sample_name"])

def filtered_seqs():
    for seq in SeqIO.parse(args.fasta, "fasta"):
        if seq.id in samples_to_keep:
            yield seq

SeqIO.write(filtered_seqs(), args.out_prefix + ".fa", "fasta")

if args.vcf:
    with open(args.out_prefix + ".vcf", "w") as f:
        subprocess.run(["bcftools", "view",
                        "-s", ",".join(samples_to_keep),
                        "-c", "1",
                        args.vcf],
                       stdout=f)

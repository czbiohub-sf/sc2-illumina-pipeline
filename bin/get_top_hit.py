#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd
import shutil
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument('--sequences', help='FASTA file of reference sequences')
parser.add_argument('--sampleName')
parser.add_argument('--assembly')
parser.add_argument('--default', help='reference sequence if no hits')
args = parser.parse_args()

if len(SeqIO.read(args.assembly, 'fasta')) > 0:
    subprocess.run(" ".join(["blastn", "-db", "gisaid.nt", "-query",
                    args.assembly, "-num_threads", "32", "-out",
                    f"{args.sampleName}.blast.tsv", "-outfmt", "'6",
                    "sacc", "nident", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore'"]), shell=True)

    df = pd.read_csv(f'{args.sampleName}.blast.tsv', sep='\t', names=['sacc', 'nident', 'pident',
                                                          'length', 'mismatch', 'gapopen', 'qstart',
                                                          'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    try:
        top_hit = df['sacc'][0]
        sequences = SeqIO.index(args.sequences, 'fasta')
        top_hit_seq = sequences[top_hit]
        SeqIO.write(top_hit_seq, 'nearest_gisaid.fasta', 'fasta')
    except IndexError:
        shutil.copyfile(args.default, 'nearest_gisaid.fasta')
else:
    shutil.copyfile(args.default, 'nearest_gisaid.fasta')
    with open(f'{args.sampleName}.blast.tsv', 'w+') as f:
        f.write('Empty query sequence.')

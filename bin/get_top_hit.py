#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd
import shutil


parser = argparse.ArgumentParser()
parser.add_argument('--tsv', help='blast output TSV')
parser.add_argument('--sequences', help='FASTA file of reference sequences')
parser.add_argument('--default', help='reference sequence if no hits')
args = parser.parse_args()

df = pd.read_csv(args.tsv, sep='\t', names=['sacc', 'nident', 'pident',
                                                      'length', 'mismatch', 'gapopen', 'qstart',
                                                      'qend', 'sstart', 'send', 'evalue', 'bitscore'])
try:
    top_hit = df['sacc'][0]
    sequences = SeqIO.index(args.sequences, 'fasta')
    top_hit_seq = sequences[top_hit]
    SeqIO.write(top_hit_seq, 'nearest_gisaid.fasta', 'fasta')
except IndexError:
    shutil.copyfile(args.default, 'nearest_gisaid.fasta')


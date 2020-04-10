#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--csv', help='sourmash search output CSV')
parser.add_argument('--sequences', help='FASTA file of reference sequences')
args = parser.parse_args()

sourmash_out = pd.read_csv(args.csv)
top_hit = sourmash_out['name'][0]

sequences = SeqIO.index(args.sequences, 'fasta')

top_hit_seq = sequences[top_hit]

SeqIO.write(top_hit_seq, 'nearest_gisaid.fasta', 'fasta')
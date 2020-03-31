#!/usr/bin/env python3


import pandas as pd
import argparse
from Bio import SeqIO
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--file', help='quast unaligned_report.tsv file')
parser.add_argument('-n', help='max number of Ns to allow', type=int, default=100)
parser.add_argument('--assembly', help='assembly FASTA')
parser.add_argument('-l', help='min length for passing QC', type=int, default=29000)
args = parser.parse_args()

df = pd.read_csv(args.file, sep='\t')
num_Ns = df.iloc[4, 1]

seqrecords = SeqIO.parse(args.assembly, 'fasta')
seq_length = sum(len(s) for s in seqrecords)

if num_Ns <= args.n and seq_length >= args.l:
	shutil.copy(args.assembly, 'passed_QC/')
else:
	shutil.copy(args.assembly, 'failed_QC/')

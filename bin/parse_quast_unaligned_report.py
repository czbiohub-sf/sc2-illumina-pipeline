#!/usr/bin/env python3


import pandas as pd
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--unaligned_report', help='quast unaligned_report.tsv')
parser.add_argument('--report', help='quast report.tsv')
parser.add_argument('-n', help='max number of Ns to allow', type=int, default=100)
parser.add_argument('--assembly', help='assembly FASTA')
parser.add_argument('-l', help='min length for passing QC', type=int, default=29000)
args = parser.parse_args()

df = pd.read_csv(args.unaligned_report, sep='\t')
num_Ns = df.iloc[4, 1]

report_df = pd.read_csv(args.report, sep='\t')
seq_length = int(report_df.iloc[6,1])

if num_Ns <= args.n and seq_length >= args.l:
	shutil.copy(args.assembly, 'passed_QC/')
else:
	shutil.copy(args.assembly, 'failed_QC/')

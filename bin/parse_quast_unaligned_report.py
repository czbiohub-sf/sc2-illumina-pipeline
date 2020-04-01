#!/usr/bin/env python3


import pandas as pd
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--report', help='quast report.tsv')
parser.add_argument('-n', help='max number of Ns to allow', type=int, default=100)
parser.add_argument('--assembly', help='assembly FASTA')
parser.add_argument('-l', help='min length for passing QC', type=int, default=29000)
args = parser.parse_args()

report_df = pd.read_csv(args.report, sep='\t')
total_length = int(report_df.loc[report_df['Assembly']=="Total length (>= 0 bp)"].iloc[0,1])
num_Ns = round(float(report_df.loc[report_df['Assembly']=="# N's per 100 kbp"].values[0][1])*(total_length/100000))
seq_length = total_length - num_Ns

if num_Ns <= args.n and seq_length >= args.l:
	shutil.copy(args.assembly, 'passed_QC/')
else:
    shutil.copy(args.assembly, 'failed_QC/')


#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import re
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--public_identifiers')
parser.add_argument('--sampleName')
parser.add_argument('--in_fa')
args = parser.parse_args()

df = pd.read_csv(args.public_identifiers, sep='\t')
try:
  submission_id = df[df['sample_name']==args.sampleName]['gisaid_name'].values[0]
  submission_id = re.search(pattern='hCoV-19/(.*)$', string=submission_id).group(1)
except IndexError:
  submission_id = args.sampleName
seq = SeqIO.read(args.in_fa, 'fasta')
seq.id = submission_id
seq.description = submission_id
seq.name = submission_id
try:
  filename = re.search(pattern='USA/(.*)/2020', string=submission_id).group(1)
except AttributeError:
  filename = submission_id
SeqIO.write(seq, f'{filename}.fasta', 'fasta')
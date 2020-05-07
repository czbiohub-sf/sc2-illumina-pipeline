#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import re
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--public_identifiers')
parser.add_argument('--in_fa')
args = parser.parse_args()

df = pd.read_csv(args.public_identifiers, sep='\t')

sequences = list(SeqIO.parse(args.in_fa, 'fasta'))
new_sequences = []
for seq in sequences:
    try:
        publicID = df[df['sample_name']==seq.id]['gisaid_name'].values[0]
        publicID = re.search(pattern='hCoV-19/(.*)$', string=publicID).group(1)
    except IndexError:
        publicID = seq.id
    seq.id = publicID
    seq.name = publicID
    seq.description = publicID
    new_sequences.append(seq)
SeqIO.write(new_sequences, 'renamed_sequences.fasta', 'fasta')
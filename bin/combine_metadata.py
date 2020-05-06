#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--sample_metadata')
parser.add_argument('--global_metadata')
args = parser.parse_args()

sample_meta = pd.read_csv(args.sample_metadata, sep='\t')
global_meta = pd.read_csv(args.global_metadata, sep='\t')

df = pd.concat([global_meta, sample_meta], join='inner', ignore_index=True)
df.to_csv('metadata.tsv', sep='\t', index=False)
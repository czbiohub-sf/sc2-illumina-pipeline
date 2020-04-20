#!/usr/bin/env python3


import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--prev_metadata', '-pm', help='Previous metadata to include')
parser.add_argument('--prev_sequences', '-ps', help='Previous sequences to include')
parser.add_argument('--new_sequences', '-ns', help='New sample sequences to include')
parser.add_argument('--date', '-d', help='Date for the samples (e.g. 2020-03-03)')
parser.add_argument('--date_tsv', help='TSV with per sample dates')
parser.add_argument('--region', '-r', help='Region for the samples (e.g. North America)')
parser.add_argument('--country', '-c', help='Country for the samples (e.g. USA)')
parser.add_argument('--division', '-div', help='Division for the samples (e.g. California)')
parser.add_argument('--location', '-loc', help='Location for the samples (e.g. San Francisco County)')
parser.add_argument('--country_exposure', '-cexp', help='Country exposure (default: ?)', default='?')
parser.add_argument('--div_exposure', '-divexp', help='Division exposure (default: ?)', default='?')
parser.add_argument('--originating_lab', '-origlab', help='Originating lab (e.g. Biohub)')
parser.add_argument('--submitting_lab', '-sublab', help='Submitting lab (e.g. Biohub)')
parser.add_argument('--date_submitted', '-subdate', help='Date of submission')
args = parser.parse_args()

# Build the new metadata.tsv
old_metadata = pd.read_csv(args.prev_metadata, sep='\t')
new_metadata = pd.DataFrame()

sequences = list(SeqIO.parse(args.new_sequences, 'fasta'))

strains = [s.id for s in sequences]
num_strains = len(strains)
unknown_col = ['?']*num_strains
new_metadata['strain'] = strains
new_metadata['virus'] = ['ncov']*num_strains
new_metadata['gisaid_epi_isl'] = unknown_col
new_metadata['genbank_accession'] = unknown_col
new_metadata['region'] = [args.region]*num_strains
new_metadata['country'] = [args.country]*num_strains
new_metadata['division'] = [args.division]*num_strains
new_metadata['location'] = [args.location]*num_strains
new_metadata['country_exposure'] = [args.country_exposure]*num_strains
new_metadata['division_exposure'] = [args.div_exposure]*num_strains
new_metadata['segment'] = ['genome']*num_strains
new_metadata['length'] = [len(s) for s in sequences]
new_metadata['host'] =  ['Human']*num_strains
new_metadata['age'] = unknown_col
new_metadata['sex'] = unknown_col
new_metadata['originating_lab'] = [args.originating_lab]*num_strains
new_metadata['submitting_lab'] = [args.submitting_lab]*num_strains
new_metadata['authors'] = unknown_col
new_metadata['url'] = unknown_col
new_metadata['title'] = unknown_col
new_metadata['date_submitted'] = [args.date_submitted]*num_strains
if args.date_tsv:
    date_df = pd.read_csv(args.date_tsv, sep='\t')
    new_metadata = new_metadata.merge(date_df, how='left', on='strain')
    new_metadata = new_metadata.fillna({'date': args.date})
else:
    new_metadata['date'] = [args.date]*num_strains
df = pd.concat([old_metadata, new_metadata], sort=False)
df.to_csv('metadata.tsv', sep='\t', index=False)

# Build the new sequences.fasta
old_sequences = list(SeqIO.parse(args.prev_sequences, 'fasta'))
SeqIO.write(old_sequences + sequences, 'all_sequences.fasta', 'fasta')

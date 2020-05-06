#!/usr/bin/env python3


import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sample_sequences', '-s', help='Sample sequences to include')
parser.add_argument('--output', '-o', help='output file')
parser.add_argument('--date', '-d', help='Date for the samples (e.g. 2020-03-03)')
parser.add_argument('--region', '-r', help='Region for the samples (e.g. North America)')
parser.add_argument('--country', '-c', help='Country for the samples (e.g. USA)')
parser.add_argument('--division', '-div', help='Division for the samples (e.g. California)')
parser.add_argument('--location', '-loc', help='Location for the samples (e.g. San Francisco)')
parser.add_argument('--submitting_lab', '-sublab', help='Submitting lab (e.g. Biohub)')
parser.add_argument('--date_submitted', '-subdate', help='Date of submission')
args = parser.parse_args()

sequences = list(SeqIO.parse(args.sample_sequences, 'fasta'))

# Build the new metadata.tsv
new_metadata = pd.DataFrame()

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
new_metadata['region_exposure'] = [args.region]*num_strains
new_metadata['country_exposure'] = [args.country]*num_strains
new_metadata['division_exposure'] = [args.division]*num_strains
new_metadata['segment'] = ['genome']*num_strains
new_metadata['length'] = [len(s) for s in sequences]
new_metadata['host'] =  ['Human']*num_strains
new_metadata['age'] = unknown_col
new_metadata['sex'] = unknown_col
new_metadata['originating_lab'] = unknown_col
new_metadata['submitting_lab'] = [args.submitting_lab]*num_strains
new_metadata['authors'] = unknown_col
new_metadata['url'] = unknown_col
new_metadata['title'] = unknown_col
new_metadata['date_submitted'] = [args.date_submitted]*num_strains
new_metadata['date'] = [args.date]*num_strains

new_metadata.to_csv(args.output, sep='\t', index=False)

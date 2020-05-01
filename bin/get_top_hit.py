#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd
import shutil
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument('--sequences', help='FASTA file of reference sequences')
parser.add_argument('--sampleName')
parser.add_argument('--assembly')
parser.add_argument('--minLength', type=int)
parser.add_argument('--default', help='reference sequence if no hits')
parser.add_argument('--metadata', help='metadata TSV with strain names and divisions')
args = parser.parse_args()

if args.metadata:
    division_metadata = pd.read_csv(args.metadata, sep='\t', usecols=['strain', 'country', 'division', 'date'])
    # only keep strains that have a full date
    division_metadata = division_metadata[(division_metadata['date'].str.len() > 7)]
    division_metadata = division_metadata[~division_metadata['date'].str.contains('X')]
    division_metadata['date'] = division_metadata['date'].astype('datetime64[ns]')
    def check_level(row):
        if row['division'] == 'California':
            return 'CA'
        elif row['country'] == 'USA':
            return row['division']
        else:
            return row['country']
    division_metadata['level'] = division_metadata.apply(check_level, axis=1)

# check that the query is not empty
if len(SeqIO.read(args.assembly, 'fasta')) > 0:
    subprocess.run(" ".join(["blastn", "-db", "blast_seqs.nt", "-query",
                    args.assembly, "-num_threads", "32", "-out",
                    f"{args.sampleName}.blast.tsv", "-outfmt", "'6",
                    "sacc", "nident", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore'"]), shell=True)

    df = pd.read_csv(f'{args.sampleName}.blast.tsv', sep='\t', names=['sacc', 'nident', 'pident',
                                                          'length', 'mismatch', 'gapopen', 'qstart',
                                                          'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    if df['length'].max() < args.minLength/2:
        # One HSP should extend to at least half the minimum length since we only expect a few SNPS
        # Otherwise the alignment is probably bad
        shutil.copyfile(args.default, f'{args.sampleName}_nearest_blast.fasta')
    else:
        try:
            # In case the TSV is empty
            top_hit = df[df['bitscore']==df['bitscore'].max()]['sacc']
            sequences = SeqIO.index(args.sequences, 'fasta')
            if args.metadata:
                top_hit_metadata = division_metadata[division_metadata['strain'].isin(top_hit)]
                # sort by date
                top_hit_metadata = top_hit_metadata.sort_values(by='date', ascending=False)
                levels = top_hit_metadata['level'].unique()
                to_concat = []
                for l in levels:
                    if l == 'CA':
                        CA_hits = top_hit_metadata[top_hit_metadata['level']==l]
                        to_concat.append(CA_hits)
                    else:
                        # take the most recent sample
                        sampled_hit = top_hit_metadata[top_hit_metadata['level']==l][:1]
                        to_concat.append(sampled_hit)
                top_hit_metadata = pd.concat(to_concat)
                top_hit_seqs = [sequences[h] for h in top_hit_metadata['strain']]
            else:
                top_hit_seqs = [sequences[h] for h in top_hit]
            SeqIO.write(top_hit_seqs, f'{args.sampleName}_nearest_blast.fasta', 'fasta')
        except IndexError:
            shutil.copyfile(args.default, f'{args.sampleName}_nearest_blast.fasta')
else:
    shutil.copyfile(args.default, f'{args.sampleName}_nearest_blast.fasta')
    with open(f'{args.sampleName}.blast.tsv', 'w+') as f:
        f.write('Empty query sequence.')

#!/usr/bin/env python3

import argparse
from Bio import SeqIO, SeqFeature
import pandas as pd
import numpy as np
import pysam


parser = argparse.ArgumentParser()
parser.add_argument('--reference', help='Genbank reference')
parser.add_argument('--vcf', help='Nucleotide VCF for sample')
parser.add_argument('--sample', help='Sample name')
parser.add_argument('--clades', help='TSV with clade information from Nextstrain')
args = parser.parse_args()

reference = SeqIO.read(args.reference, 'genbank')
CDS_features = filter(lambda x: x.type=='CDS', reference.features)
vcf = pysam.VariantFile(args.vcf)
# copy the reference sequence to a MutableSeq object and change to alternate alleles
alternate = reference.seq.tomutable()

for rec in vcf.fetch():
    pos = rec.pos - 1
    ref_allele = rec.ref
    # TODO: do we need to handle case with more than 1 alt?
    alt_allele = rec.alts[0]
    if len(ref_allele) != 1 or len(alt_allele) != 1:
        # skip non-SNPs
        continue
    else:
        alternate[pos] = alt_allele

sample_aa = {}
for cds in CDS_features:
    if cds.qualifiers['gene'][0] == 'orf1ab':
        orf1a = SeqFeature.SeqFeature(
            location=cds.location.parts[0],
            strand=cds.strand,
            qualifiers={'gene': 'ORF1a'},
            type=cds.type).extract(alternate).translate(gap='-')
        orf1b = SeqFeature.SeqFeature(
            location=cds.location.parts[0],
            strand=cds.strand,
            qualifiers={'gene': 'ORF1b'},
            type=cds.type).extract(alternate).translate(gap='-')
        sample_aa['ORF1a'] = orf1a
        sample_aa['ORF1b'] = orf1b
    else:
        sample_aa[cds.qualifiers['gene'][0]] = cds.extract(alternate).translate(gap='-')
# Create dictionary of clades defined by alleles
clades_df = pd.read_csv(args.clades, sep='\t')
clades = {}

for index, row in clades_df.iterrows():
    allele = (row['gene'], row['site']-1, row['alt'])
    if row['clade'] in clades:
        clades[row['clade']].append(allele)
    else:
        clades[row['clade']] = [allele]

# Store clades in a set in case more than one clade has all conditions fulfilled
sample_clade = set()

for clade in clades:
    for allele in clades[clade]:
        gene, pos, alt = allele
        if gene == 'nuc':
            if (alternate[pos] == alt) and (allele == clades[clade][-1]):
                sample_clade.add(clade)
            elif alternate[pos] == alt:
                continue
            else:
                break
        else:
            if (sample_aa[gene][pos] == alt) and (allele == clades[clade][-1]):
                sample_clade.add(clade)
            elif sample_aa[gene][pos] == alt:
                continue
            else:
                break

if not sample_clade:
    sample_clade = ['unknown']

with open(f'{args.sample}.clades', 'w+') as f:
    for clade in sorted(sample_clade):
        print(clade, file=f)

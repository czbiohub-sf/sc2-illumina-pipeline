#!/usr/bin/env python3

from Bio import SeqIO
import subprocess
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--assembly')
parser.add_argument('--noreads', action='store_true')
parser.add_argument('--sample')
parser.add_argument('--ref')
parser.add_argument('--threads')
parser.add_argument('--bam')
parser.add_argument('--R1')
parser.add_argument('--R2')
args = parser.parse_args()

sequence = SeqIO.read(args.assembly, 'fasta')

if len(sequence) > 0:
	if args.noreads:
		subprocess.run(f"quast --min-contig 0 -o {args.sample} -r {args.ref} -t {args.threads} --ref-bam {args.bam} {args.assembly}", shell=True)
	else:
		subprocess.run(f"quast --min-contig 0 -o {args.sample} -r {args.ref} -t {args.threads} --ref-bam {args.bam} {args.assembly} -1 {args.R1} -2 {args.R2}", shell=True)
else:
	os.makedirs(args.sample, exist_ok=True)
	with open(f"{args.sample}/{args.sample}.txt", 'w+') as f:
		f.write(f"{args.sample} is empty.")

#!/usr/bin/env python3

import argparse
import pandas as pd
import pysam

parser = argparse.ArgumentParser()
parser.add_argument("--in_vcf")
parser.add_argument("--out_csv")
args = parser.parse_args()

snps_df = []
for i, rec in enumerate(pysam.VariantFile(args.in_vcf)):
    if len(rec.alleles) == 2 and all(a in list("ACGT") for a in rec.alleles):
        for s in rec.samples.values():
            a, = s["GT"]
            if a is None:
                ad0, ad1 = 0, 0
            else:
                ad0, ad1 = s["AD"]
                if ad1 is None:
                    ad1 = 0
            snps_df.append({
                "sample": s.name,
                "index": i,
                "pos": rec.pos,
                "gt": a,
                "ad0": ad0,
                "ad1": ad1
            })

snps_df = pd.DataFrame(
    snps_df, columns=["sample", "index", "pos", "gt", "ad0", "ad1"])
snps_df.to_csv(args.out_csv, index=False, na_rep="NA")

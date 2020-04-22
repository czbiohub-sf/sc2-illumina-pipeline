#!/usr/bin/env python3

import sys
import json
import collections
import pandas as pd

rows = []

if sys.argv[1] == 'core':
    col_keys = ["sample_name", "depth_avg", "mapped_reads",
                "total_reads", "n_actg", "n_missing", "n_gap",
                "n_ambiguous"]
elif sys.argv[1] == 'analysis':
    col_keys = ["sample_name", "clade", "n_actg", "n_missing", "n_gap",
                "n_ambiguous", "nearest_sequence",
                "ref_snps", "ref_mnps", "ref_indels"]
elif sys.argv[1] == 'all':
    col_keys = ["sample_name", "clade", "depth_avg", "mapped_reads",
                "total_reads", "n_actg", "n_missing", "n_gap",
                "n_ambiguous", "nearest_sequence",
                "new_snps", "new_mnps", "new_indels"]

for fname in sys.argv[2:]:
    with open(fname) as f:
        row = json.load(f)
    allele_counts = collections.Counter()
    for k, v in row.pop("allele_counts").items():
        allele_counts[k.upper()] += v
    row["n_actg"] = sum(v for k, v in allele_counts.items() if k in "ACTGU")
    row["n_missing"] = allele_counts["N"]
    row["n_gap"] = allele_counts["-"]
    row["n_ambiguous"] = sum(v for k, v in allele_counts.items() if k not in "ACTGUN-")
    if (sys.argv[1] == 'analysis') or (sys.argv[1] == 'all'):
        row["clade"] = ";".join(row["clade"])

    reordered_row = {}
    for key in col_keys:
        reordered_row[key] = row.pop(key)
    reordered_row.update(row)
    rows.append(reordered_row)

pd.DataFrame(rows).to_csv(sys.stdout, index=False, sep="\t", float_format="%.3f")

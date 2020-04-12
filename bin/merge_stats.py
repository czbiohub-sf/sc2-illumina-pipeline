#!/usr/bin/env python3

import sys
import json
import collections
import pandas as pd

rows = []
for fname in sys.argv[1:]:
    with open(fname) as f:
        row = json.load(f)
    allele_counts = collections.Counter()
    for k, v in row.pop("allele_counts").items():
        allele_counts[k.upper()] += v
    row["n_actg"] = sum(v for k, v in allele_counts.items() if k in "ACTGU")
    row["n_missing"] = allele_counts["N"]
    row["n_gap"] = allele_counts["-"]
    row["n_ambiguous"] = sum(v for k, v in allele_counts.items() if k not in "ACTGUN-")
    row["clade"] = ";".join(row["clade"])

    reordered_row = {}
    for key in [
            "sample_name", "clade", "depth_avg", "mapped_reads",
            "total_reads", "n_actg", "n_missing", "n_gap",
            "n_ambiguous", "nearest_sequence",
            "new_snps", "new_mnps", "new_indels"]:
        reordered_row[key] = row.pop(key)
    reordered_row.update(row)
    rows.append(reordered_row)

pd.DataFrame(rows).to_csv(sys.stdout, index=False, sep="\t", float_format="%.3f")

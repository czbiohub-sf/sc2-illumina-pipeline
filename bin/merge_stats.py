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
    rows.append(row)

pd.DataFrame(rows).to_csv(sys.stdout, index=False, sep="\t", float_format="%.1f")

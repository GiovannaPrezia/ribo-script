#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) != 3:
    sys.exit(
        "Usage: python 02_filter_smorfs.py <ribotricer_tsv> <output_tsv>"
    )

input_file = sys.argv[1]
output_file = sys.argv[2]

print(f"Loading: {input_file}")

df = pd.read_csv(
    input_file,
    sep="\t"
)

# ==========================================================
# FIND ORF LENGTH COLUMN
# ==========================================================

possible_cols = [
    "length",
    "orf_length",
    "ORF_length",
    "length_aa"
]

length_col = None

for col in possible_cols:
    if col in df.columns:
        length_col = col
        break

if length_col is None:
    raise ValueError(
        f"Could not identify ORF length column.\n"
        f"Columns found:\n{list(df.columns)}"
    )

print(f"Using length column: {length_col}")

# ==========================================================
# FILTER smORFs
# ==========================================================

smorfs = df[
    (df[length_col] >= 20) &
    (df[length_col] <= 150)
].copy()

print(f"Total ORFs      : {len(df):,}")
print(f"smORFs (20-150) : {len(smorfs):,}")

smorfs.to_csv(
    output_file,
    sep="\t",
    index=False
)

print(f"Saved: {output_file}")

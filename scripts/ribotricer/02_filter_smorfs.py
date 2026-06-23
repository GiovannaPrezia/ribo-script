#!/usr/bin/env python3

import sys
import pandas as pd


def find_column(columns, candidates):
    for col in candidates:
        if col in columns:
            return col
    return None


def main():
    if len(sys.argv) != 3:
        sys.exit(
            "Usage: python 02_filter_smorfs.py <ribotricer_orfs.tsv> <output_smorfs.tsv>"
        )

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    print(f"Loading Ribotricer ORFs: {input_file}")

    df = pd.read_csv(input_file, sep="\t")

    length_col = find_column(
        df.columns,
        [
            "ORF_length",
            "orf_length",
            "length",
            "length_nt",
            "orf_len",
        ],
    )

    if length_col is None:
        raise ValueError(
            "Could not identify ORF length column. "
            f"Columns found: {list(df.columns)}"
        )

    df[length_col] = pd.to_numeric(df[length_col], errors="coerce")

    # Ribotricer usually reports ORF length in nucleotides.
    # If values look nucleotide-scale, convert to amino acids.
    median_length = df[length_col].dropna().median()

    if median_length > 300:
        df["ORF_length_nt"] = df[length_col]
        df["ORF_length_aa"] = (df[length_col] / 3).round(0).astype("Int64")
    else:
        df["ORF_length_aa"] = df[length_col].round(0).astype("Int64")
        df["ORF_length_nt"] = df["ORF_length_aa"] * 3

    smorfs = df[
        (df["ORF_length_aa"] >= 20) &
        (df["ORF_length_aa"] <= 150)
    ].copy()

    print(f"Total ORFs       : {len(df):,}")
    print(f"smORFs 20-150 aa : {len(smorfs):,}")

    smorfs.to_csv(output_file, sep="\t", index=False)

    print(f"Saved: {output_file}")


if __name__ == "__main__":
    main()

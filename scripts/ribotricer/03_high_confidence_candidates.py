#!/usr/bin/env python3

import sys
import os
import glob
import numpy as np
import pandas as pd


MIN_REPLICATES = 2
MIN_PHASE_SCORE = 0.50
MIN_READ_COUNT = 10


def find_column(columns, candidates):
    for col in candidates:
        if col in columns:
            return col
    return None


def main():
    if len(sys.argv) != 3:
        sys.exit(
            "Usage: python 03_high_confidence_candidates.py "
            "<ribotricer_mode_dir> <output_tsv>"
        )

    input_dir = sys.argv[1]
    output_file = sys.argv[2]

    files = glob.glob(
        os.path.join(input_dir, "*_ranked_lncrna_smorfs.tsv")
    )

    if len(files) == 0:
        sys.exit(f"No ranked lncRNA-smORF files found in: {input_dir}")

    dfs = []

    for file in files:
        df = pd.read_csv(file, sep="\t")

        if "sample" not in df.columns:
            sample = os.path.basename(file).replace("_ranked_lncrna_smorfs.tsv", "")
            df["sample"] = sample

        dfs.append(df)

    all_df = pd.concat(dfs, ignore_index=True)

    id_col = find_column(
        all_df.columns,
        ["ORF_ID", "orf_id", "orf_name", "id"]
    )

    if id_col is None:
        id_col = "transcript_id"

    phase_col = find_column(
        all_df.columns,
        [
            "phase_score",
            "phaseScore",
            "periodicity_score",
            "phase_score_valid_codons",
            "phase_component"
        ]
    )

    read_col = find_column(
        all_df.columns,
        [
            "read_count",
            "reads",
            "count",
            "RPF_count",
            "valid_codons"
        ]
    )

    if phase_col is None:
        all_df["phase_for_summary"] = 0
    else:
        all_df["phase_for_summary"] = pd.to_numeric(
            all_df[phase_col],
            errors="coerce"
        ).fillna(0)

    if read_col is None:
        all_df["reads_for_summary"] = 0
    else:
        all_df["reads_for_summary"] = pd.to_numeric(
            all_df[read_col],
            errors="coerce"
        ).fillna(0)

    keep_cols = [
        "transcript_id",
        "gene_id_y",
        "gene_name_y",
        "gene_type_y",
        "transcript_type_y",
        "ORF_length_aa",
        "ORF_length_nt",
        "chrom",
        "strand",
        "start_codon"
    ]

    available_keep_cols = [
        col for col in keep_cols
        if col in all_df.columns and col != id_col
    ]

    summary = (
        all_df
        .groupby(id_col, dropna=False)
        .agg(
            n_replicates=("sample", "nunique"),
            samples=("sample", lambda x: ",".join(sorted(set(map(str, x))))),
            mean_phase_score=("phase_for_summary", "mean"),
            max_phase_score=("phase_for_summary", "max"),
            mean_read_count=("reads_for_summary", "mean"),
            max_read_count=("reads_for_summary", "max"),
            n_calls=("sample", "count")
        )
        .reset_index()
    )

    annotation = (
        all_df[[id_col] + available_keep_cols]
        .drop_duplicates(subset=[id_col])
    )

    summary = summary.merge(
        annotation,
        on=id_col,
        how="left"
    )

    summary["ranking_score"] = (
        summary["max_phase_score"] +
        np.log10(summary["max_read_count"] + 1) +
        summary["n_replicates"]
    )

    summary["high_confidence"] = (
        (summary["n_replicates"] >= MIN_REPLICATES) &
        (summary["max_phase_score"] >= MIN_PHASE_SCORE) &
        (summary["max_read_count"] >= MIN_READ_COUNT)
    )

    high_confidence = summary[
        summary["high_confidence"]
    ].copy()

    high_confidence = high_confidence.sort_values(
        "ranking_score",
        ascending=False
    )

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    summary_file = output_file.replace(
        "_high_confidence_candidates.tsv",
        "_lncrna_smorfs_reproducibility_summary.tsv"
    )

    summary.to_csv(
        summary_file,
        sep="\t",
        index=False
    )

    high_confidence.to_csv(
        output_file,
        sep="\t",
        index=False
    )

    print(f"Total unique lncRNA-smORFs : {len(summary):,}")
    print(f"High-confidence candidates : {len(high_confidence):,}")
    print(f"Criteria: n_replicates>={MIN_REPLICATES}, max_phase_score>={MIN_PHASE_SCORE}, max_read_count>={MIN_READ_COUNT}")
    print(f"Saved summary: {summary_file}")
    print(f"Saved high-confidence: {output_file}")


if __name__ == "__main__":
    main()

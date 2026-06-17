#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def load_length_tables(input_dir: Path, suffix: str, source: str) -> pd.DataFrame:
    rows = []

    for file in sorted(input_dir.glob(f"*{suffix}")):
        name = file.name.replace(suffix, "")

        # Ex: iPSC_DAY21_Ger_rep1.all_lengths
        parts = name.split(".")
        sample = parts[0]
        mode = parts[1] if len(parts) > 1 else "unknown"

        df = pd.read_csv(file, sep="\t", names=["length", "count"])
        df["sample"] = sample
        df["mode"] = mode
        df["source"] = source
        rows.append(df)

    if not rows:
        return pd.DataFrame(columns=["length", "count", "sample", "mode", "source"])

    return pd.concat(rows, ignore_index=True)


def plot_read_lengths(df: pd.DataFrame, out_png: Path, title: str):
    if df.empty:
        print(f"[WARN] No data for {title}")
        return

    samples = sorted(df["sample"].unique())
    modes = sorted(df["mode"].unique())

    n_panels = len(samples) * len(modes)
    fig, axes = plt.subplots(
        n_panels,
        1,
        figsize=(10, max(3, 2.8 * n_panels)),
        sharex=True
    )

    if n_panels == 1:
        axes = [axes]

    i = 0
    for sample in samples:
        for mode in modes:
            ax = axes[i]
            sub = df[(df["sample"] == sample) & (df["mode"] == mode)]

            ax.bar(sub["length"], sub["count"])
            ax.set_title(f"{sample} — {mode}")
            ax.set_ylabel("Reads")

            if not sub.empty:
                peak = sub.loc[sub["count"].idxmax()]
                ax.axvline(peak["length"], linestyle="--", linewidth=1)
                ax.text(
                    peak["length"],
                    peak["count"],
                    f"peak: {int(peak['length'])} nt",
                    ha="left",
                    va="bottom",
                    fontsize=9
                )

            i += 1

    axes[-1].set_xlabel("Read length (nt)")
    fig.suptitle(title, y=1.002, fontsize=14)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_summary(df: pd.DataFrame, out_tsv: Path):
    if df.empty:
        return

    summary = (
        df.sort_values("count", ascending=False)
        .groupby(["sample", "mode", "source"])
        .first()
        .reset_index()
        .rename(columns={"length": "dominant_length", "count": "dominant_length_count"})
    )

    total = (
        df.groupby(["sample", "mode", "source"])["count"]
        .sum()
        .reset_index()
        .rename(columns={"count": "total_reads"})
    )

    summary = summary.merge(total, on=["sample", "mode", "source"])
    summary["dominant_length_fraction"] = (
        summary["dominant_length_count"] / summary["total_reads"]
    )

    summary.to_csv(out_tsv, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc_dir", required=True)
    parser.add_argument("--project", required=True)
    args = parser.parse_args()

    qc_dir = Path(args.qc_dir)
    fig_dir = qc_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    fastq_df = load_length_tables(
        qc_dir / "read_lengths_fastq",
        ".fastq_read_lengths.tsv",
        "FASTQ_clean"
    )

    bam_df = load_length_tables(
        qc_dir / "read_lengths_bam",
        ".bam_read_lengths.tsv",
        "BAM_aligned"
    )

    all_df = pd.concat([fastq_df, bam_df], ignore_index=True)

    plot_read_lengths(
        fastq_df,
        fig_dir / f"{args.project}_RPF_read_lengths_FASTQ_clean.png",
        "Clean FASTQ read length distribution"
    )

    plot_read_lengths(
        bam_df,
        fig_dir / f"{args.project}_RPF_read_lengths_BAM_aligned.png",
        "Aligned BAM read length distribution"
    )

    make_summary(
        all_df,
        qc_dir / f"{args.project}_RPF_read_length_summary.tsv"
    )

    print("RPF QC plots generated:")
    print(fig_dir)


if __name__ == "__main__":
    main()

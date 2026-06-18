#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("--sample", required=True)
parser.add_argument("--qc_dir", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

sample = args.sample
qc_dir = Path(args.qc_dir)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

files = {
    "trimmed": qc_dir / f"{sample}.length_distribution_trimmed.txt",
    "noN_noPolyG_all_lengths": qc_dir / f"{sample}.length_distribution_noN_noPolyG_all_lengths.txt",
    "noN_noPolyG_28_36": qc_dir / f"{sample}.length_distribution_noN_noPolyG_28_36.txt",
}

for label, file in files.items():
    if not file.exists():
        print(f"[WARN] Missing file: {file}")
        continue

    df = pd.read_csv(file, sep=r"\s+", names=["count", "length"])

    plt.figure(figsize=(7, 4))
    plt.bar(df["length"], df["count"])

    peak = df.loc[df["count"].idxmax()]
    plt.axvline(peak["length"], linestyle="--", linewidth=1)
    plt.text(
        peak["length"],
        peak["count"],
        f"peak: {int(peak['length'])} nt",
        fontsize=9,
        va="bottom"
    )

    plt.title(f"{sample} — {label}")
    plt.xlabel("Read length (nt)")
    plt.ylabel("Reads")
    plt.tight_layout()

    output = outdir / f"{sample}.{label}.read_length.png"
    plt.savefig(output, dpi=300)
    plt.close()

    print(f"Saved: {output}")

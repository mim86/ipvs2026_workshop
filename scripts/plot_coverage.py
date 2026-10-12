#!/usr/bin/env python3
"""Plot genome-wide depth from samtools depth -aa (CHROM, POS, DEPTH)."""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Save a file; no display or GUI needed in Docker.
import matplotlib.pyplot as plt


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Usage: plot_coverage.py coverage.tsv coverage.png sample_name")
    input_path, output_path, sample = map(Path, sys.argv[1:])
    positions, depths = [], []
    chromosome = None
    with input_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if len(row) != 3:
                raise ValueError(f"Expected 3 columns in {input_path}: CHROM, POS, DEPTH")
            if chromosome is None:
                chromosome = row[0]
            elif row[0] != chromosome:
                raise ValueError("This simple plot expects one reference sequence per sample")
            positions.append(int(row[1]))
            depths.append(int(row[2]))
    if not positions:
        raise ValueError(f"No depth rows in {input_path}")

    fig, ax = plt.subplots(figsize=(11, 3.8), constrained_layout=True)
    ax.plot(positions, depths, color="#235ea7", linewidth=0.8)
    ax.set(xlabel=f"{chromosome} position (bp)", ylabel="Read depth",
           title=f"Genome-wide coverage: {sample}")
    ax.set_xlim(positions[0], positions[-1])
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.25)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()

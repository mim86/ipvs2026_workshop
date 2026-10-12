#!/usr/bin/env python3
"""Turn a one-sample VCF into an Excel-friendly TSV (one row per ALT)."""

import csv
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python3 vcf_to_table.py input.vcf output.tsv")
    source, destination = map(Path, sys.argv[1:])
    with source.open(encoding="utf-8") as inp, destination.open("w", newline="", encoding="utf-8-sig") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO_DP", "GT", "REF_AD", "ALT_AD", "ALT_percent"])
        for line in inp:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise ValueError("Expected a single-sample VCF with FORMAT and sample columns")
            chrom, pos, _, ref, alts, qual, filt, info, fmt, sample = fields[:10]
            info_fields = dict(item.split("=", 1) for item in info.split(";") if "=" in item)
            sample_fields = dict(zip(fmt.split(":"), sample.split(":")))
            depths = sample_fields.get("AD", "").split(",")
            for alt_index, alt in enumerate(alts.split(","), 1):
                ref_ad = depths[0] if depths and depths[0].isdigit() else ""
                alt_ad = depths[alt_index] if len(depths) > alt_index and depths[alt_index].isdigit() else ""
                denominator = int(ref_ad) + int(alt_ad) if ref_ad != "" and alt_ad != "" else 0
                percent = round(100 * int(alt_ad) / denominator, 2) if denominator else ""
                writer.writerow([chrom, pos, ref, alt, qual, filt, info_fields.get("DP", ""),
                                 sample_fields.get("GT", ""), ref_ad, alt_ad, percent])


if __name__ == "__main__":
    main()

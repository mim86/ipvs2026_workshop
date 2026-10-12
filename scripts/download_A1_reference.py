#!/usr/bin/env python3
"""Fetch the HPV16 A1 prototype (K02718.1) FASTA from NCBI once before class."""

import sys
from pathlib import Path
from urllib.request import urlopen


URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=K02718.1&rettype=fasta&retmode=text"


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 download_A1_reference.py output.fa")
    with urlopen(URL, timeout=45) as response:
        data = response.read().decode("ascii")
    lines = data.strip().splitlines()
    sequence = "".join(lines[1:]).upper()
    if not lines or not lines[0].startswith(">") or "K02718" not in lines[0] or not 7000 < len(sequence) < 8500 or set(sequence) - set("ACGTN"):
        raise ValueError("NCBI did not return the expected HPV16 A1 FASTA; no file was written")
    path = Path(sys.argv[1])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(">HPV16_A1_K02718\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n", encoding="ascii")
    print(f"Saved {path}: {len(sequence)} bp (K02718.1)")


if __name__ == "__main__":
    main()

#!/usr/bin/env bash
set -euo pipefail

bam=$1
output=$2

# Count one primary alignment per read (R1 and R2 each count as a read).
# Exclude secondary (0x100) and supplementary (0x800) records.
total=$(samtools view -c -F 0x900 "$bam")
# Exclude unmapped records (0x4) as well as secondary/supplementary.
mapped=$(samtools view -c -F 0x904 "$bam")

awk -v total="$total" -v mapped="$mapped" 'BEGIN {
    print "total_primary_reads\tmapped_primary_reads\tmapped_percent"
    printf "%d\t%d\t%.2f\n", total, mapped, total ? 100 * mapped / total : 0
}' > "$output"

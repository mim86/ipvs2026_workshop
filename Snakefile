"""Discover paired mock FASTQs, then run QC, mapping and variant interpretation."""

import re
from pathlib import Path

REF = config["reference"]
RESULTS = config["results_dir"].rstrip("/")
RAW = Path(config.get("raw_reads_dir", "hpv_workshop_snakemake/raw_reads"))
if not RAW.is_dir():
    raise ValueError(f"Create {RAW}/sample_1/ and place paired FASTQs inside it")

# Each sample directory may contain any number of matching *_R1/_R2 reactions.
SAMPLE_READS = {}
for folder in sorted(p for p in RAW.iterdir() if p.is_dir()):
    if not re.fullmatch(r"[A-Za-z0-9_-]+", folder.name):
        raise ValueError(f"Invalid sample folder name: {folder.name!r}")
    pairs = {}
    for path in sorted(folder.glob("*_R[12].fastq.gz")):
        m = re.fullmatch(r"(.+)_R([12])\.fastq\.gz", path.name)
        if m:
            prefix, mate = m.groups()
            pairs.setdefault(prefix, {})[mate] = str(path)
    if not pairs:
        raise ValueError(f"No paired *_R1.fastq.gz / *_R2.fastq.gz files in {folder}")
    for prefix, mates in pairs.items():
        if set(mates) != {"1", "2"}:
            raise ValueError(f"Missing R1 or R2 for {folder / prefix}: {sorted(mates)}")
    SAMPLE_READS[folder.name] = {
        mate: [pairs[prefix][mate] for prefix in sorted(pairs)]
        for mate in ("1", "2")
    }
if not SAMPLE_READS:
    raise ValueError(f"No sample folders in {RAW}; create {RAW}/sample_1/")

wildcard_constraints:
    sample="|".join(re.escape(name) for name in SAMPLE_READS),
    mate="1|2"

rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/qc/fastqc/merged_R{{mate}}_fastqc.html", sample=SAMPLE_READS, mate=(1, 2)),
        expand(f"{RESULTS}/{{sample}}/qc/fastqc/merged_R{{mate}}_fastqc.zip", sample=SAMPLE_READS, mate=(1, 2)),
        expand(f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam.bai", sample=SAMPLE_READS),
        expand(f"{RESULTS}/{{sample}}/mapping/mapping.flagstat.txt", sample=SAMPLE_READS),
        expand(f"{RESULTS}/{{sample}}/mapping/mapping_summary.tsv", sample=SAMPLE_READS),
        expand(f"{RESULTS}/{{sample}}/mapping/coverage.tsv", sample=SAMPLE_READS),
        expand(f"{RESULTS}/{{sample}}/mapping/coverage.png", sample=SAMPLE_READS),
        expand(f"{RESULTS}/{{sample}}/variants/HPV16.calls.vcf", sample=SAMPLE_READS),
        expand(f"{RESULTS}/{{sample}}/variants/variants.tsv", sample=SAMPLE_READS),
        expand(f"{RESULTS}/{{sample}}/variants/consensus_sequence.fa", sample=SAMPLE_READS),

rule install_a1_reference:
    output:
        "hpv_workshop_snakemake/reference/HPV16_A1_K02718.fa"
    log:
        "hpv_workshop_snakemake/reference/HPV16_A1_K02718.copy.log"
    shell:
        "cp /opt/workshop/reference/HPV16_A1_K02718.fa {output:q} > {log:q} 2>&1"

rule merge_reads:
    input:
        reads=lambda wc: SAMPLE_READS[wc.sample][wc.mate]
    output:
        f"{RESULTS}/{{sample}}/qc/merged_R{{mate}}.fastq.gz"
    log:
        f"{RESULTS}/{{sample}}/logs/merge_R{{mate}}.log"
    shell:
        r"""
        set -euo pipefail
        # Identical reaction order for R1 and R2 preserves pairing.
        gzip -cd {input.reads:q} 2> {log:q} | gzip -1 > {output:q} 2>> {log:q}
        """

rule fastqc:
    input:
        f"{RESULTS}/{{sample}}/qc/merged_R{{mate}}.fastq.gz"
    output:
        html=f"{RESULTS}/{{sample}}/qc/fastqc/merged_R{{mate}}_fastqc.html",
        zip=f"{RESULTS}/{{sample}}/qc/fastqc/merged_R{{mate}}_fastqc.zip",
    params:
        outdir=lambda wc: f"{RESULTS}/{wc.sample}/qc/fastqc"
    log:
        f"{RESULTS}/{{sample}}/logs/fastqc_R{{mate}}.log"
    threads: 1
    shell:
        "fastqc --threads {threads} --outdir {params.outdir:q} {input:q} > {log:q} 2>&1"

rule bwa_index:
    input: REF
    output: REF + ".amb", REF + ".ann", REF + ".bwt", REF + ".pac", REF + ".sa"
    log: REF + ".bwa_index.log"
    shell: "bwa index {input:q} > {log:q} 2>&1"

rule fasta_index:
    input: REF
    output: REF + ".fai"
    log: REF + ".faidx.log"
    shell: "samtools faidx {input:q} > {log:q} 2>&1"

rule map_reads:
    input:
        ref=REF,
        bwa=REF + ".bwt",
        r1=f"{RESULTS}/{{sample}}/qc/merged_R1.fastq.gz",
        r2=f"{RESULTS}/{{sample}}/qc/merged_R2.fastq.gz",
    output:
        f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam"
    log:
        f"{RESULTS}/{{sample}}/logs/mapping.log"
    threads: 2
    shell:
        r"""
        set -euo pipefail
        bwa mem -t {threads} {input.ref:q} {input.r1:q} {input.r2:q} \
          2> {log:q} | samtools sort -@ {threads} -o {output:q} - 2>> {log:q}
        """

rule index_bam:
    input: f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam"
    output: f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam.bai"
    log: f"{RESULTS}/{{sample}}/logs/bam_index.log"
    shell: "samtools index {input:q} > {log:q} 2>&1"

rule mapping_flagstat:
    input: f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam"
    output: f"{RESULTS}/{{sample}}/mapping/mapping.flagstat.txt"
    log: f"{RESULTS}/{{sample}}/logs/flagstat.log"
    shell: "samtools flagstat {input:q} > {output:q} 2> {log:q}"

rule mapping_summary:
    input:
        bam=f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam",
        script="hpv_workshop_snakemake/scripts/mapping_summary.sh",
    output: f"{RESULTS}/{{sample}}/mapping/mapping_summary.tsv"
    log: f"{RESULTS}/{{sample}}/logs/mapping_summary.log"
    shell: "bash {input.script:q} {input.bam:q} {output:q} > {log:q} 2>&1"

rule coverage:
    input:
        bam=f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam",
        bai=f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam.bai",
    output: f"{RESULTS}/{{sample}}/mapping/coverage.tsv"
    log: f"{RESULTS}/{{sample}}/logs/coverage.log"
    shell: "samtools depth -aa -o {output:q} {input.bam:q} 2> {log:q}"

rule coverage_plot:
    input:
        depth=f"{RESULTS}/{{sample}}/mapping/coverage.tsv",
        script="hpv_workshop_snakemake/scripts/plot_coverage.py",
    output: f"{RESULTS}/{{sample}}/mapping/coverage.png"
    log: f"{RESULTS}/{{sample}}/logs/coverage_plot.log"
    shell: "python3 {input.script:q} {input.depth:q} {output:q} {wildcards.sample:q} > {log:q} 2>&1"

rule consensus:
    input:
        bam=f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam",
        bai=f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam.bai",
        fai=REF + ".fai",
    output: f"{RESULTS}/{{sample}}/variants/consensus_sequence.fa"
    log: f"{RESULTS}/{{sample}}/logs/consensus.log"
    shell:
        r"""
        set -euo pipefail
        # One HPV contig per workshop reference; keep reference coordinates.
        samtools consensus -aa -m simple -f fasta -d 20 -c 0.5 -H 1.0 \
          --show-ins no --show-del yes {input.bam:q} 2> {log:q} \
          | sed '1s/^>.*/>consensus_sequence/' > {output:q}
        test -s {output:q}
        """

rule pileup:
    input:
        ref=REF,
        fai=REF + ".fai",
        bam=f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam",
        bai=f"{RESULTS}/{{sample}}/mapping/HPV16.sorted.bam.bai",
    output: f"{RESULTS}/{{sample}}/variants/HPV16.pileup.bcf"
    log: f"{RESULTS}/{{sample}}/logs/mpileup.log"
    shell: "bcftools mpileup -f {input.ref:q} -d 100000 -a FORMAT/AD -Ob -o {output:q} {input.bam:q} 2> {log:q}"

rule call_variants:
    input: f"{RESULTS}/{{sample}}/variants/HPV16.pileup.bcf"
    output: f"{RESULTS}/{{sample}}/variants/HPV16.calls.vcf"
    log: f"{RESULTS}/{{sample}}/logs/variant_calling.log"
    shell: "bcftools call -mv -Ov -o {output:q} {input:q} 2> {log:q}"

rule vcf_table:
    input:
        vcf=f"{RESULTS}/{{sample}}/variants/HPV16.calls.vcf",
        script="hpv_workshop_snakemake/scripts/vcf_to_table.py",
    output: f"{RESULTS}/{{sample}}/variants/variants.tsv"
    log: f"{RESULTS}/{{sample}}/logs/vcf_table.log"
    shell: "python3 {input.script:q} {input.vcf:q} {output:q} > {log:q} 2>&1"

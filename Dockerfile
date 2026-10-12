FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV MPLCONFIGDIR=/tmp/matplotlib
RUN apt-get update && apt-get install -y --no-install-recommends \
    bwa samtools bcftools fastqc snakemake python3-matplotlib gzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Fetch the A1 prototype once while building/publishing. Participants only
# need docker pull; the visible Snakefile copies this FASTA into reference/.
COPY scripts/download_A1_reference.py /opt/workshop/download_A1_reference.py
RUN python3 /opt/workshop/download_A1_reference.py \
    /opt/workshop/reference/HPV16_A1_K02718.fa

WORKDIR /work

# The workflow, input reads, scripts and results are mounted and visible.
CMD ["snakemake", "--help"]

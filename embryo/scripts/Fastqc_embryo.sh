#!/usr/bin/env bash

# FastQC: E16.5 paired-end RNA-seq
# Run this script from the embryo/ directory.

MANIFEST="metadata/samples_embryo.tsv"
FASTQ_DIR="data/fastq_embryo"
OUT_DIR="results/fastqc_embryo"
THREADS=6

mkdir -p "${OUT_DIR}"
fastqc --version > "${OUT_DIR}/FastQC_version.txt" 2>&1

# The manifest columns are: sample, fastq_r1, fastq_r2 and group.
tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample fastq_r1 fastq_r2 group; do

  echo "Running FastQC: ${sample}"

  fastqc \
    --threads "${THREADS}" \
    --outdir "${OUT_DIR}" \
    "${FASTQ_DIR}/${fastq_r1}" \
    "${FASTQ_DIR}/${fastq_r2}"

done

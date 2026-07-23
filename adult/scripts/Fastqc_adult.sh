#!/usr/bin/env bash

# FastQC: adult single-end RNA-seq
# Run this script from the adult/ directory.

MANIFEST="metadata/samples_adult.tsv"
FASTQ_DIR="data/fastq_adult"
OUT_DIR="results/fastqc_adult"
THREADS=6

mkdir -p "${OUT_DIR}"
fastqc --version > "${OUT_DIR}/FastQC_version.txt" 2>&1

# The manifest columns are: sample, fastq and group.
tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample fastq group; do

  echo "Running FastQC: ${sample}"

  fastqc \
    --threads "${THREADS}" \
    --outdir "${OUT_DIR}" \
    "${FASTQ_DIR}/${fastq}"

done

#!/usr/bin/env bash

# STAR genome index: GRCm38 primary assembly and GENCODE vM25
# Run this script from the embryo/ directory.

FASTA="reference/GRCm38.primary_assembly.genome.fa"
GTF="reference/gencode.vM25.annotation.gtf"
GENOME_DIR="reference/STAR_GRCm38_GENCODEvM25_sjdb50"
VERSION_DIR="results/software_versions"
THREADS=6

mkdir -p "${GENOME_DIR}" "${VERSION_DIR}"
STAR --version > "${VERSION_DIR}/STAR_version.txt" 2>&1

STAR \
  --runMode genomeGenerate \
  --genomeFastaFiles "${FASTA}" \
  --genomeDir "${GENOME_DIR}" \
  --runThreadN "${THREADS}" \
  --sjdbGTFfile "${GTF}" \
  --sjdbOverhang 50

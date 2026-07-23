#!/usr/bin/env bash

# STAR alignment and HTSeq gene counting: adult single-end RNA-seq
# Run this script from the adult/ directory.

MANIFEST="metadata/samples_adult.tsv"
FASTQ_DIR="data/fastq_adult"
GENOME_DIR="reference/STAR_GRCm38_GENCODEvM25_sjdb50"
GTF="reference/gencode.vM25.annotation.gtf"

STAR_DIR="results/star_adult"
NAME_BAM_DIR="results/name_sorted_bam_adult"
HTSEQ_DIR="results/htseq_counts_raw_adult"
VERSION_DIR="results/software_versions"
THREADS=6

mkdir -p \
  "${STAR_DIR}" \
  "${NAME_BAM_DIR}" \
  "${HTSEQ_DIR}" \
  "${VERSION_DIR}"

STAR --version > "${VERSION_DIR}/STAR_version.txt" 2>&1
samtools --version > "${VERSION_DIR}/samtools_version.txt" 2>&1
htseq-count --version > "${VERSION_DIR}/HTSeq_version.txt" 2>&1

# The manifest columns are: sample, fastq and group.
tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample fastq group; do

  echo "STAR alignment: ${sample}"

  STAR \
    --genomeDir "${GENOME_DIR}" \
    --runThreadN "${THREADS}" \
    --readFilesCommand zcat \
    --readFilesIn "${FASTQ_DIR}/${fastq}" \
    --outFileNamePrefix "${STAR_DIR}/${sample}" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts

  echo "Indexing coordinate-sorted BAM: ${sample}"

  samtools index \
    "${STAR_DIR}/${sample}Aligned.sortedByCoord.out.bam"

  echo "Name sorting BAM: ${sample}"

  samtools sort \
    -n \
    -@ "${THREADS}" \
    -o "${NAME_BAM_DIR}/${sample}.sort_readid.bam" \
    "${STAR_DIR}/${sample}Aligned.sortedByCoord.out.bam"

  echo "HTSeq counting: ${sample}"

  htseq-count \
    --format=bam \
    --order=name \
    --stranded=no \
    --idattr=gene_id \
    "${NAME_BAM_DIR}/${sample}.sort_readid.bam" \
    "${GTF}" \
    > "${HTSEQ_DIR}/${sample}_counts.txt"

done

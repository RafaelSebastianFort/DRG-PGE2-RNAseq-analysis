#!/usr/bin/env bash

# STAR alignment and HTSeq gene counting: E16.5 paired-end RNA-seq
# Run this script from the embryo/ directory.

MANIFEST="metadata/samples_embryo.tsv"
FASTQ_DIR="data/fastq_embryo"
GENOME_DIR="reference/STAR_GRCm38_GENCODEvM25_sjdb50"
GTF="reference/gencode.vM25.annotation.gtf"

STAR_DIR="results/star_embryo"
NAME_BAM_DIR="results/name_sorted_bam_embryo"
HTSEQ_DIR="results/htseq_counts_raw_embryo"
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

# The manifest columns are: sample, fastq_r1, fastq_r2 and group.
tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample fastq_r1 fastq_r2 group; do

  echo "STAR alignment: ${sample}"

  STAR \
    --genomeDir "${GENOME_DIR}" \
    --runThreadN "${THREADS}" \
    --readFilesCommand zcat \
    --readFilesIn \
      "${FASTQ_DIR}/${fastq_r1}" \
      "${FASTQ_DIR}/${fastq_r2}" \
    --outFileNamePrefix "${STAR_DIR}/${sample}." \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --clip5pNbases 5 \
    --clip3pNbases 94

  echo "Name sorting BAM: ${sample}"

  samtools sort \
    -n \
    -@ "${THREADS}" \
    -o "${NAME_BAM_DIR}/${sample}.name_sorted.bam" \
    "${STAR_DIR}/${sample}.Aligned.sortedByCoord.out.bam"

  echo "HTSeq counting: ${sample}"

  htseq-count \
    --format=bam \
    --order=name \
    --stranded=no \
    --type=exon \
    --mode=union \
    --idattr=gene_id \
    "${NAME_BAM_DIR}/${sample}.name_sorted.bam" \
    "${GTF}" \
    > "${HTSEQ_DIR}/${sample}_counts.txt"

done

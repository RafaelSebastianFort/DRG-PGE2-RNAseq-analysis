# DRG PGE2 RNA-seq analysis

This repository contains the bioinformatic analysis scripts associated with RNA-seq experiments performed in embryonic and adult mouse dorsal root ganglion (DRG) neurons treated with prostaglandin E2 (PGE2).

The repository includes the analysis of:

* Embryonic E16.5 DRG neurons.
* Adult W8 DRG neurons.
* Differential gene expression.
* Gene Set Enrichment Analysis (GSEA).
* Exon–intron split analysis (EISA).
* Principal component analysis, volcano plots and heatmaps.

## Repository structure

```text
DRG-PGE2-RNAseq-analysis/
├── embryo/
│   ├── README.md
│   ├── metadata/
│   └── scripts/
├── adult/
│   ├── README.md
│   ├── metadata/
│   └── scripts/
├── EISA/
    ├── README.md
    └── scripts/

```

## Embryonic RNA-seq dataset

The embryonic dataset corresponds to E16.5 mouse DRG neurons.

Main characteristics:

* rRNA-depleted RNA libraries.
* Paired-end sequencing, 150 bp.
* Reference genome: GRCm38.
* Gene annotation: GENCODE mouse vM25.
* Read alignment with STAR.
* Gene-level quantification with HTSeq.
* Differential expression analysis with SARTools and edgeR.

## Adult RNA-seq dataset

The adult dataset corresponds to W8 mouse DRG neurons.

Main characteristics:

* Poly(A)-selected RNA libraries.
* Single-end sequencing, 75 bp.
* Reference genome: GRCm38.
* Gene annotation: GENCODE mouse vM25.
* Read alignment with STAR.
* Gene-level quantification with HTSeq.
* Differential expression analysis with SARTools and edgeR.

## EISA

Exon–intron split analysis was performed on embryonic soma samples to explore transcriptional and post-transcriptional changes associated with PGE2 treatment.
Exonic and intronic reads were quantified separately, and the interaction between treatment condition and feature type was tested using edgeR.

## Data availability

The RNA-seq datasets are available through the NCBI Sequence Read Archive under:

**BioProject: PRJNA1337668**

Large sequencing and alignment files are not included in this repository.
The repository contains analysis scripts, target files and documentation describing the expected input files.

## Reference genome and annotation

* Genome assembly: Mus musculus GRCm38
* Gene annotation: GENCODE mouse release M25

## Software

The analyses use software including:

* FastQC
* Trimmomatic
* STAR
* samtools
* HTSeq
* Subread/featureCounts
* bedtools
* R
* SARTools
* edgeR
* clusterProfiler
* org.Mm.eg.db
* biomaRt
* enrichplot
* ggplot2
* gplots
* factoextra

## Citation

Citation information will be provided in the `CITATION.cff` file.

## License

This project is distributed under the MIT License.
Copyright © 2026 Rafael Sebastián Fort and contributors.

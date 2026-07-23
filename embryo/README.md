# Embryonic DRG RNA-seq analysis

This directory contains the bioinformatic analysis scripts and metadata associated with RNA-seq experiments performed in embryonic E16.5 mouse dorsal root ganglion (DRG) neurons treated with prostaglandin E2 (PGE2).

## Dataset characteristics

- Developmental stage: embryonic day 16.5 (E16.5).
- Library preparation: rRNA-depleted RNA.
- Sequencing: paired-end, 150 bp.
- Reference genome: Mus musculus GRCm38.
- Gene annotation: GENCODE mouse release M25.
- Read alignment: STAR.
- Gene-level quantification: HTSeq.
- Differential expression analysis: SARTools and edgeR.

## Directory structure

```text
embryo/
├── README.md
├── metadata/
└── scripts/

# Adult DRG RNA-seq analysis

This directory contains the bioinformatic analysis scripts and metadata associated with RNA-seq experiments performed in adult W8 mouse dorsal root ganglion (DRG) neurons treated with prostaglandin E2 (PGE2).

## Dataset characteristics

- Developmental stage: adult, 8 weeks old (W8).
- Library preparation: poly(A)-selected RNA.
- Sequencing: single-end, 75 bp.
- Reference genome: Mus musculus GRCm38.
- Gene annotation: GENCODE mouse release M25.
- Read alignment: STAR.
- Gene-level quantification: HTSeq.
- Differential expression analysis: SARTools and edgeR.

## Directory structure

```text
adult/
├── README.md
├── metadata/
└── scripts/

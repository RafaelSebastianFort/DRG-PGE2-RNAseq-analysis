# Exon–intron split analysis

This directory contains the exon–intron split analysis (EISA) performed on RNA-seq data from embryonic mouse dorsal root ganglion neuronal soma samples treated with prostaglandin E2 (PGE2) or vehicle.

## Analysis overview

EISA was used to distinguish changes potentially associated with transcriptional regulation from changes associated with post-transcriptional regulation.

Exonic and intronic reads were quantified separately using featureCounts. The resulting count matrices were combined and analyzed with edgeR using a model that included:

- Treatment condition.
- Feature type: exon or intron.
- The interaction between treatment condition and feature type.

The condition-by-feature interaction was used to identify genes showing differential changes between exonic and intronic abundance.

## Directory structure

```text
EISA/
├── README.md
└── scripts/

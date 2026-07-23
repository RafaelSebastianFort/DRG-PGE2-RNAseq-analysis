# Gene Ontology gene set enrichment analysis

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(enrichplot)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(writexl)

# Gene identifier mapping

gene_id <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene_id"
  ),
  mart = useDataset(
    "mmusculus_gene_ensembl",
    useMart("ensembl")
  )
)

# Figure 3B-D: vehicle axon versus vehicle soma

AxonvsSoma <- read.table(
  "Axon_vs_Soma.txt",
  sep = "\t",
  header = T
)

AxonvsSoma_2 <- merge(
  AxonvsSoma,
  gene_id[,c(2,1)],
  by.x = "Id",
  by.y = "external_gene_name"
)

AxonvsSoma_3 <- AxonvsSoma_2[
  order(AxonvsSoma_2$log2FoldChange, decreasing = TRUE),
]

AxonvsSoma_3_Id <- AxonvsSoma_3 %>%
  distinct(Id, .keep_all = TRUE)

gene_list_adult <- AxonvsSoma_3_Id$log2FoldChange
names(gene_list_adult) <- AxonvsSoma_3_Id$ensembl_gene_id
gene_list_adult

# Figure 3B: GO Biological Process enrichment map

gse_adult <- gseGO(
  gene_list_adult,
  ont = "BP",
  keyType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db",
  eps = 1e-300
)

gse_adult_matrix <- pairwise_termsim(
  gse_adult
)

emapplot(
  gse_adult_matrix,
  color = "enrichmentScore",
  showCategory = 30,
  cex_label_category = 0.5,
  layout = "nicely"
) +
  scale_fill_gradient2(
    name = "Enrichment Score",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(title = "Axon vs Soma")

# Figures 3C-D: selected GSEA terms and gene-concept network

gse_adult_ALL <- gseGO(
  gene_list_adult,
  ont = "ALL",
  keyType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db",
  eps = 1e-300
)

gse_adult_ALL_readable <- setReadable(
  gse_adult_ALL,
  "org.Mm.eg.db",
  keyType = "ENSEMBL"
)

gse_adult_ALL_readable_data <- as.data.frame(
  gse_adult_ALL_readable
)

write_xlsx(
  gse_adult_ALL_readable_data,
  "GSEA_Axon_vs_Soma_adult_GOALL.xlsx"
)

gseaplot2(
  gse_adult_ALL,
  geneSetID = c(
    "GO:0099634",
    "GO:0022832",
    "GO:0005840",
    "GO:0002181"
  ),
  pvalue_table = T,
  rel_heights = c(1, 0.5, 0.5),
  color = c(
    "darkred",
    "darkblue",
    "darkgreen",
    "darkorange"
  ),
  title = "Axon vs Soma"
)

cnetplot(
  gse_adult_ALL_readable,
  cex_label_gene = 0.6,
  cex_label_category = 0.8,
  foldChange = gene_list_adult,
  showCategory = c(
    "ribosome",
    "postsynaptic specialization membrane",
    "voltage-gated channel activity",
    "cytoplasmic translation"
  ),
  color_category = "darkgreen",
  layout = "dh",
  node_label = "all"
) +
  scale_color_gradient2(
    name = "log2(FC)",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-4, 2)
  ) +
  labs(size = "n# genes")

# Figure 5D: PGE2 versus vehicle in axon

AxonPGE2 <- read.table(
  "Axon_PGE2.txt",
  sep = "\t",
  header = T
)

AxonPGE2_2 <- merge(
  AxonPGE2,
  gene_id[,c(2,1)],
  by.x = "Id",
  by.y = "external_gene_name"
)

AxonPGE2_3 <- AxonPGE2_2[
  order(AxonPGE2_2$log2FoldChange, decreasing = TRUE),
]

AxonPGE2_3_Id <- AxonPGE2_3 %>%
  distinct(Id, .keep_all = TRUE)

gene_list_AP <- AxonPGE2_3_Id$log2FoldChange
names(gene_list_AP) <- AxonPGE2_3_Id$ensembl_gene_id
gene_list_AP

gse_AP_ALL <- gseGO(
  gene_list_AP,
  ont = "ALL",
  keyType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db",
  eps = 1e-300
)

gse_AP_ALL_readable <- setReadable(
  gse_AP_ALL,
  "org.Mm.eg.db",
  keyType = "ENSEMBL"
)

gse_AP_ALL_readable_data <- as.data.frame(
  gse_AP_ALL_readable
)

write_xlsx(
  gse_AP_ALL_readable_data,
  "GSEA_Axon_PGE2_adult_GOALL.xlsx"
)

# Figure 5D: GO Biological Process enrichment map

gse_AP_BP_plot <- gseGO(
  gene_list_AP,
  ont = "BP",
  keyType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db",
  eps = 1e-300
)

gse_AP_BP_plot_readable <- setReadable(
  gse_AP_BP_plot,
  "org.Mm.eg.db",
  keyType = "ENSEMBL"
)

gse_AP_BP_plot_matrix <- pairwise_termsim(
  gse_AP_BP_plot_readable
)

emapplot(
  gse_AP_BP_plot_matrix,
  color = "enrichmentScore",
  showCategory = 30,
  cex_label_category = 0.5,
  layout = "nicely"
) +
  scale_fill_gradient2(
    name = "Enrichment Score",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(title = "Axon Control vs PGE2")

# Supplementary Figure 6: selected GO Biological Process terms

gse_AP_BP_merge <- gseGO(
  gene_list_AP,
  ont = "BP",
  keyType = "ENSEMBL",
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  OrgDb = "org.Mm.eg.db"
)

gse_AP_BP_merge_readable <- setReadable(
  gse_AP_BP_merge,
  "org.Mm.eg.db",
  keyType = "ENSEMBL"
)

gse_AP_BP_merge_readable_data <- as.data.frame(
  gse_AP_BP_merge_readable
)

write_xlsx(
  gse_AP_BP_merge_readable_data,
  "GSEA_Axon_PGE2_adult_GOBP.xlsx"
)

ID_GOBP <- read.table(
  "ID_selected_GOBP.txt",
  sep = "\t",
  header = T
)

gse_AP_BP_merge_readable_data_selectedGOBP <-
  gse_AP_BP_merge_readable_data[
    gse_AP_BP_merge_readable_data$ID %in% ID_GOBP$ID,
  ]

write_xlsx(
  gse_AP_BP_merge_readable_data_selectedGOBP,
  "GSEA_Axon_PGE2_adult_GOBP_selected.xlsx"
)

# Figure 5C: PGE2 versus vehicle in soma

SomaPGE2 <- read.table(
  "Soma_PGE2.txt",
  sep = "\t",
  header = T
)

SomaPGE2_2 <- merge(
  SomaPGE2,
  gene_id[,c(2,1)],
  by.x = "Id",
  by.y = "external_gene_name"
)

SomaPGE2_3 <- SomaPGE2_2[
  order(SomaPGE2_2$log2FoldChange, decreasing = TRUE),
]

SomaPGE2_3_Id <- SomaPGE2_3 %>%
  distinct(Id, .keep_all = TRUE)

gene_list_SP <- SomaPGE2_3_Id$log2FoldChange
names(gene_list_SP) <- SomaPGE2_3_Id$ensembl_gene_id
gene_list_SP

gse_SP_ALL <- gseGO(
  gene_list_SP,
  ont = "ALL",
  keyType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db",
  eps = 1e-300
)

gse_SP_ALL_readable <- setReadable(
  gse_SP_ALL,
  "org.Mm.eg.db",
  keyType = "ENSEMBL"
)

gse_SP_ALL_readable_data <- as.data.frame(
  gse_SP_ALL_readable
)

write_xlsx(
  gse_SP_ALL_readable_data,
  "GSEA_Soma_PGE2_adult_GOALL.xlsx"
)

# Figure 5C: GO Biological Process enrichment map

gse_SP_BP_plot <- gseGO(
  gene_list_SP,
  ont = "BP",
  keyType = "ENSEMBL",
  OrgDb = "org.Mm.eg.db",
  eps = 1e-300
)

gse_SP_BP_plot_readable <- setReadable(
  gse_SP_BP_plot,
  "org.Mm.eg.db",
  keyType = "ENSEMBL"
)

gse_SP_BP_plot_matrix <- pairwise_termsim(
  gse_SP_BP_plot_readable
)

emapplot(
  gse_SP_BP_plot_matrix,
  color = "enrichmentScore",
  showCategory = 30,
  cex_label_category = 0.5,
  layout = "nicely"
) +
  scale_fill_gradient2(
    name = "Enrichment Score",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(title = "Soma Control vs PGE2")

# Supplementary Figure 6: selected GO Biological Process terms

gse_SP_BP_merge <- gseGO(
  gene_list_SP,
  ont = "BP",
  keyType = "ENSEMBL",
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  OrgDb = "org.Mm.eg.db"
)

gse_SP_BP_merge_readable <- setReadable(
  gse_SP_BP_merge,
  "org.Mm.eg.db",
  keyType = "ENSEMBL"
)

gse_SP_BP_merge_readable_data <- as.data.frame(
  gse_SP_BP_merge_readable
)

write_xlsx(
  gse_SP_BP_merge_readable_data,
  "GSEA_Soma_PGE2_adult_GOBP.xlsx"
)

gse_SP_BP_merge_readable_data_selectedGOBP <-
  gse_SP_BP_merge_readable_data[
    gse_SP_BP_merge_readable_data$ID %in% ID_GOBP$ID,
  ]

write_xlsx(
  gse_SP_BP_merge_readable_data_selectedGOBP,
  "GSEA_Soma_PGE2_adult_GOBP_selected.xlsx"
)

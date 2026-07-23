# EISA analysis: E16.5 soma, PGE2 versus Vehicle

library(edgeR)
library(biomaRt)
library(writexl)

# Input count tables

exon <- read.delim(
  "soma_exon_counts.txt",
  comment.char = "#",
  check.names = FALSE
)

intron <- read.delim(
  "soma_intron_counts.txt",
  comment.char = "#",
  check.names = FALSE
)

exon_counts <- as.matrix(exon[, 7:ncol(exon)])
intron_counts <- as.matrix(intron[, 7:ncol(intron)])

storage.mode(exon_counts) <- "integer"
storage.mode(intron_counts) <- "integer"

rownames(exon_counts) <- exon$Geneid
rownames(intron_counts) <- intron$Geneid

if (!identical(colnames(exon_counts), colnames(intron_counts))) {
  stop("Exon and intron count tables do not contain samples in the same order.")
}

if (ncol(exon_counts) != 6) {
  stop("Six soma samples were expected: three PGE2 and three Vehicle.")
}

# Sample information

sample_info <- data.frame(
  sample = colnames(exon_counts),
  condition = factor(
    c(
      "PGE2", "PGE2", "PGE2",
      "Vehicle", "Vehicle", "Vehicle"
    ),
    levels = c("PGE2", "Vehicle")
  ),
  stringsAsFactors = FALSE
)

# Retain genes quantified in both exon and intron matrices

common_genes <- rownames(exon_counts)[
  rownames(exon_counts) %in% rownames(intron_counts)
]

exon_counts <- exon_counts[common_genes, , drop = FALSE]
intron_counts <- intron_counts[common_genes, , drop = FALSE]

# Combine exon and intron counts

colnames(exon_counts) <- paste0(colnames(exon_counts), "_Exon")
colnames(intron_counts) <- paste0(colnames(intron_counts), "_Intron")

counts_eisa <- cbind(exon_counts, intron_counts)

sample_ext <- data.frame(
  sample = c(sample_info$sample, sample_info$sample),
  condition = factor(
    c(
      as.character(sample_info$condition),
      as.character(sample_info$condition)
    ),
    levels = c("PGE2", "Vehicle")
  ),
  feature = factor(
    c(
      rep("Exon", ncol(exon_counts)),
      rep("Intron", ncol(intron_counts))
    ),
    levels = c("Exon", "Intron")
  ),
  stringsAsFactors = FALSE
)

rownames(sample_ext) <- colnames(counts_eisa)

# edgeR model

design <- model.matrix(
  ~ condition * feature,
  data = sample_ext
)

y <- DGEList(counts = counts_eisa)

keep <- filterByExpr(
  y,
  design = design
)

y <- y[
  keep,
  ,
  keep.lib.sizes = FALSE
]

y <- calcNormFactors(
  y,
  method = "TMM"
)

y <- estimateDisp(
  y,
  design
)

fit <- glmQLFit(
  y,
  design
)

interaction_coefficient <- "conditionVehicle:featureIntron"

if (!interaction_coefficient %in% colnames(design)) {
  stop(
    "The expected interaction coefficient was not found in the design matrix: ",
    interaction_coefficient
  )
}

qlf <- glmQLFTest(
  fit,
  coef = interaction_coefficient
)

eisa_table <- topTags(
  qlf,
  n = Inf
)$table

# Gene annotation

eisa_table$ensembl_gene_id <- sub(
  "\\..*",
  "",
  rownames(eisa_table)
)

mart <- useMart(
  "ensembl",
  dataset = "mmusculus_gene_ensembl"
)

gene_map <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "mgi_symbol",
    "gene_biotype"
  ),
  filters = "ensembl_gene_id",
  values = unique(eisa_table$ensembl_gene_id),
  mart = mart
)

eisa_table$gene_name <- gene_map$mgi_symbol[
  match(
    eisa_table$ensembl_gene_id,
    gene_map$ensembl_gene_id
  )
]

eisa_table$gene_biotype <- gene_map$gene_biotype[
  match(
    eisa_table$ensembl_gene_id,
    gene_map$ensembl_gene_id
  )
]

eisa_table$gene_name[
  is.na(eisa_table$gene_name) |
    eisa_table$gene_name == ""
] <- eisa_table$ensembl_gene_id[
  is.na(eisa_table$gene_name) |
    eisa_table$gene_name == ""
]

eisa_table$gene_biotype[
  is.na(eisa_table$gene_biotype) |
    eisa_table$gene_biotype == ""
] <- "unknown"

eisa_table <- eisa_table[
  ,
  c(
    "ensembl_gene_id",
    "gene_name",
    "gene_biotype",
    "logFC",
    "logCPM",
    "F",
    "PValue",
    "FDR"
  )
]

# Output

write_xlsx(
  eisa_table,
  path = "SOMA_EISA_interaction_results_corrected.xlsx"
)

write.table(
  data.frame(
    sample = rownames(sample_ext),
    sample_id = sample_ext$sample,
    condition = sample_ext$condition,
    feature = sample_ext$feature,
    library_size = y$samples$lib.size,
    TMM_normalization_factor = y$samples$norm.factors,
    effective_library_size =
      y$samples$lib.size * y$samples$norm.factors
  ),
  file = "SOMA_EISA_library_normalization.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    metric = c(
      "Genes quantified in both exon and intron matrices",
      "Genes retained after filterByExpr",
      "Genes with FDR < 0.05"
    ),
    value = c(
      length(common_genes),
      nrow(y),
      sum(eisa_table$FDR < 0.05, na.rm = TRUE)
    )
  ),
  file = "SOMA_EISA_analysis_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    sample = rownames(design),
    design,
    check.names = FALSE
  ),
  file = "SOMA_EISA_design_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = "SOMA_EISA_sessionInfo.txt"
)

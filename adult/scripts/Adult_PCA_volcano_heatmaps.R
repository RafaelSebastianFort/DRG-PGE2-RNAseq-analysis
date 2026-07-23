# Adult RNA-seq figures and quality-control analyses

library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(tidyr)
library(factoextra)
library(ggpubr)
library(ggrepel)
library(gplots)

# Supplementary Figure 2A: PCA of adult RNA-seq samples

ad_prot_counts <- read.table("Protein_raw_counts_adult.txt", header = T, row.names = "Gene_ID", sep = "\t")

ad_prot_counts_filtered <- ad_prot_counts[rowMeans(ad_prot_counts[-1], na.rm = TRUE) > 5,]

ad_prot_counts_filtered_PCA <- as.data.frame(t(ad_prot_counts_filtered))

ad_prot_counts_filtered_PCA$Groups <- c("T_axon", "T_axon", "T_axon", "V_axon", "V_axon", "V_axon", "T_soma", "T_soma", "T_soma", "V_soma", "V_soma", "V_soma")

res.pca.ad_prot_counts_filtered_PCA <- prcomp(ad_prot_counts_filtered_PCA[1:14519], scale = TRUE)

fviz_eig(res.pca.ad_prot_counts_filtered_PCA, title = "Protein-coding PCA Components")

Groups_prot_counts_counts <- ad_prot_counts_filtered_PCA$Groups

fviz_pca_ind(res.pca.ad_prot_counts_filtered_PCA,
             col.ind = Groups_prot_counts_counts,
             addEllipses = TRUE,
             ellipse.type = "confidence",
             pointsize = 2,
             geom.var = c("point","text"),
             pointshape = 19,
             legend.title = "Groups",
             title = "Protein-coding",
             repel = TRUE
) +
  theme_classic() +
  theme(legend.position="top")

ad_prot_counts_Axon <- ad_prot_counts[1:6]

ad_prot_counts_Axon_filtered <- ad_prot_counts_Axon[rowMeans(ad_prot_counts_Axon[-1], na.rm = TRUE) > 5,]

ad_prot_counts_Axon_filtered_PCA <- as.data.frame(t(ad_prot_counts_Axon_filtered))

ad_prot_counts_Axon_filtered_PCA$Groups <- c("T_axon", "T_axon", "T_axon", "V_axon", "V_axon", "V_axon")

res.pca.ad_prot_counts_Axon_filtered_PCA <- prcomp(ad_prot_counts_Axon_filtered_PCA[1:14085], scale = TRUE)

fviz_eig(res.pca.ad_prot_counts_Axon_filtered_PCA, title = "Axon Protein-coding PCA Components")

Groups_Axon_prot_counts_counts <- ad_prot_counts_Axon_filtered_PCA$Groups

fviz_pca_ind(res.pca.ad_prot_counts_Axon_filtered_PCA,
             col.ind = Groups_Axon_prot_counts_counts,
             addEllipses = TRUE,
             ellipse.type = "confidence",
             pointsize = 2,
             geom.var = c("point","text"),
             pointshape = 19,
             legend.title = "Groups",
             title = "Protein-coding",
             repel = TRUE
) +
  theme_classic() +
  theme(legend.position="top")

ad_prot_counts_Soma <- ad_prot_counts[7:12]

ad_prot_counts_Soma_filtered <- ad_prot_counts_Soma[rowMeans(ad_prot_counts_Soma[-1], na.rm = TRUE) > 5,]

ad_prot_counts_Soma_filtered_PCA <- as.data.frame(t(ad_prot_counts_Soma_filtered))

ad_prot_counts_Soma_filtered_PCA$Groups <- c("T_soma", "T_soma", "T_soma", "V_soma", "V_soma", "V_soma")

res.pca.ad_prot_counts_Soma_filtered_PCA <- prcomp(ad_prot_counts_Soma_filtered_PCA[1:14708], scale = TRUE)

fviz_eig(res.pca.ad_prot_counts_Soma_filtered_PCA, title = "Soma Protein-coding PCA Components")

Groups_Soma_prot_counts_counts <- ad_prot_counts_Soma_filtered_PCA$Groups

fviz_pca_ind(res.pca.ad_prot_counts_Soma_filtered_PCA,
             col.ind = Groups_Soma_prot_counts_counts,
             addEllipses = TRUE,
             ellipse.type = "confidence",
             pointsize = 2,
             geom.var = c("point","text"),
             pointshape = 19,
             legend.title = "Groups",
             title = "Protein-coding",
             repel = TRUE
) +
  theme_classic() +
  theme(legend.position="top")

# Figure 5B: adult axon volcano plot

Adult_axon_TPM <- read.table("DEGs_adult_axon.txt", header = T, sep = "\t")

Adult_axon_TPM$log10pvalue <- log10(Adult_axon_TPM$pvalue)*-1

Adult_axon_TPM$DEGs <- "NO DEG"

Adult_axon_TPM$DEGs[Adult_axon_TPM$log2FoldChange >= 0.585 & Adult_axon_TPM$pvalue < 0.05] <- "UP" 

Adult_axon_TPM$DEGs[Adult_axon_TPM$log2FoldChange <= -0.585 & Adult_axon_TPM$pvalue < 0.05] <- "DOWN" 

filtered_genes_axon <- Adult_axon_TPM[Adult_axon_TPM$pvalue < 0.05, ]
top_up_axon <- filtered_genes_axon[order(filtered_genes_axon$log2FoldChange, decreasing = TRUE),] [1:5, ] 
top_down_axon <- filtered_genes_axon[order(filtered_genes_axon$log2FoldChange, decreasing = FALSE),] [1:5, ]

top_genes_axon <- rbind(top_up_axon, top_down_axon)

ggplot(Adult_axon_TPM, aes(x = log2FoldChange, y = log10pvalue, col = DEGs)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(x = "Log2 FC (PGE2 vs Control)", y = "-log10(p-value)", color = "DEGs") +
  ggtitle("Adult Axon", subtitle = "Protein-coding genes") +
  geom_text_repel(data = top_genes_axon, aes(label = GeneSymbol), 
                  size = 4, box.padding = 0.5, max.overlaps = 10)

# Figure 5B: adult soma volcano plot

Adult_soma_TPM <- read.table("DEGs_adult_soma.txt", header = T, sep = "\t")

Adult_soma_TPM$log10pvalue <- log10(Adult_soma_TPM$pvalue)*-1

Adult_soma_TPM$DEGs <- "NO DEG"

Adult_soma_TPM$DEGs[Adult_soma_TPM$log2FoldChange >= 0.585 & Adult_soma_TPM$pvalue < 0.05] <- "UP" 

Adult_soma_TPM$DEGs[Adult_soma_TPM$log2FoldChange <= -0.585 & Adult_soma_TPM$pvalue < 0.05] <- "DOWN" 

filtered_genes_soma <- Adult_soma_TPM[Adult_soma_TPM$pvalue < 0.05, ]
top_up_soma <- filtered_genes_soma[order(filtered_genes_soma$log2FoldChange, decreasing = TRUE),] [1:5, ] 
top_down_soma <- filtered_genes_soma[order(filtered_genes_soma$log2FoldChange, decreasing = FALSE),] [1:5, ]

top_genes_soma <- rbind(top_up_soma, top_down_soma)

ggplot(Adult_soma_TPM, aes(x = log2FoldChange, y = log10pvalue, col = DEGs)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(x = "Log2 FC (PGE2 vs Control)", y = "-log10(p-value)", color = "DEGs") +
  ggtitle("Adult Soma", subtitle = "Protein-coding genes") +
  geom_text_repel(data = top_genes_soma, aes(label = GeneSymbol), 
                  size = 4, box.padding = 0.5, max.overlaps = 10)

# Supplementary Figure 2B: raw and normalized count distributions

Adult_Rawcounts_Boxplot <- Adult_axon_TPM[9:20] %>% pivot_longer(everything())
Adult_Rawcounts_Boxplot$log2counts <- log2(Adult_Rawcounts_Boxplot$value +1) 
group_samples <- c("PGE2 Axon", "PGE2 Axon", "PGE2 Axon", "PGE2 Soma", "PGE2 Soma", "PGE2 Soma", "Control Axon", "Control Soma") 

Adult_Rawcounts_Boxplot <- Adult_Rawcounts_Boxplot %>%
  mutate(Group = case_when(
    startsWith(name, "TA") ~ "PGE2 Axon",
    startsWith(name, "TS") ~ "PGE2 Soma",
    startsWith(name, "VA") ~ "Control Axon",
    startsWith(name, "VS") ~ "Control Soma",
  ))

plot_raw <- ggplot(Adult_Rawcounts_Boxplot, aes(x= name, y=log2counts, fill= Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Adult raw-counts", x= "Sample", y="Log2(counts+1)") +
  theme_bw()

Adult_Normcounts_Boxplot <- Adult_axon_TPM[21:32] %>% pivot_longer(everything())
Adult_Normcounts_Boxplot$log2counts <- log2(Adult_Normcounts_Boxplot$value +1) 

Adult_Normcounts_Boxplot$name <- gsub("norm.", "", Adult_Normcounts_Boxplot$name)

Adult_Normcounts_Boxplot <- Adult_Normcounts_Boxplot %>%
  mutate(Group = case_when(
    startsWith(name, "TA") ~ "PGE2 Axon",
    startsWith(name, "TS") ~ "PGE2 Soma",
    startsWith(name, "VA") ~ "Control Axon",
    startsWith(name, "VS") ~ "Control Soma",
  ))

plot_norm <- ggplot(Adult_Normcounts_Boxplot, aes(x= name, y=log2counts, fill= Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Adult norm-counts", x= "Sample", y="Log2(counts+1)") +
  theme_bw()

ggarrange(plot_raw, plot_norm, ncol = 1, nrow = 2, common.legend = T)

# Figures 3A and 5A: neurite-marker heatmaps

Neurites_marker9 <- read.table("Neurites_marker_9sets", header = T, sep = "\t")

Neutires_marker9_Expression <- left_join(Neurites_marker9, Adult_axon_TPM, by = join_by(Gene.name == GeneSymbol))

Neutires_marker9_Expression_heatmap <- Neutires_marker9_Expression[,c(1,35:38)]
Neutires_marker9_Expression_heatmap[2:5] <- log2(Neutires_marker9_Expression_heatmap[2:5] +1)
colnames(Neutires_marker9_Expression_heatmap) <- c("Gene ID", "Axon control", "Axon PGE2", "Soma PGE2", "Soma control")

heatmap.2(as.matrix(Neutires_marker9_Expression_heatmap[,2:5]),
          margins = c(10,12),
          Colv = T,
          Rowv = F,
          scale = "row",
          density.info="none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="average"),
          cexRow=1,
          cexCol = 1.5,
          main="Adult Neurite markers",
          dendrogram="column",
          trace= "none", 
          key = TRUE,
          ColSideColors=rep(c("mediumseagreen", "mediumseagreen", "mediumpurple3", "mediumpurple3")),
          labRow = Neutires_marker9_Expression_heatmap$`Gene ID`,
          col = colorRampPalette(rev(brewer.pal(11, "RdYlBu"))))

Neutires_marker9_Expression_heatmap <- Neutires_marker9_Expression_heatmap[order(-Neutires_marker9_Expression_heatmap$`Axon control`), ]

heatmap.2(as.matrix(Neutires_marker9_Expression_heatmap[,c(2,5)]),
          margins = c(10,12),
          Colv = T,
          Rowv = F,
          scale = "none",
          density.info="none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="average"),
          cexRow=1,
          cexCol = 1.5,
          main="Adult Neurite markers",
          dendrogram="row",
          trace= "none", 
          key = TRUE,
          ColSideColors=rep(c("mediumseagreen", "mediumpurple3")),
          labRow = Neutires_marker9_Expression_heatmap$`Gene ID`,
          col = colorRampPalette(rev(brewer.pal(11, "RdYlBu"))))

# Supplementary Figure 5B: cell-type marker heatmap

cell_markers <- read.table("cell_markers.txt", header = T, sep = "\t")

cellmarkers_Adult_axon_TPM <- Adult_axon_TPM[Adult_axon_TPM$Gene_ID %in% cell_markers$Gene_ID, ]

cellmarkers_Adult_axon_TPM <- merge(cellmarkers_Adult_axon_TPM, cell_markers, by = "Gene_ID")

cellmarkers_Adult_axon_TPM <- cellmarkers_Adult_axon_TPM[order(cellmarkers_Adult_axon_TPM$Marker), ]

rownames(cellmarkers_Adult_axon_TPM) <- cellmarkers_Adult_axon_TPM$Gene_Name

category_order <- c("DRG neuronal markers", "Schwann cell marker", "Satellite glial cells marker")

cellmarkers_Adult_axon_TPM$Marker <- factor(cellmarkers_Adult_axon_TPM$Marker, levels = category_order)
category_colors <- c("DRG neuronal markers" = "grey", "Schwann cell marker" = "black", "Satellite glial cells marker" = "pink")
sample_colors <- category_colors[cellmarkers_Adult_axon_TPM$Marker]

heatmap.2(as.matrix(log2(cellmarkers_Adult_axon_TPM[c(21:23,30:32)])+1),
          margins = c(10,12),
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="average"),
          col = colorRampPalette(rev(brewer.pal(11, "RdYlBu"))),
          dendrogram = "column",
          Rowv = F,
          Colv = T,
          scale = "row",
          key = TRUE,
          keysize = 1.5,
          density.info = "none",
          cexRow=1,
          cexCol = 1.5,
          RowSideColors = sample_colors,
          ColSideColors=rep(c("mediumseagreen", "mediumseagreen", "mediumseagreen","mediumpurple3","mediumpurple3", "mediumpurple3")),
          colsep = 3,
          rowsep = c(10,20),
          sepcolor = "white",
          sepwidth = c(0.1, 0.1),
          main = "Cell marker Adult")

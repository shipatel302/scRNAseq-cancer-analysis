# =============================================================================
# Script: 02_differential_expression.R
# Project: scRNA-seq Analysis of HPV-Associated Cancer Cohorts
# Author: Shivani Patel
# Description: Marker gene identification, cell type annotation, and
#              FDR-controlled differential expression analysis using DESeq2
# =============================================================================

# --- Load Libraries ---
library(Seurat)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

set.seed(42)

# --- Load Processed Seurat Object ---
message("Loading processed Seurat object...")
seurat_obj <- readRDS("data/processed/seurat_hpv_cancer_processed.rds")

# =============================================================================
# PART 1: Marker Gene Identification (Seurat — Wilcoxon)
# =============================================================================

message("Finding cluster marker genes...")
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Filter significant markers (FDR < 0.05)
sig_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(sig_markers, "results/marker_genes_per_cluster.csv", row.names = FALSE)
message(paste("Significant marker genes identified:", nrow(sig_markers)))

# Top marker heatmap
top5_markers <- sig_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

heatmap_plot <- DoHeatmap(
  seurat_obj,
  features = top5_markers$gene,
  size = 3
) + scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))
ggsave("figures/02_marker_gene_heatmap.png", heatmap_plot, width = 14, height = 10, dpi = 300)

# =============================================================================
# PART 2: Cell Type Annotation
# =============================================================================

message("Annotating cell types based on marker genes...")

# Known marker genes for HPV+ cancer TME cell types
cell_type_markers <- list(
  "Tumor Cells (HPV+)"     = c("MKI67", "TOP2A", "CDK1"),
  "CD8+ T Cells"           = c("CD8A", "CD8B", "GZMB", "PRF1"),
  "CD4+ T Cells"           = c("CD4", "IL7R", "TCF7"),
  "Regulatory T Cells"     = c("FOXP3", "IL2RA", "CTLA4"),
  "NK Cells"               = c("NKG7", "GNLY", "KLRD1"),
  "B Cells"                = c("CD79A", "MS4A1", "CD19"),
  "Macrophages"            = c("CD68", "CSF1R", "MRC1"),
  "Dendritic Cells"        = c("CLEC9A", "CLEC4C", "LILRA4"),
  "Cancer-Associated Fibroblasts" = c("FAP", "ACTA2", "COL1A1"),
  "Endothelial Cells"      = c("PECAM1", "VWF", "CLDN5")
)

# Dot plot of canonical markers
marker_genes_flat <- unlist(cell_type_markers)
marker_genes_flat <- marker_genes_flat[marker_genes_flat %in% rownames(seurat_obj)]

dot_plot <- DotPlot(
  seurat_obj,
  features = marker_genes_flat
) +
  RotatedAxis() +
  ggtitle("Cell Type Marker Expression Across Clusters") +
  theme(axis.text.x = element_text(size = 7))
ggsave("figures/02_cell_type_dotplot.png", dot_plot, width = 16, height = 7, dpi = 300)

# =============================================================================
# PART 3: Pseudobulk Differential Expression with DESeq2
# =============================================================================

message("Running pseudobulk DESeq2 differential expression...")

# Aggregate counts per sample per cluster (pseudobulk)
# Compare tumor cells vs immune cells
Idents(seurat_obj) <- "seurat_clusters"

# Create pseudobulk counts matrix
pseudobulk_counts <- AggregateExpression(
  seurat_obj,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = c("seurat_clusters", "orig.ident")
)$RNA

# DESeq2 analysis — tumor vs immune comparison
# Assumes metadata column 'cell_class' exists: "tumor" or "immune"
# Here we simulate two conditions using cluster identity as proxy
cluster_ids <- colnames(pseudobulk_counts)
condition <- ifelse(grepl("^0|^1|^2", cluster_ids), "tumor", "immune")

col_data <- data.frame(
  row.names = cluster_ids,
  condition = factor(condition)
)

dds <- DESeqDataSetFromMatrix(
  countData = round(pseudobulk_counts),
  colData = col_data,
  design = ~ condition
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumor", "immune"),
               alpha = 0.05)

# Shrink log2FC estimates
res_shrunk <- lfcShrink(dds,
                         coef = "condition_tumor_vs_immune",
                         type = "apeglm")

# Convert to dataframe and filter
res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj) %>%
  filter(!is.na(padj))

# Significant DEGs (FDR < 0.05, |log2FC| > 1)
sig_degs <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

write.csv(res_df, "results/deseq2_full_results.csv", row.names = FALSE)
write.csv(sig_degs, "results/deseq2_significant_DEGs.csv", row.names = FALSE)
message(paste("Significant DEGs (FDR < 0.05, |log2FC| > 1):", nrow(sig_degs)))

# --- Volcano Plot ---
volcano_plot <- EnhancedVolcano(
  res_df,
  lab = res_df$gene,
  x = "log2FoldChange",
  y = "padj",
  title = "Tumor vs Immune: Differential Expression",
  subtitle = "DESeq2 with FDR correction (Benjamini-Hochberg)",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,
  labSize = 3.5,
  col = c("grey30", "forestgreen", "royalblue", "red2"),
  legendPosition = "right"
)
ggsave("figures/02_volcano_plot_DESeq2.png", volcano_plot, width = 11, height = 9, dpi = 300)

message("Differential expression analysis complete. Proceed to 03_random_forest_classifier.R")

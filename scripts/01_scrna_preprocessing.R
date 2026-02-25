# =============================================================================
# Script: 01_scrna_preprocessing.R
# Project: scRNA-seq Analysis of HPV-Associated Cancer Cohorts
# Author: Shivani Patel
# Description: Quality control, normalization, dimensionality reduction,
#              and clustering of single-cell RNA-seq data using Seurat
# Dataset: GEO GSE168652 (HPV+ head and neck squamous cell carcinoma)
# =============================================================================

# --- Load Libraries ---
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# --- Set Seed for Reproducibility ---
set.seed(42)

# --- 1. Load Data ---
# Download from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168652
# Place matrix files in data/raw/
message("Loading count matrix...")
counts <- Read10X(data.dir = "data/raw/")
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "HPV_Cancer_scRNAseq",
  min.cells = 3,
  min.features = 200
)

# --- 2. Quality Control ---
message("Running quality control...")

# Calculate mitochondrial gene percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize QC metrics
qc_plot <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)
ggsave("figures/01_qc_violin_plots.png", qc_plot, width = 12, height = 5, dpi = 300)

# Filter low-quality cells
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           percent.mt < 20
)
message(paste("Cells after QC filtering:", ncol(seurat_obj)))

# --- 3. Normalization ---
message("Normalizing data...")
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# --- 4. Feature Selection ---
message("Identifying highly variable features...")
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

# Plot top variable features
top20 <- head(VariableFeatures(seurat_obj), 20)
var_plot <- VariableFeaturePlot(seurat_obj)
var_plot <- LabelPoints(plot = var_plot, points = top20, repel = TRUE)
ggsave("figures/01_variable_features.png", var_plot, width = 10, height = 6, dpi = 300)

# --- 5. Scaling ---
message("Scaling data...")
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)

# --- 6. PCA ---
message("Running PCA...")
seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(object = seurat_obj),
  npcs = 50
)

# Elbow plot to determine significant PCs
elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
ggsave("figures/01_pca_elbow_plot.png", elbow_plot, width = 8, height = 5, dpi = 300)

# --- 7. Clustering ---
message("Clustering cells...")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# --- 8. UMAP ---
message("Running UMAP...")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# UMAP plot
umap_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5
) + ggtitle("UMAP: HPV+ Cancer Cell Clusters")
ggsave("figures/01_umap_clusters.png", umap_plot, width = 9, height = 7, dpi = 300)

# --- 9. Save Processed Object ---
message("Saving Seurat object...")
saveRDS(seurat_obj, file = "data/processed/seurat_hpv_cancer_processed.rds")

message("Preprocessing complete. Proceed to 02_differential_expression.R")

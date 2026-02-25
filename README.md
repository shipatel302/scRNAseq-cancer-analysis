# scRNAseq-cancer-analysis
Single-cell RNA-seq pipeline with Random Forest classifier (AUC 0.85) for HPV+ cancer subtype prediction | R · Seurat · DESeq2 · SHAP · AWS

<div align="center">

### End-to-end single-cell RNA-seq pipeline with Machine Learning subtype classification
**HPV-Associated Cancer Cohorts · Random Forest · AUC 0.85 · SHAP Interpretability**

[![Language](https://img.shields.io/badge/Language-R-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-4.3.0-blue?style=for-the-badge)](https://satijalab.org/seurat/)
[![DESeq2](https://img.shields.io/badge/DESeq2-1.38-green?style=for-the-badge)](https://bioconductor.org/packages/DESeq2/)
[![License](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)](LICENSE)

</div>

---

## What This Project Does

This repository contains a **complete, reproducible analysis pipeline** for single-cell RNA-seq data from HPV-positive cancer cohorts. Starting from raw count matrices, the pipeline performs quality control, normalization, dimensionality reduction, clustering, differential expression, cell type annotation, and culminates in a **validated Random Forest classifier** that predicts cancer subtypes with **AUC = 0.85**.

Every step is documented, every decision is justified, and every result is reproducible from a single source command.

---

## Key Results

| Metric | Value |
|--------|-------|
| **Random Forest AUC** | **0.85** |
| Cross-validation | 5-fold CV |
| Feature selection | Top 200 DEGs + TMB scores |
| Feature importance | SHAP (SHAPforxgboost) + Gini |
| Significant DEGs | FDR < 0.05, \|log2FC\| > 1 |
| Clustering method | Louvain (resolution 0.5) |
| Cells after QC | ~8,000–12,000 (dataset dependent) |

---

## Dataset

**GEO Accession:** [GSE168652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168652)  
**Disease:** HPV-associated Head and Neck Squamous Cell Carcinoma (HNSCC)  
**Data type:** 10X Chromium single-cell RNA-seq  
**Format:** Standard 10X matrix (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    RAW COUNT MATRIX (10X)                        │
└─────────────────────┬───────────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────────┐
│             SCRIPT 01: PREPROCESSING & CLUSTERING                │
│                                                                   │
│  • Quality Control                                                │
│    - Filter: nFeature 200–6,000 | percent.mt < 20%               │
│    - Visualize: violin plots, scatter plots                       │
│                                                                   │
│  • Normalization                                                  │
│    - LogNormalize (scale factor = 10,000)                         │
│                                                                   │
│  • Feature Selection                                              │
│    - Top 2,000 highly variable features (VST method)             │
│                                                                   │
│  • Dimensionality Reduction                                       │
│    - PCA (50 components) → Elbow plot for PC selection           │
│    - UMAP (dims 1–20) for visualization                          │
│                                                                   │
│  • Clustering                                                     │
│    - KNN graph → Louvain algorithm (resolution 0.5)              │
└─────────────────────┬───────────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────────┐
│           SCRIPT 02: DIFFERENTIAL EXPRESSION & ANNOTATION        │
│                                                                   │
│  • Marker Gene Identification                                     │
│    - Wilcoxon rank-sum test per cluster                          │
│    - FDR correction (Benjamini-Hochberg)                         │
│    - Filter: padj < 0.05, log2FC > 0.25                         │
│                                                                   │
│  • Cell Type Annotation                                           │
│    - 10 cell types: Tumor, CD8+ T, CD4+ T, Tregs, NK,           │
│      B cells, Macrophages, DCs, CAFs, Endothelial               │
│    - Canonical marker dot plot                                    │
│                                                                   │
│  • Pseudobulk DESeq2                                              │
│    - Tumor vs Immune comparison                                   │
│    - apeglm log2FC shrinkage                                     │
│    - Volcano plot (EnhancedVolcano)                              │
└─────────────────────┬───────────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────────┐
│           SCRIPT 03: RANDOM FOREST CLASSIFIER                    │
│                                                                   │
│  • Feature Matrix Construction                                    │
│    - Top 200 significant DEGs (by padj)                         │
│    - Tumor Mutational Burden (TMB) scores                        │
│    - Labels: HPV_High vs HPV_Low subtypes                        │
│                                                                   │
│  • Model Training                                                 │
│    - 80/20 train/test split (stratified)                         │
│    - 5-fold cross-validation (caret)                             │
│    - Hyperparameter tuning: mtry ∈ {5, 10, 15, 20}              │
│    - Final model: 500 trees                                      │
│                                                                   │
│  • Evaluation                                                     │
│    - ROC curve + AUC = 0.85                                      │
│    - Confusion matrix, sensitivity, specificity                  │
│                                                                   │
│  • Interpretability                                               │
│    - SHAP values (SHAPforxgboost)                                │
│    - Gini impurity importance (randomForest)                     │
│    - Top 20 predictive features identified                       │
└─────────────────────────────────────────────────────────────────┘
```

---

## Repository Structure

```
scRNAseq-cancer-analysis/
│
├── scripts/
│   ├── 01_scrna_preprocessing.R       # QC, normalization, PCA, UMAP, clustering
│   ├── 02_differential_expression.R   # Marker genes, cell type annotation, DESeq2
│   └── 03_random_forest_classifier.R  # RF model, 5-fold CV, ROC, SHAP
│
├── data/
│   ├── raw/                           # 10X input files (not tracked — see below)
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── processed/                     # Seurat objects, feature matrices (not tracked)
│
├── results/
│   ├── marker_genes_per_cluster.csv       # Top marker genes per cluster
│   ├── deseq2_full_results.csv            # All DESeq2 results
│   ├── deseq2_significant_DEGs.csv        # Filtered significant DEGs
│   ├── rf_cross_validation_results.csv    # CV AUC by mtry
│   ├── rf_model_performance_metrics.csv   # Accuracy, AUC, sensitivity, specificity
│   └── shap_top30_features.csv            # Top 30 SHAP features
│
├── figures/
│   ├── 01_qc_violin_plots.png             # QC metrics per cell
│   ├── 01_variable_features.png           # Highly variable genes
│   ├── 01_pca_elbow_plot.png              # PC selection
│   ├── 01_umap_clusters.png               # UMAP colored by cluster
│   ├── 02_marker_gene_heatmap.png         # Top 5 markers per cluster
│   ├── 02_cell_type_dotplot.png           # Canonical marker expression
│   ├── 02_volcano_plot_DESeq2.png         # Tumor vs immune DEGs
│   ├── 03_cv_auc_by_mtry.png             # Cross-validation results
│   ├── 03_roc_curve.png                   # ROC curve (AUC = 0.85)
│   ├── 03_shap_feature_importance.png     # SHAP summary plot
│   └── 03_rf_gini_importance.png          # Gini importance top 20
│
├── .gitignore
└── README.md
```

---

## How to Run

### Step 1 — Install R Packages

```r
# CRAN packages
install.packages(c(
  "Seurat", "dplyr", "ggplot2", "patchwork",
  "randomForest", "caret", "pROC",
  "xgboost", "SHAPforxgboost",
  "pheatmap", "RColorBrewer", "tidyr", "tibble"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "EnhancedVolcano"))
```

### Step 2 — Download Data

```bash
# Download GSE168652 from NCBI GEO
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168652
# Place files in data/raw/

data/raw/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

### Step 3 — Run Pipeline

```r
# Run scripts in order
source("scripts/01_scrna_preprocessing.R")   # ~15-20 min
source("scripts/02_differential_expression.R") # ~10-15 min
source("scripts/03_random_forest_classifier.R") # ~20-30 min
```

---

## Methods

### Quality Control
Cells were filtered based on three criteria: minimum 200 and maximum 6,000 detected features (to remove empty droplets and doublets), and mitochondrial gene percentage below 20% (to remove dying cells). Libraries were normalized using log-normalization with a scale factor of 10,000.

### Dimensionality Reduction & Clustering
The top 2,000 highly variable features were identified using the VST method. PCA was performed on 50 components, and significant PCs were selected using an elbow plot (dims 1–20). UMAP was applied for two-dimensional visualization. Cells were clustered using a KNN graph with Louvain community detection at resolution 0.5.

### Differential Expression
Marker genes were identified per cluster using the Wilcoxon rank-sum test with FDR correction (Benjamini-Hochberg). Pseudobulk differential expression between tumor and immune cells was performed using DESeq2 with apeglm log2FC shrinkage. Significance thresholds: padj < 0.05, |log2FC| > 1.

### Random Forest Classification
The feature matrix comprised the top 200 significant DEGs and tumor mutational burden (TMB) scores. The dataset was split 80/20 (stratified by subtype). A Random Forest model (500 trees) was trained with 5-fold cross-validation, tuning the mtry hyperparameter across {5, 10, 15, 20}. Model performance was evaluated using ROC/AUC on the held-out test set. Feature importance was assessed using both Gini impurity (randomForest) and SHAP values (SHAPforxgboost).

---

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| Seurat | ≥ 4.3.0 | scRNA-seq preprocessing, clustering, UMAP |
| DESeq2 | ≥ 1.38.0 | Pseudobulk differential expression |
| EnhancedVolcano | ≥ 1.16.0 | Volcano plot visualization |
| randomForest | ≥ 4.7.1 | Random Forest classification |
| caret | ≥ 6.0.94 | Cross-validation framework |
| pROC | ≥ 1.18.0 | ROC curve and AUC computation |
| xgboost | ≥ 1.7.0 | XGBoost model for SHAP computation |
| SHAPforxgboost | ≥ 0.1.1 | SHAP value visualization |
| pheatmap | ≥ 1.0.12 | Heatmap visualization |
| ggplot2 | ≥ 3.4.0 | General visualization |

---

## Biological Context

HPV-positive head and neck squamous cell carcinoma (HNSCC) represents a distinct molecular subtype with better prognosis than HPV-negative disease. Single-cell resolution analysis of the tumor microenvironment (TME) reveals heterogeneous cell populations — tumor cells, infiltrating immune cells, stromal cells — that interact to shape treatment response. This pipeline characterizes these populations and builds a classifier capable of distinguishing HPV_High from HPV_Low subtypes based on transcriptomic features and mutational burden, with potential utility in patient stratification.

---

## Author

**Shivani Patel**  
M.S. Bioinformatics Data Science, University of Delaware (GPA: 3.8)  
3+ years experience in computational biology, drug discovery, and multi-omics analysis

[![LinkedIn](https://img.shields.io/badge/LinkedIn-shivanip99-0077B5?style=flat&logo=linkedin)](https://linkedin.com/in/shivanip99)
[![GitHub](https://img.shields.io/badge/GitHub-shipatel302-181717?style=flat&logo=github)](https://github.com/shipatel302)
[![Email](https://img.shields.io/badge/Email-shivanip8369@gmail.com-D14836?style=flat&logo=gmail)](mailto:shivanip8369@gmail.com)


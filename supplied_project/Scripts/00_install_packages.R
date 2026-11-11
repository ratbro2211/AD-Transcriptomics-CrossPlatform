# =============================================================================
# 00_install_packages.R
# Project: Cross-Platform Transcriptomic Analysis of Alzheimer's Disease
# Author:  Suraj Rathi | SynBioAI Lab, Atria University
# Date:    2024
#
# Description:
#   Run this script ONCE before any other script in the project.
#   Installs all required CRAN and Bioconductor packages.
# =============================================================================

# ── Bioconductor core ──────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c(
  "DESeq2",          # RNA-seq differential expression
  "ashr",            # Adaptive shrinkage (apeglm alternative)
  "limma",           # Microarray differential expression
  "GEOquery",        # Fetch datasets from NCBI GEO
  "Biobase",         # ExpressionSet infrastructure
  "AnnotationDbi",   # Annotation interface
  "org.Hs.eg.db",    # Human gene annotation database
  "clusterProfiler", # GO and KEGG enrichment analysis
  "enrichplot",      # Enrichment visualisation
  "DOSE"             # Disease ontology enrichment
)

BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

# ── CRAN packages ──────────────────────────────────────────────────────────
cran_pkgs <- c(
  # Data wrangling
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "stringr",

  # Visualisation
  "ggplot2",
  "ggrepel",         # Non-overlapping labels on volcano plots
  "pheatmap",        # Heatmaps
  "RColorBrewer",
  "patchwork",       # Combine multiple ggplots
  "VennDiagram",     # Venn diagrams

  # Machine learning
  "caret",           # Unified ML interface + cross-validation
  "randomForest",    # Random Forest classifier
  "glmnet",          # LASSO logistic regression
  "e1071",           # SVM and Naive Bayes
  "class",           # KNN

  # Utilities
  "ggfortify",
  "scales"
)

install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

# ── Verify installation ────────────────────────────────────────────────────
all_pkgs <- c(bioc_pkgs, cran_pkgs)
loaded   <- sapply(all_pkgs, requireNamespace, quietly = TRUE)

cat("\n=== Installation Summary ===\n")
cat("Successfully installed:", sum(loaded), "/", length(all_pkgs), "packages\n")

if (any(!loaded)) {
  cat("FAILED to install:\n")
  cat(paste(" -", all_pkgs[!loaded]), sep = "\n")
} else {
  cat("All packages installed successfully. You are ready to run the analysis.\n")
}

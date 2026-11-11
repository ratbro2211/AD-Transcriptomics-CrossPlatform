# Run every entry point from the repository root.
config <- list(
  seed = 42L,
  padj_cutoff = 0.05,
  lfc_cutoff = 0.5,
  top_features = 100L,
  train_fraction = 0.8,
  rna_counts = "data/GSE104704/raw_counts.tsv",
  rna_provenance = "data/GSE104704/counts_provenance.tsv",
  # Baseline disease contrast. Add justified, non-confounded covariates only
  # after checking metadata, missingness, and design rank.
  rna_design = ~ condition,
  microarray_design = ~ condition
)

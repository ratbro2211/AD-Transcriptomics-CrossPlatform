# Run from repository root, not from inside Cross_platform/.
# Supplied historical files in results/ are NEVER overwritten.
source("R/common.R")
source("R/comparison_helpers.R")
run_comparison <- function() {
  rna_dir <- file.path("RNA_seq_analysis", "results")
  ma_dir <- file.path("Microarray_analysis", "results")
  threshold_files <- c(file.path(rna_dir, "thresholds.csv"), file.path(ma_dir, "thresholds.csv"))
  if (!all(file.exists(threshold_files))) stop("Run both corrected DE pipelines first.")
  thresholds <- lapply(threshold_files, read.csv)
  for (t in thresholds) {
    if (nrow(t) != 1 || t$padj_cutoff != config$padj_cutoff || t$lfc_cutoff != config$lfc_cutoff)
      stop("Thresholds changed: rerun both DE pipelines with the current config.R.")
  }
  rna <- read.csv(file.path(rna_dir, "de_all.csv"))
  ma <- read.csv(file.path(ma_dir, "de_all.csv"))
  comp <- compare_de(rna, ma, config$padj_cutoff, config$lfc_cutoff)
  out <- file.path("results", "rerun")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(comp$shared, file.path(out, "common_tested_genes.csv"), row.names = FALSE)
  write.csv(comp$overlap, file.path(out, "overlap_DE_genes.csv"), row.names = FALSE)
  write.csv(comp$summary, file.path(out, "cross_platform_summary.csv"), row.names = FALSE)
  # New ML scores are optional; do not fall back to historical leaky scores.
  ml_paths <- c(file.path(rna_dir, "ml_metrics.csv"), file.path(ma_dir, "ml_metrics.csv"))
  if (all(file.exists(ml_paths))) {
    ml <- lapply(ml_paths, read.csv)
    ml[[1]]$Platform <- "RNAseq"
    ml[[2]]$Platform <- "Microarray"
    write.csv(do.call(rbind, ml), file.path(out, "ml_comparison.csv"), row.names = FALSE)
  } else {
    message("New ML comparison skipped: both corrected ML runs are required.")
  }
  save_session(out)
  print(comp$summary)
}
run_comparison()

collapse_gene_results <- function(d) {
  required <- c("gene_symbol", "log2FoldChange", "padj")
  if (!all(required %in% names(d))) stop("Missing standardized DE columns.")
  d <- d[!is.na(d$gene_symbol) & nzchar(d$gene_symbol) &
           is.finite(d$log2FoldChange) & !is.na(d$padj), , drop = FALSE]
  # Pick one representative using abundance, NOT disease-associated p-values.
  abundance <- if ("baseMean" %in% names(d)) d$baseMean else
    if ("AveExpr" %in% names(d)) d$AveExpr else rep(0, nrow(d))
  abundance[!is.finite(abundance)] <- -Inf
  d <- d[order(d$gene_symbol, -abundance), , drop = FALSE]
  d[!duplicated(d$gene_symbol), , drop = FALSE]
}

compare_de <- function(rna, array, padj_cutoff, lfc_cutoff) {
  rna <- collapse_gene_results(rna)
  array <- collapse_gene_results(array)
  rna <- rna[, c("gene_symbol", "log2FoldChange", "padj")]
  array <- array[, c("gene_symbol", "log2FoldChange", "padj")]
  names(rna) <- c("Gene", "rnaseq_LFC", "rnaseq_padj")
  names(array) <- c("Gene", "array_logFC", "array_adjP")
  shared <- merge(rna, array, by = "Gene")
  hit_rna <- shared$rnaseq_padj < padj_cutoff & abs(shared$rnaseq_LFC) > lfc_cutoff
  hit_array <- shared$array_adjP < padj_cutoff & abs(shared$array_logFC) > lfc_cutoff
  overlap <- shared[hit_rna & hit_array, , drop = FALSE]
  overlap$concordant <- sign(overlap$rnaseq_LFC) == sign(overlap$array_logFC)
  safe_cor <- function(a, b) {
    if (length(a) < 3 || sd(a) == 0 || sd(b) == 0) return(NA_real_)
    cor(a, b, method = "pearson")
  }
  summary <- data.frame(
    metric = c("RNAseq_unique_DEGs", "Microarray_unique_DEGs", "Common_tested_genes",
               "Shared_DEGs", "Concordant_shared_DEGs", "Reversed_shared_DEGs",
               "Pearson_r_common_tested_genes", "Pearson_r_shared_DEGs",
               "padj_cutoff", "absolute_log2FC_cutoff"),
    value = c(
      sum(rna$rnaseq_padj < padj_cutoff & abs(rna$rnaseq_LFC) > lfc_cutoff),
      sum(array$array_adjP < padj_cutoff & abs(array$array_logFC) > lfc_cutoff),
      nrow(shared), nrow(overlap), sum(overlap$concordant), sum(!overlap$concordant),
      safe_cor(shared$rnaseq_LFC, shared$array_logFC),
      safe_cor(overlap$rnaseq_LFC, overlap$array_logFC), padj_cutoff, lfc_cutoff))
  list(shared = shared, overlap = overlap, summary = summary)
}

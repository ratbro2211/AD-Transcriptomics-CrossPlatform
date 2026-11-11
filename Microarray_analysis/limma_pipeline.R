# GSE132903: Illumina HumanHT-12 v4; 97 AD vs 98 non-demented controls.
# Analyze the GEO processed log-expression matrix. No neqc/voom/raw
# background-correction claim is made for this input.
source("R/common.R")
require_packages(c("GEOquery", "Biobase", "limma", "AnnotationDbi", "illuminaHumanv4.db"))
run_microarray <- function() {
  paths <- output_dirs("Microarray_analysis")
  dir.create("data/GSE132903", recursive = TRUE, showWarnings = FALSE)
  gse <- GEOquery::getGEO("GSE132903", GSEMatrix = TRUE, getGPL = FALSE,
                          destdir = "data/GSE132903")[[1]]
  if (Biobase::annotation(gse) != "GPL10558") stop("Unexpected platform: expected GPL10558.")
  meta <- Biobase::pData(gse)
  title <- as.character(meta$title)
  group <- ifelse(grepl("_AD_", title), "AD", ifelse(grepl("_ND_", title), "Control", NA_character_))
  meta$condition <- validate_groups(group, c(Control = 98, AD = 97))
  expr <- Biobase::exprs(gse)
  if (!is.numeric(expr) || any(!is.finite(expr))) stop("Non-finite expression matrix.")
  if (!identical(colnames(expr), rownames(meta))) stop("Expression and metadata not aligned.")
  # This processed series matrix is already on a log scale. Reject unexpected
  # values instead of silently log-transforming a different input.
  if (max(expr) > 50) stop("Unexpected non-log input: inspect submitter preprocessing.")
  valid <- intersect(rownames(expr), AnnotationDbi::keys(illuminaHumanv4.db::illuminaHumanv4.db,
                                                        keytype = "PROBEID"))
  if (!length(valid)) stop("No Illumina probe IDs could be annotated.")
  symbols <- AnnotationDbi::mapIds(illuminaHumanv4.db::illuminaHumanv4.db,
                                  keys = valid, keytype = "PROBEID", column = "SYMBOL",
                                  multiVals = function(x) if (length(unique(x)) == 1) x[1] else NA_character_)
  mapped <- unname(symbols[match(rownames(expr), names(symbols))])
  ok <- !is.na(mapped) & nzchar(mapped)
  write.csv(data.frame(probe_id = rownames(expr), gene_symbol = mapped),
            file.path(paths$results, "probe_annotation.csv"), row.names = FALSE)
  # Average probes by symbol without ranking them by disease-associated p-values.
  # This fixed, label-independent map is also used for ML input.
  gene_expr <- limma::avereps(expr[ok, , drop = FALSE], ID = mapped[ok])
  norm <- limma::normalizeBetweenArrays(gene_expr, method = "quantile")
  design <- model.matrix(config$microarray_design, meta)
  if (qr(design)$rank < ncol(design)) stop("Microarray design is not full rank.")
  coef <- which(colnames(design) == "conditionAD")
  if (length(coef) != 1) stop("Disease coefficient not available.")
  fit <- limma::eBayes(limma::lmFit(norm, design), trend = TRUE)
  df <- limma::topTable(fit, coef = coef, number = Inf, adjust.method = "BH")
  df$gene_id <- rownames(df)
  df$gene_symbol <- rownames(df)
  df$log2FoldChange <- df$logFC
  df$padj <- df$adj.P.Val
  export_de(df, paths)
  plot_de(df, paths, "GSE132903: AD vs non-demented Control")
  # Use the source log-expression matrix, NOT the full-cohort re-normalized
  # matrix. Additional ML scaling/filtering is fit within training folds.
  # The submitter's own preprocessing may still be cohort-wide; see README.
  saveRDS(list(expression = gene_expr, labels = meta$condition,
               sample_id = rownames(meta)), file.path(paths$cache, "ml_input.rds"))
  saveRDS(fit, file.path(paths$cache, "limma_fit.rds"))
  meta$sample_id <- rownames(meta)
  write.csv(meta, file.path(paths$results, "sample_metadata.csv"), row.names = FALSE)
  save_session(paths$results)
}
run_microarray()

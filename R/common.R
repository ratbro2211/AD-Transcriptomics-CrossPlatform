if (!file.exists("config.R")) {
  stop("Run from the repository root, e.g. Rscript RNA_seq_analysis/deseq2_pipeline.R")
}
source("config.R")
options(stringsAsFactors = FALSE)

require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install dependencies first: ", paste(missing, collapse = ", "))
}

output_dirs <- function(platform) {
  paths <- list(
    results = file.path(platform, "results"),
    figures = file.path(platform, "figures"),
    cache = file.path(platform, "RDS")
  )
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}

save_session <- function(path) {
  writeLines(capture.output(sessionInfo()), file.path(path, "session_info.txt"))
}

validate_counts <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x))
    stop("Counts must be a nonempty numeric matrix.")
  if (any(!is.finite(x)) || any(x < 0) || any(x != floor(x)))
    stop("DESeq2 needs nonnegative integer RAW counts: do not round normalized values.")
  if (any(x > .Machine$integer.max)) stop("Counts exceed the integer range.")
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
    stop("Counts require unique, nonempty gene and sample IDs.")
  if (any(colSums(x) <= 0)) stop("Empty count library.")
  invisible(TRUE)
}

validate_groups <- function(groups, expected) {
  observed <- table(factor(groups, levels = names(expected)))
  if (anyNA(groups) || length(groups) != sum(expected) ||
      !identical(as.integer(observed), as.integer(expected)))
    stop("Unexpected sample groups. Inspect metadata; never infer diagnosis by a broad 'AD' substring.")
  factor(groups, levels = c("Control", "AD"))
}

gene_symbols <- function(ids) {
  require_packages(c("AnnotationDbi", "org.Hs.eg.db"))
  clean <- sub("\\.[0-9]+$", "", ids)
  type <- if (all(grepl("^ENSG[0-9]+$", clean))) "ENSEMBL" else
    if (all(grepl("^[0-9]+$", clean))) "ENTREZID" else "SYMBOL"
  valid <- intersect(clean, AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = type))
  result <- rep(NA_character_, length(ids))
  if (length(valid)) {
    mapped <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = valid,
                                   column = "SYMBOL", keytype = type,
                                   multiVals = function(x) if (length(unique(x)) == 1) x[1] else NA_character_)
    result <- unname(mapped[match(clean, names(mapped))])
  }
  result
}

export_de <- function(df, paths) {
  stopifnot(all(c("gene_id", "gene_symbol", "log2FoldChange", "padj") %in% names(df)))
  df <- df[order(df$padj, na.last = TRUE), , drop = FALSE]
  sig <- df[!is.na(df$padj) & df$padj < config$padj_cutoff &
              abs(df$log2FoldChange) > config$lfc_cutoff, , drop = FALSE]
  write.csv(df, file.path(paths$results, "de_all.csv"), row.names = FALSE)
  write.csv(sig, file.path(paths$results, "de_significant.csv"), row.names = FALSE)
  write.csv(data.frame(padj_cutoff = config$padj_cutoff, lfc_cutoff = config$lfc_cutoff),
            file.path(paths$results, "thresholds.csv"), row.names = FALSE)
  message("Significant rows: ", nrow(sig), "; inspect unique mapped gene counts separately.")
  invisible(sig)
}

plot_de <- function(df, paths, title) {
  png(file.path(paths$figures, "volcano.png"), width = 1200, height = 850, res = 130)
  on.exit(dev.off())
  ok <- is.finite(df$log2FoldChange) & !is.na(df$padj)
  d <- df[ok, , drop = FALSE]
  hit <- d$padj < config$padj_cutoff & abs(d$log2FoldChange) > config$lfc_cutoff
  plot(d$log2FoldChange, -log10(pmax(d$padj, .Machine$double.xmin)),
       pch = 16, cex = 0.45, col = ifelse(hit, "#4477AA", "#888888"),
       xlab = "log2 fold change (AD / Control)", ylab = "-log10 BH-adjusted p",
       main = title)
  abline(v = c(-1, 1) * config$lfc_cutoff, h = -log10(config$padj_cutoff),
         lty = 2, col = "grey50")
}

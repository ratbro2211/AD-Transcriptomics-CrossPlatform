# GSE104704: 12 AD vs 10 older controls. Eight Young samples are excluded.
# The submitter's supplementary table contains library-size-adjusted values,
# NOT raw counts. Integer coercion/rounding is not a valid recovery of counts.
source("R/common.R")
require_packages(c("GEOquery", "Biobase", "DESeq2", "ashr", "AnnotationDbi", "org.Hs.eg.db"))
run_rnaseq <- function() {
  paths <- output_dirs("RNA_seq_analysis")
  dir.create("data/GSE104704", recursive = TRUE, showWarnings = FALSE)
  gse <- GEOquery::getGEO("GSE104704", GSEMatrix = TRUE,
                          getGPL = FALSE, destdir = "data/GSE104704")[[1]]
  meta <- Biobase::pData(gse)
  title <- as.character(meta$title)
  group <- ifelse(grepl("\\.RNA\\.AD$", title), "AD",
                  ifelse(grepl("\\.RNA\\.Old$", title), "Control",
                         ifelse(grepl("\\.RNA\\.Young$", title), "Young", NA_character_)))
  if (anyNA(group)) stop("Unrecognized GEO title. Inspect metadata before proceeding.")
  selected <- group != "Young"
  meta <- meta[selected, , drop = FALSE]
  meta$condition <- validate_groups(group[selected], c(Control = 10, AD = 12))
  meta$sample_id <- rownames(meta)
  write.csv(meta, file.path(paths$results, "sample_metadata.csv"), row.names = FALSE)

  if (!file.exists(config$rna_counts) || !file.exists(config$rna_provenance)) {
    stop("Metadata downloaded, but verified raw counts and provenance are required. ",
         "See data/README.md. Do NOT use the normalized supplementary table or round its values.")
  }
  provenance <- read.delim(config$rna_provenance, check.names = FALSE)
  fields <- c("dataset", "input_type", "source_url", "processing")
  if (nrow(provenance) != 1 || !all(fields %in% names(provenance)) ||
      anyNA(provenance[fields]) || any(!nzchar(unlist(provenance[fields]))) ||
      provenance$dataset != "GSE104704" || provenance$input_type != "raw_gene_counts")
    stop("Complete the raw-count provenance record; see data/README.md.")
  # A provenance declaration documents origin; integer validation alone does
  # not prove a file is raw counts. Review its upstream processing independently.
  tab <- read.delim(config$rna_counts, check.names = FALSE)
  if (names(tab)[1] != "gene_id") stop("First counts column must be gene_id.")
  ids <- tab$gene_id
  counts <- as.matrix(tab[, -1, drop = FALSE])
  rownames(counts) <- ids
  validate_counts(counts)
  if (!all(meta$sample_id %in% colnames(counts)))
    stop("Counts columns must be GEO GSM IDs; supply an explicit mapping upstream if necessary.")
  counts <- counts[, meta$sample_id, drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(meta)))
  design <- model.matrix(config$rna_design, meta)
  if (qr(design)$rank < ncol(design)) stop("RNA-seq design is not full rank.")
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, meta, design = config$rna_design)
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  if (nrow(dds) < 10) stop("Too few expressed genes after prefiltering.")
  dds <- DESeq2::DESeq(dds)
  coef <- "condition_AD_vs_Control"
  if (!coef %in% DESeq2::resultsNames(dds)) stop("Disease coefficient not available.")
  res <- DESeq2::results(dds, name = coef, alpha = config$padj_cutoff)
  shrunk <- DESeq2::lfcShrink(dds, coef = coef, res = res, type = "ashr")
  df <- as.data.frame(shrunk)
  df$gene_id <- rownames(df)
  df$gene_symbol <- gene_symbols(df$gene_id)
  # ML does not use full-cohort DEGs, fitted dispersions, or a full-cohort VST.
  # Export ALL raw genes for training-fold feature filtering later.
  saveRDS(list(raw_counts = counts, labels = meta$condition,
               sample_id = meta$sample_id), file.path(paths$cache, "ml_input.rds"))
  saveRDS(dds, file.path(paths$cache, "dds.rds"))
  export_de(df, paths)
  plot_de(df, paths, "GSE104704: AD vs older Control")
  file.copy(config$rna_provenance, file.path(paths$results, "counts_provenance.tsv"), overwrite = TRUE)
  save_session(paths$results)
}
run_rnaseq()

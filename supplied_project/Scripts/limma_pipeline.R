# =============================================================================
# limma_pipeline.R
# Project: Cross-Platform Transcriptomic Analysis of Alzheimer's Disease
# Dataset: GSE132903 (Microarray, Affymetrix, n = 195)
# Author:  Suraj Rathi | SynBioAI Lab, Atria University
#
# Description:
#   Performs differential expression analysis on the GSE132903 microarray dataset
#   comparing Alzheimer's Disease (AD) vs. healthy Control samples.
#   Uses the limma package with eBayes moderated t-statistics.
#   Includes background correction, quantile normalisation, and log2 transformation.
#
# Output:
#   - figures/microarray_volcano.png       Volcano plot of DEGs
#   - results/microarray_limma_all.csv     Full limma results table
#   - results/microarray_limma_sig.csv     Significant DEGs
#   - RDS/eset_norm.rds                   Normalised ExpressionSet
#   - RDS/limma_fit.rds                   limma fit object for downstream use
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(Biobase)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
})

# ── 1. Configuration ──────────────────────────────────────────────────────
GEO_ID       <- "GSE132903"
PADJ_CUTOFF  <- 0.05
LFC_CUTOFF   <- 1
N_LABEL      <- 15
OUT_DIR      <- "figures"
RES_DIR      <- "results"
RDS_DIR      <- "RDS"

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(RES_DIR, showWarnings = FALSE)
dir.create(RDS_DIR, showWarnings = FALSE)

set.seed(42)

# ── 2. Download Data from GEO ─────────────────────────────────────────────
cat("Downloading", GEO_ID, "from GEO...\n")
gse <- getGEO(GEO_ID, GSEMatrix = TRUE, AnnotGPL = TRUE)
gse <- gse[[1]]

cat("Dataset dimensions:", nrow(gse), "probes ×", ncol(gse), "samples\n")

# ── 3. Build Condition Vector ─────────────────────────────────────────────
meta <- pData(gse)
cat("Sample metadata columns:\n"); print(colnames(meta))

condition_col <- grep("characteristics|disease|status|title",
                      colnames(meta), ignore.case = TRUE, value = TRUE)[1]
cat("Using column:", condition_col, "\n")

condition_raw <- as.character(meta[[condition_col]])
cat("Unique values:\n"); print(unique(condition_raw))

# ADJUST PATTERNS to match your GEO metadata
condition <- case_when(
  grepl("alzheimer|AD|disease", condition_raw, ignore.case = TRUE) ~ "AD",
  grepl("control|normal|healthy|CTL", condition_raw, ignore.case = TRUE) ~ "Control",
  TRUE ~ NA_character_
)

if (any(is.na(condition))) {
  stop("Could not assign all samples to AD/Control. ",
       "Inspect 'condition_raw' and adjust grepl patterns.")
}

condition <- factor(condition, levels = c("Control", "AD"))
cat("\nCondition assignments:\n"); print(table(condition))
gse$condition <- condition

# ── 4. Normalisation ──────────────────────────────────────────────────────
# Log2 transform if not already done
expr <- exprs(gse)

if (max(expr, na.rm = TRUE) > 50) {
  # Data appears to be on linear scale — log2 transform
  expr <- log2(expr + 1)
  cat("Log2 transformation applied\n")
} else {
  cat("Data already on log scale\n")
}

# Quantile normalisation via limma
expr_norm <- normalizeBetweenArrays(expr, method = "quantile")
exprs(gse) <- expr_norm

saveRDS(gse, file.path(RDS_DIR, "eset_norm.rds"))
cat("Normalised ExpressionSet saved\n")

# ── 5. Probe Annotation ───────────────────────────────────────────────────
# Map probe IDs to gene symbols using the platform annotation in gse
feat <- fData(gse)
cat("Feature data columns:\n"); print(colnames(feat))

# Common column names for gene symbols in Affymetrix annotations
symbol_col <- grep("Gene.Symbol|GENE_SYMBOL|Symbol|Gene symbol",
                   colnames(feat), ignore.case = TRUE, value = TRUE)[1]
cat("Using gene symbol column:", symbol_col, "\n")

feat$gene_symbol <- as.character(feat[[symbol_col]])

# Remove probes without a gene symbol or mapping to multiple genes
feat$gene_symbol <- gsub(" ///.*", "", feat$gene_symbol)   # keep first symbol
feat$gene_symbol[feat$gene_symbol == "" | feat$gene_symbol == "---"] <- NA

# Assign annotation back
fData(gse) <- feat

# ── 6. Fit limma Model ────────────────────────────────────────────────────
cat("\nFitting limma model...\n")

design <- model.matrix(~ 0 + condition, data = pData(gse))
colnames(design) <- levels(condition)

contrast_matrix <- makeContrasts(
  AD_vs_Control = AD - Control,
  levels        = design
)

# Fit linear model
fit  <- lmFit(exprs(gse), design)
fit2 <- contrasts.fit(fit, contrast_matrix)

# eBayes moderation — borrowing information across probes
fit2 <- eBayes(fit2, trend = TRUE)   # trend = TRUE recommended for microarrays

saveRDS(fit2, file.path(RDS_DIR, "limma_fit.rds"))

# ── 7. Extract Results ────────────────────────────────────────────────────
top_all <- topTable(fit2,
                    coef       = "AD_vs_Control",
                    number     = Inf,
                    adjust.method = "BH",
                    sort.by    = "P")

# Merge gene symbols
top_all$probe_id    <- rownames(top_all)
top_all$gene_symbol <- fData(gse)[rownames(top_all), "gene_symbol"]

# Use probe_id as fallback when symbol is NA
top_all$label <- ifelse(is.na(top_all$gene_symbol),
                        top_all$probe_id,
                        top_all$gene_symbol)

write.csv(top_all,
          file = file.path(RES_DIR, "microarray_limma_all.csv"),
          row.names = FALSE)

# Significant DEGs
sig_df <- top_all %>%
  filter(!is.na(adj.P.Val),
         adj.P.Val < PADJ_CUTOFF,
         abs(logFC) > LFC_CUTOFF)

write.csv(sig_df,
          file = file.path(RES_DIR, "microarray_limma_sig.csv"),
          row.names = FALSE)

cat("\n=== limma Summary ===\n")
cat("Total probes tested:", nrow(top_all), "\n")
cat("Significant DEGs   :", nrow(sig_df), "\n")
cat("  Upregulated      :", sum(sig_df$logFC > 0), "\n")
cat("  Downregulated    :", sum(sig_df$logFC < 0), "\n")
cat("Top upregulated    :", sig_df$label[which.max(sig_df$logFC)], "\n")

# ── 8. Volcano Plot ───────────────────────────────────────────────────────
plot_df <- top_all %>%
  filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
  mutate(
    neg_log10_padj = -log10(adj.P.Val),
    significance   = case_when(
      adj.P.Val < PADJ_CUTOFF & logFC >  LFC_CUTOFF ~ "Upregulated",
      adj.P.Val < PADJ_CUTOFF & logFC < -LFC_CUTOFF ~ "Downregulated",
      TRUE                                           ~ "Not Significant"
    )
  )

label_genes <- plot_df %>%
  filter(significance != "Not Significant") %>%
  arrange(adj.P.Val) %>%
  slice_head(n = N_LABEL)

volcano_plot <- ggplot(plot_df, aes(x = logFC,
                                    y = neg_log10_padj,
                                    colour = significance)) +
  geom_point(alpha = 0.4, size = 1.0) +
  scale_colour_manual(
    values = c(
      "Upregulated"     = "#E41A1C",
      "Downregulated"   = "#377EB8",
      "Not Significant" = "grey70"
    )
  ) +
  geom_text_repel(
    data         = label_genes,
    aes(label    = label),
    size         = 3,
    max.overlaps = 20,
    colour       = "black",
    fontface     = "italic"
  ) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
             linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(PADJ_CUTOFF),
             linetype = "dashed", colour = "grey40") +
  labs(
    title    = "Differential Probe Expression — GSE132903 (Microarray)",
    subtitle = paste0("AD vs. Control | limma + eBayes | ",
                      nrow(sig_df), " DEGs (adj.P < 0.05, |logFC| > 1)"),
    x        = "log2 Fold Change (AD / Control)",
    y        = expression(-log[10](adj.P.Val)),
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(OUT_DIR, "microarray_volcano.png"),
       plot  = volcano_plot,
       width = 8, height = 6, dpi = 300)

cat("Volcano plot saved\n")
cat("\nDone. Next: run pathway_enrichment.R\n")

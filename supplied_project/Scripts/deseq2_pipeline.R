# =============================================================================
# deseq2_pipeline.R
# Project: Cross-Platform Transcriptomic Analysis of Alzheimer's Disease
# Dataset: GSE104704 (RNA-seq, Illumina, n = 22)
# Author:  Suraj Rathi | SynBioAI Lab, Atria University
#
# Description:
#   Performs differential expression analysis on the GSE104704 RNA-seq dataset
#   comparing Alzheimer's Disease (AD) vs. healthy Control brain samples.
#   Uses DESeq2 with ashr adaptive shrinkage for log2 fold-change estimation.
#
# Output:
#   - figures/rnaseq_volcano.png      Volcano plot of DEGs
#   - results/rnaseq_deseq2_all.csv   Full DESeq2 results table
#   - results/rnaseq_deseq2_sig.csv   Significant DEGs (padj < 0.05, |lfc| > 1)
#   - RDS/dds.rds                     DESeqDataSet object for downstream use
#   - RDS/res_shrunk.rds              Shrinkage results for ML pipeline
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ashr)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
})

# ── 1. Configuration ──────────────────────────────────────────────────────
GEO_ID       <- "GSE104704"
PADJ_CUTOFF  <- 0.05
LFC_CUTOFF   <- 1        # |log2FoldChange| threshold
N_LABEL      <- 15       # number of top genes to label on volcano plot
OUT_DIR      <- "figures"
RES_DIR      <- "results"
RDS_DIR      <- "RDS"

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(RES_DIR, showWarnings = FALSE)
dir.create(RDS_DIR, showWarnings = FALSE)

set.seed(42)

# ── 2. Download Data from GEO ─────────────────────────────────────────────
cat("Downloading", GEO_ID, "from GEO...\n")
gse <- getGEO(GEO_ID, GSEMatrix = TRUE, AnnotGPL = FALSE)
gse <- gse[[1]]

# Extract count matrix (rows = genes, cols = samples)
count_matrix <- exprs(gse)
mode(count_matrix) <- "integer"   # DESeq2 requires integer counts

# Extract sample metadata
meta <- pData(gse)

# ── 3. Build Condition Vector ─────────────────────────────────────────────
# Inspect column names to identify the condition field
# Typical GEO fields: "characteristics_ch1", "title", "source_name_ch1"
cat("Sample metadata columns:\n")
print(colnames(meta))

# Adjust this selector based on actual GEO metadata structure
# Common pattern: grepl "Alzheimer" or "control" in characteristics_ch1
condition_col <- grep("characteristics|disease|status|title",
                      colnames(meta), ignore.case = TRUE, value = TRUE)[1]
cat("Using condition column:", condition_col, "\n")

condition_raw <- as.character(meta[[condition_col]])
cat("Unique values:\n"); print(unique(condition_raw))

# Map to binary factor — ADJUST PATTERNS to match your GEO metadata
condition <- case_when(
  grepl("alzheimer|AD|disease", condition_raw, ignore.case = TRUE) ~ "AD",
  grepl("control|normal|healthy|CTL", condition_raw, ignore.case = TRUE) ~ "Control",
  TRUE ~ NA_character_
)

if (any(is.na(condition))) {
  stop("Could not assign condition for all samples. ",
       "Check 'condition_raw' and adjust the grepl patterns above.")
}

condition <- factor(condition, levels = c("Control", "AD"))
cat("\nCondition assignments:\n"); print(table(condition))

# ── 4. Build DESeqDataSet ─────────────────────────────────────────────────
col_data <- data.frame(
  condition = condition,
  row.names = colnames(count_matrix)
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = col_data,
  design    = ~ condition
)

# Pre-filter: keep genes with >= 10 total counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("\nGenes after pre-filtering:", nrow(dds), "\n")

# ── 5. Run DESeq2 ─────────────────────────────────────────────────────────
cat("Running DESeq2...\n")
dds <- DESeq(dds)
saveRDS(dds, file.path(RDS_DIR, "dds.rds"))

# ── 6. Apply ashr Shrinkage ───────────────────────────────────────────────
# ashr is recommended over apeglm for datasets with small n
res_shrunk <- lfcShrink(dds,
                        coef    = "condition_AD_vs_Control",
                        type    = "ashr")
saveRDS(res_shrunk, file.path(RDS_DIR, "res_shrunk.rds"))

# ── 7. Export Full Results ────────────────────────────────────────────────
res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(padj)

write.csv(res_df,
          file = file.path(RES_DIR, "rnaseq_deseq2_all.csv"),
          row.names = FALSE)

# Significant DEGs
sig_df <- res_df %>%
  filter(!is.na(padj),
         padj  < PADJ_CUTOFF,
         abs(log2FoldChange) > LFC_CUTOFF)

write.csv(sig_df,
          file = file.path(RES_DIR, "rnaseq_deseq2_sig.csv"),
          row.names = FALSE)

cat("\n=== DESeq2 Summary ===\n")
cat("Total genes tested:", nrow(res_df), "\n")
cat("Significant DEGs  :", nrow(sig_df), "\n")
cat("  Upregulated     :", sum(sig_df$log2FoldChange > 0), "\n")
cat("  Downregulated   :", sum(sig_df$log2FoldChange < 0), "\n")
cat("Top downregulated gene:", sig_df$gene_id[which.min(sig_df$log2FoldChange)], "\n")

# ── 8. Volcano Plot ───────────────────────────────────────────────────────
plot_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance   = case_when(
      padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Upregulated",
      padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Downregulated",
      TRUE                                               ~ "Not Significant"
    )
  )

# Pick top genes to label (most significant up and down)
label_genes <- plot_df %>%
  filter(significance != "Not Significant") %>%
  arrange(padj) %>%
  slice_head(n = N_LABEL)

volcano_plot <- ggplot(plot_df, aes(x = log2FoldChange,
                                    y = neg_log10_padj,
                                    colour = significance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_colour_manual(
    values = c(
      "Upregulated"     = "#E41A1C",
      "Downregulated"   = "#377EB8",
      "Not Significant" = "grey70"
    )
  ) +
  geom_text_repel(
    data        = label_genes,
    aes(label   = gene_id),
    size        = 3,
    max.overlaps= 20,
    colour      = "black",
    fontface    = "italic"
  ) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
             linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(PADJ_CUTOFF),
             linetype = "dashed", colour = "grey40") +
  labs(
    title    = "Differential Gene Expression — GSE104704 (RNA-seq)",
    subtitle = paste0("AD vs. Control | DESeq2 + ashr | ",
                      nrow(sig_df), " DEGs (padj < 0.05, |log2FC| > 1)"),
    x        = "Shrinkage log2 Fold Change (AD / Control)",
    y        = expression(-log[10](padj)),
    colour   = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(OUT_DIR, "rnaseq_volcano.png"),
       plot   = volcano_plot,
       width  = 8, height = 6, dpi = 300)

cat("Volcano plot saved to", file.path(OUT_DIR, "rnaseq_volcano.png"), "\n")
cat("\nDone. Next: run pathway_enrichment.R\n")

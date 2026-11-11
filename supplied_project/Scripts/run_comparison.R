# =============================================================================
# run_comparison.R
# Project: Cross-Platform Transcriptomic Analysis of Alzheimer's Disease
# Author:  Suraj Rathi | SynBioAI Lab, Atria University
#
# Description:
#   Cross-platform comparison of RNA-seq (GSE104704) and microarray (GSE132903)
#   differential expression results. Includes:
#     1. Venn diagram of DEG overlap
#     2. Log2FC concordance scatter plot (with Pearson r)
#     3. Side-by-side ML performance comparison
#     4. Summary statistics table
#
# Prerequisites:
#   Both RNA_seq_analysis/ and Microarray_analysis/ pipelines must be complete.
#   Required files:
#     - RNA_seq_analysis/results/rnaseq_deseq2_sig.csv
#     - Microarray_analysis/results/microarray_limma_sig.csv
#     - RNA_seq_analysis/results/rnaseq_deseq2_all.csv
#     - Microarray_analysis/results/microarray_limma_all.csv
#     - RNA_seq_analysis/results/rnaseq_ml_metrics.csv
#     - Microarray_analysis/results/microarray_ml_metrics.csv
#
# Output:
#   - figures/venn_diagram.png
#   - figures/lfc_concordance_scatter.png
#   - figures/ml_cross_platform_comparison.png
#   - results/shared_degs.csv
#   - results/cross_platform_summary.csv
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(VennDiagram)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(grid)
})

# ── 1. Configuration & Paths ──────────────────────────────────────────────
RNA_RES  <- "../RNA_seq_analysis/results"
MA_RES   <- "../Microarray_analysis/results"
OUT_DIR  <- "figures"
RES_DIR  <- "results"

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(RES_DIR, showWarnings = FALSE)

# ── 2. Load Results ───────────────────────────────────────────────────────
rna_sig  <- read.csv(file.path(RNA_RES, "rnaseq_deseq2_sig.csv"),
                     stringsAsFactors = FALSE)
ma_sig   <- read.csv(file.path(MA_RES,  "microarray_limma_sig.csv"),
                     stringsAsFactors = FALSE)

rna_all  <- read.csv(file.path(RNA_RES, "rnaseq_deseq2_all.csv"),
                     stringsAsFactors = FALSE)
ma_all   <- read.csv(file.path(MA_RES,  "microarray_limma_all.csv"),
                     stringsAsFactors = FALSE)

rna_ml   <- read.csv(file.path(RNA_RES, "rnaseq_ml_metrics.csv"),
                     stringsAsFactors = FALSE)
ma_ml    <- read.csv(file.path(MA_RES,  "microarray_ml_metrics.csv"),
                     stringsAsFactors = FALSE)

cat("RNA-seq DEGs   :", nrow(rna_sig), "\n")
cat("Microarray DEGs:", nrow(ma_sig),  "\n")

# ── 3. Gene Sets for Venn Diagram ─────────────────────────────────────────
rna_genes <- unique(rna_sig$gene_id)

# Microarray: use gene_symbol (remove NAs and blanks)
ma_genes <- unique(na.omit(ma_sig$gene_symbol))
ma_genes <- ma_genes[ma_genes != "" & ma_genes != "---"]

# RNA-seq gene IDs are already symbols; find overlap
shared_genes <- intersect(rna_genes, ma_genes)

cat("\nRNA-seq unique genes  :", length(rna_genes),   "\n")
cat("Microarray unique genes:", length(ma_genes),     "\n")
cat("Shared genes           :", length(shared_genes), "\n")

# Save shared genes table with LFC from both platforms
rna_lfc <- rna_all %>%
  filter(gene_id %in% shared_genes) %>%
  select(gene_id, rna_lfc = log2FoldChange, rna_padj = padj)

ma_lfc <- ma_all %>%
  filter(gene_symbol %in% shared_genes) %>%
  group_by(gene_symbol) %>%
  slice_min(adj.P.Val, n = 1) %>%          # one row per gene (best probe)
  ungroup() %>%
  select(gene_id = gene_symbol, ma_lfc = logFC, ma_padj = adj.P.Val)

shared_df <- inner_join(rna_lfc, ma_lfc, by = "gene_id") %>%
  mutate(direction_concordant = sign(rna_lfc) == sign(ma_lfc))

write.csv(shared_df,
          file = file.path(RES_DIR, "shared_degs.csv"),
          row.names = FALSE)

cat("Concordant direction   :", sum(shared_df$direction_concordant),
    "/", nrow(shared_df), "\n")
cat("Pearson r (LFC)        :",
    round(cor(shared_df$rna_lfc, shared_df$ma_lfc, method = "pearson"), 3), "\n")

# ── 4. Venn Diagram ───────────────────────────────────────────────────────
venn_path <- file.path(OUT_DIR, "venn_diagram.png")

venn.plot <- venn.diagram(
  x          = list(`RNA-seq\n(GSE104704)` = rna_genes,
                    `Microarray\n(GSE132903)` = ma_genes),
  filename   = NULL,
  fill       = c("#4E79A7", "#F28E2B"),
  alpha      = 0.5,
  cex        = 1.6,
  cat.cex    = 1.3,
  cat.fontface = "bold",
  fontface   = "bold",
  main       = "Shared DEGs: RNA-seq vs. Microarray",
  main.cex   = 1.4,
  main.fontface = "bold"
)

png(venn_path, width = 1800, height = 1600, res = 300)
grid.draw(venn.plot)
dev.off()
cat("Venn diagram saved\n")

# ── 5. LFC Concordance Scatter ────────────────────────────────────────────
r_val <- cor(shared_df$rna_lfc, shared_df$ma_lfc, method = "pearson")

# Label notable genes
label_genes <- shared_df %>%
  arrange(desc(abs(rna_lfc) + abs(ma_lfc))) %>%
  slice_head(n = 12)

lfc_scatter <- ggplot(shared_df, aes(x = rna_lfc, y = ma_lfc)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = direction_concordant), size = 2.5, alpha = 0.8) +
  scale_colour_manual(
    values = c("TRUE" = "#2CA02C", "FALSE" = "#D62728"),
    labels = c("TRUE" = "Concordant", "FALSE" = "Reversed"),
    name   = "Direction"
  ) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30",
              linetype = "solid", linewidth = 0.8) +
  geom_text_repel(
    data        = label_genes,
    aes(label   = gene_id),
    size        = 3,
    max.overlaps= 15,
    fontface    = "italic",
    colour      = "black"
  ) +
  annotate("text",
           x = min(shared_df$rna_lfc) * 0.85,
           y = max(shared_df$ma_lfc)  * 0.90,
           label = paste0("Pearson r = ", round(r_val, 3)),
           hjust = 0, size = 4, fontface = "italic") +
  labs(
    title    = "Cross-Platform LFC Concordance — 47 Shared DEGs",
    subtitle = "GSE104704 (RNA-seq) vs. GSE132903 (Microarray) | AD vs. Control",
    x        = "RNA-seq log2FC (DESeq2 + ashr)",
    y        = "Microarray log2FC (limma + eBayes)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(OUT_DIR, "lfc_concordance_scatter.png"),
       plot = lfc_scatter, width = 8, height = 7, dpi = 300)
cat("LFC concordance scatter saved\n")

# ── 6. ML Cross-Platform Comparison ──────────────────────────────────────
# Combine both metric tables
rna_plot <- rna_ml %>%
  mutate(Platform = "RNA-seq\n(LOOCV, n=22)") %>%
  select(Model, Accuracy, AUC, Platform)

ma_plot  <- ma_ml %>%
  mutate(Platform = "Microarray\n(80/20 split, n=195)") %>%
  select(Model, Accuracy, AUC, Platform)

combined <- bind_rows(rna_plot, ma_plot) %>%
  pivot_longer(cols = c(Accuracy, AUC),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Model    = factor(Model, levels = rna_ml$Model[order(rna_ml$AUC)]),
         Platform = factor(Platform))

ml_compare <- ggplot(combined,
                     aes(x = Model, y = Value,
                         fill = Metric, alpha = Platform)) +
  geom_col(position = position_dodge(0.75), width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", Value)),
            position = position_dodge(0.75),
            vjust = -0.4, size = 2.8) +
  scale_fill_manual(values  = c(Accuracy = "#4E79A7", AUC = "#F28E2B")) +
  scale_alpha_manual(values = c(0.65, 1.0)) +
  facet_wrap(~ Platform, ncol = 2) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2)) +
  coord_flip() +
  labs(
    title    = "ML Performance: RNA-seq vs. Microarray",
    subtitle = "Comparison of 5 classifiers across both platforms",
    x        = NULL, y = "Score", fill = NULL
  ) +
  guides(alpha = "none") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(OUT_DIR, "ml_cross_platform_comparison.png"),
       plot = ml_compare, width = 12, height = 5, dpi = 300)
cat("ML cross-platform comparison saved\n")

# ── 7. Summary Statistics Table ───────────────────────────────────────────
summary_df <- data.frame(
  Metric   = c(
    "Total DEGs (RNA-seq)",
    "Total DEGs (Microarray)",
    "Shared DEGs",
    "Concordant direction",
    "Direction reversal",
    "Pearson r (LFC)",
    "Best RNA-seq model (AUC)",
    "Best Microarray model (AUC)"
  ),
  Value    = c(
    nrow(rna_sig),
    nrow(ma_sig),
    nrow(shared_df),
    sum(shared_df$direction_concordant),
    sum(!shared_df$direction_concordant),
    round(r_val, 3),
    paste0(rna_ml$Model[which.max(rna_ml$AUC)],
           " (AUC = ", max(rna_ml$AUC), ")"),
    paste0(ma_ml$Model[which.max(ma_ml$AUC)],
           " (AUC = ", max(ma_ml$AUC), ")")
  )
)

write.csv(summary_df,
          file = file.path(RES_DIR, "cross_platform_summary.csv"),
          row.names = FALSE)

cat("\n=== Cross-Platform Summary ===\n")
print(summary_df, row.names = FALSE)
cat("\nAll cross-platform figures and tables saved.\n")
cat("Analysis complete!\n")

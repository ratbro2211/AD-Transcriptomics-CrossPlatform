# =============================================================================
# pathway_enrichment.R
# Project: Cross-Platform Transcriptomic Analysis of Alzheimer's Disease
# Dataset: GSE132903 (Microarray)
# Author:  Suraj Rathi | SynBioAI Lab, Atria University
#
# Description:
#   GO Biological Process and KEGG pathway enrichment for microarray DEGs.
#   Produces bar charts (instead of dot plots) for visual variety.
#
# Prerequisites:
#   Run limma_pipeline.R first (requires results/microarray_limma_sig.csv)
#
# Output:
#   - figures/microarray_go_bp_barchart.png
#   - figures/microarray_kegg_barchart.png
#   - results/microarray_go_results.csv
#   - results/microarray_kegg_results.csv
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(forcats)
})

# ── 1. Configuration ──────────────────────────────────────────────────────
PADJ_CUTOFF <- 0.05
MIN_GS      <- 10
MAX_GS      <- 500
TOP_N       <- 20
OUT_DIR     <- "figures"
RES_DIR     <- "results"

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(RES_DIR, showWarnings = FALSE)

# ── 2. Load DEGs ──────────────────────────────────────────────────────────
sig_df <- read.csv(file.path(RES_DIR, "microarray_limma_sig.csv"),
                   stringsAsFactors = FALSE)

cat("Loaded", nrow(sig_df), "significant DEGs (microarray)\n")

# Use gene_symbol column; fall back to probe_id if symbol is missing
gene_symbols <- unique(na.omit(sig_df$gene_symbol))
gene_symbols <- gene_symbols[gene_symbols != "" & gene_symbols != "---"]
cat("Unique gene symbols:", length(gene_symbols), "\n")

# ── 3. Gene ID Conversion ─────────────────────────────────────────────────
id_map <- bitr(
  geneID   = gene_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
entrez_ids <- unique(id_map$ENTREZID)
cat("Mapped to", length(entrez_ids), "Entrez IDs\n")

# ── 4. GO Biological Process ──────────────────────────────────────────────
cat("\nRunning GO-BP enrichment (microarray)...\n")

ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = PADJ_CUTOFF,
  qvalueCutoff  = 0.2,
  minGSSize     = MIN_GS,
  maxGSSize     = MAX_GS,
  readable      = TRUE
)
ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

go_results <- as.data.frame(ego)
write.csv(go_results,
          file = file.path(RES_DIR, "microarray_go_results.csv"),
          row.names = FALSE)
cat("GO-BP:", nrow(go_results), "significant terms\n")

# Bar chart
go_top <- go_results %>%
  arrange(p.adjust) %>%
  slice_head(n = TOP_N) %>%
  mutate(Description = fct_reorder(Description, -p.adjust))

go_bar <- ggplot(go_top, aes(x = Description,
                              y = Count,
                              fill = p.adjust)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradient(low  = "#D73027",
                      high = "#FDAE61",
                      name = "Adjusted\np-value") +
  labs(
    title    = "GO Biological Process Enrichment — GSE132903 (Microarray)",
    subtitle = paste0("Top ", TOP_N, " terms | adj. p-value < 0.05"),
    x        = NULL,
    y        = "Gene Count"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "microarray_go_bp_barchart.png"),
       plot = go_bar, width = 10, height = 7, dpi = 300)
cat("GO bar chart saved\n")

# ── 5. KEGG Pathway Enrichment ────────────────────────────────────────────
cat("\nRunning KEGG enrichment (microarray)...\n")

ekegg <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = PADJ_CUTOFF,
  qvalueCutoff  = 0.2,
  minGSSize     = MIN_GS,
  maxGSSize     = MAX_GS
)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

kegg_results <- as.data.frame(ekegg)
write.csv(kegg_results,
          file = file.path(RES_DIR, "microarray_kegg_results.csv"),
          row.names = FALSE)
cat("KEGG:", nrow(kegg_results), "significant pathways\n")

# Bar chart
n_show <- min(TOP_N, nrow(kegg_results))
kegg_top <- kegg_results %>%
  arrange(p.adjust) %>%
  slice_head(n = n_show) %>%
  mutate(Description = fct_reorder(Description, -p.adjust))

kegg_bar <- ggplot(kegg_top, aes(x = Description,
                                  y = Count,
                                  fill = p.adjust)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_gradient(low  = "#2166AC",
                      high = "#92C5DE",
                      name = "Adjusted\np-value") +
  labs(
    title    = "KEGG Pathway Enrichment — GSE132903 (Microarray)",
    subtitle = paste0("Top ", n_show, " pathways | adj. p-value < 0.05"),
    x        = NULL,
    y        = "Gene Count"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "microarray_kegg_barchart.png"),
       plot = kegg_bar, width = 10, height = 7, dpi = 300)
cat("KEGG bar chart saved\n")

cat("\n=== Pathway Enrichment Summary (Microarray) ===\n")
cat("GO-BP terms  :", nrow(go_results), "\n")
cat("KEGG pathways:", nrow(kegg_results), "\n")
cat("Top GO-BP    :", go_results$Description[1], "\n")
cat("Top KEGG     :", kegg_results$Description[1], "\n")
cat("\nDone. Next: run ml_classification.R\n")

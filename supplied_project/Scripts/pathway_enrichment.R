# =============================================================================
# pathway_enrichment.R
# Project: Cross-Platform Transcriptomic Analysis of Alzheimer's Disease
# Dataset: GSE104704 (RNA-seq)
# Author:  Suraj Rathi | SynBioAI Lab, Atria University
#
# Description:
#   Performs Gene Ontology (GO) Biological Process and KEGG pathway enrichment
#   analysis on significant DEGs identified by deseq2_pipeline.R.
#   Produces dot plots for both analyses.
#
# Prerequisites:
#   Run deseq2_pipeline.R first (requires results/rnaseq_deseq2_sig.csv)
#
# Output:
#   - figures/rnaseq_go_bp_dotplot.png   GO-BP enrichment dot plot
#   - figures/rnaseq_kegg_dotplot.png    KEGG enrichment dot plot
#   - results/rnaseq_go_results.csv      Full GO-BP results table
#   - results/rnaseq_kegg_results.csv    Full KEGG results table
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

# ── 1. Configuration ──────────────────────────────────────────────────────
PADJ_CUTOFF <- 0.05
MIN_GS      <- 10    # minimum gene set size
MAX_GS      <- 500   # maximum gene set size
TOP_N       <- 20    # top terms to show in plots
OUT_DIR     <- "figures"
RES_DIR     <- "results"

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(RES_DIR, showWarnings = FALSE)

# ── 2. Load DEG Results ───────────────────────────────────────────────────
sig_df <- read.csv(file.path(RES_DIR, "rnaseq_deseq2_sig.csv"),
                   stringsAsFactors = FALSE)

cat("Loaded", nrow(sig_df), "significant DEGs\n")

# ── 3. Gene ID Conversion (Symbol → Entrez) ───────────────────────────────
# clusterProfiler requires Entrez IDs for KEGG enrichment
gene_symbols <- unique(sig_df$gene_id)

id_map <- bitr(
  geneID    = gene_symbols,
  fromType  = "SYMBOL",
  toType    = c("ENTREZID", "ENSEMBL"),
  OrgDb     = org.Hs.eg.db
)

cat("Mapped", nrow(id_map), "out of", length(gene_symbols), "gene symbols to Entrez IDs\n")

entrez_ids <- unique(id_map$ENTREZID)

# ── 4. GO Biological Process Enrichment ───────────────────────────────────
cat("\nRunning GO Biological Process enrichment...\n")

ego <- enrichGO(
  gene         = entrez_ids,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",             # Biological Process
  pAdjustMethod= "BH",
  pvalueCutoff = PADJ_CUTOFF,
  qvalueCutoff = 0.2,
  minGSSize    = MIN_GS,
  maxGSSize    = MAX_GS,
  readable     = TRUE              # convert Entrez back to gene symbols
)

# Simplify to remove redundant GO terms
ego_simplified <- simplify(ego,
                           cutoff     = 0.7,
                           by         = "p.adjust",
                           select_fun = min)

go_results <- as.data.frame(ego_simplified)
write.csv(go_results,
          file = file.path(RES_DIR, "rnaseq_go_results.csv"),
          row.names = FALSE)

cat("GO-BP enrichment complete:", nrow(go_results), "significant terms\n")

# GO Dot Plot
go_plot <- dotplot(ego_simplified,
                   showCategory = TOP_N,
                   font.size    = 10) +
  labs(
    title    = "GO Biological Process Enrichment — GSE104704 (RNA-seq)",
    subtitle = paste0("Top ", TOP_N, " terms | Adjusted p-value < 0.05")
  ) +
  theme(plot.title = element_text(face = "bold", size = 11))

ggsave(file.path(OUT_DIR, "rnaseq_go_bp_dotplot.png"),
       plot   = go_plot,
       width  = 10, height = 8, dpi = 300)
cat("GO dot plot saved\n")

# ── 5. KEGG Pathway Enrichment ────────────────────────────────────────────
cat("\nRunning KEGG pathway enrichment...\n")

ekegg <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",           # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = PADJ_CUTOFF,
  qvalueCutoff  = 0.2,
  minGSSize     = MIN_GS,
  maxGSSize     = MAX_GS
)

# Convert Entrez IDs to readable gene symbols in results
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

kegg_results <- as.data.frame(ekegg)
write.csv(kegg_results,
          file = file.path(RES_DIR, "rnaseq_kegg_results.csv"),
          row.names = FALSE)

cat("KEGG enrichment complete:", nrow(kegg_results), "significant pathways\n")

# KEGG Dot Plot
kegg_plot <- dotplot(ekegg,
                     showCategory = min(TOP_N, nrow(kegg_results)),
                     font.size    = 10) +
  labs(
    title    = "KEGG Pathway Enrichment — GSE104704 (RNA-seq)",
    subtitle = paste0("Top ", min(TOP_N, nrow(kegg_results)),
                      " pathways | Adjusted p-value < 0.05")
  ) +
  theme(plot.title = element_text(face = "bold", size = 11))

ggsave(file.path(OUT_DIR, "rnaseq_kegg_dotplot.png"),
       plot   = kegg_plot,
       width  = 10, height = 7, dpi = 300)
cat("KEGG dot plot saved\n")

# ── 6. Summary ────────────────────────────────────────────────────────────
cat("\n=== Pathway Enrichment Summary (RNA-seq) ===\n")
cat("GO-BP significant terms :", nrow(go_results),   "\n")
cat("KEGG significant pathways:", nrow(kegg_results), "\n")
cat("Top GO-BP term           :", go_results$Description[1], "\n")
cat("Top KEGG pathway         :", kegg_results$Description[1], "\n")
cat("\nDone. Next: run ml_classification.R\n")

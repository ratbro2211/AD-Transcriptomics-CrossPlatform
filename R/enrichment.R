# Use an assay-specific tested-gene background, not every annotated human gene.
run_enrichment <- function(platform) {
  require_packages(c("clusterProfiler", "org.Hs.eg.db"))
  paths <- output_dirs(platform)
  all <- read.csv(file.path(paths$results, "de_all.csv"))
  sig <- read.csv(file.path(paths$results, "de_significant.csv"))
  symbols <- function(x) unique(x[!is.na(x) & nzchar(x)])
  tested <- symbols(all$gene_symbol[!is.na(all$padj)])
  hits <- intersect(symbols(sig$gene_symbol), tested)
  if (!length(tested)) stop("No tested gene symbols available for the enrichment universe.")
  map <- suppressMessages(clusterProfiler::bitr(tested, fromType = "SYMBOL",
                                               toType = "ENTREZID",
                                               OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  universe <- unique(map$ENTREZID)
  genes <- unique(map$ENTREZID[map$SYMBOL %in% hits])
  write.csv(map, file.path(paths$results, "enrichment_id_mapping.csv"), row.names = FALSE)
  if (length(genes) < 2) {
    writeLines("Too few mapped significant genes for enrichment; no terms were tested.",
               file.path(paths$results, "enrichment_status.txt"))
    return(invisible(NULL))
  }
  go <- clusterProfiler::enrichGO(
    gene = genes, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = config$padj_cutoff, qvalueCutoff = 0.2)
  kegg <- clusterProfiler::enrichKEGG(
    gene = genes, universe = universe, organism = "hsa", pAdjustMethod = "BH",
    pvalueCutoff = config$padj_cutoff, qvalueCutoff = 0.2)
  # Empty results are legitimate; do not fabricate a plot or enrichment claim.
  write.csv(as.data.frame(go), file.path(paths$results, "GO_enrichment.csv"), row.names = FALSE)
  write.csv(as.data.frame(kegg), file.path(paths$results, "KEGG_enrichment.csv"), row.names = FALSE)
  save_session(paths$results)
}

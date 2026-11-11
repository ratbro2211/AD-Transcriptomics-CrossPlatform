options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran <- c("ashr", "randomForest", "glmnet", "e1071", "pROC")
bioc <- c("GEOquery", "Biobase", "DESeq2", "limma", "AnnotationDbi",
          "org.Hs.eg.db", "illuminaHumanv4.db", "clusterProfiler")
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
stopifnot(all(vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)))
message("Dependencies available. Record sessionInfo() for each actual analysis.")

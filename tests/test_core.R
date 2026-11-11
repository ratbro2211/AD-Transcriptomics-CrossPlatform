# No CRAN/Bioconductor packages required. Run from the repository root.
source("R/common.R")
source("R/ml_helpers.R")
source("R/comparison_helpers.R")
passed <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAIL: ", label)
  passed <<- passed + 1L
  message("PASS: ", label)
}
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")

rfiles <- list.files(".", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
for (file in rfiles) parse(file)
check(TRUE, paste("All", length(rfiles), "R files parse"))

counts <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
check(isTRUE(validate_counts(counts)), "Accept valid integer matrix")
fractional <- counts
fractional[1, 1] <- 1.1
check(fails(validate_counts(fractional)), "Reject fractional values instead of rounding")
negative <- counts
negative[1, 1] <- -1
check(fails(validate_counts(negative)), "Reject negative counts")
bad <- counts
colnames(bad) <- c("s1", "s1")
check(fails(validate_counts(bad)), "Reject duplicate sample IDs")
check(fails(validate_groups(c("AD", "AD"), c(Control = 1, AD = 1))), "Reject wrong group counts")
check(identical(levels(validate_groups(c("AD", "Control"), c(Control = 1, AD = 1))),
                c("Control", "AD")), "Control is reference level")

set.seed(42)
train <- matrix(rnorm(12 * 8), 12, dimnames = list(paste0("s", 1:12), paste0("g", 1:8)))
y <- factor(rep(c("Control", "AD"), each = 6), levels = c("Control", "AD"))
train[, 1] <- as.integer(y == "AD") * 5 + rnorm(12, sd = 0.1)
test <- matrix(rnorm(3 * 8), 3, dimnames = list(NULL, colnames(train)))
a <- prepare_features(train, y, test, 3)
b <- prepare_features(train, y, test * 10000 + 999, 3)
check(identical(a$genes, b$genes) && identical(a$center, b$center) &&
        identical(a$scale, b$scale), "Held-out features cannot alter fitted preprocessing")
check("g1" %in% a$genes, "Training-only ranking detects synthetic signal")
check(max(abs(colMeans(a$train))) < 1e-10, "Training features are centered")
check(max(abs(apply(a$train, 2, sd) - 1)) < 1e-10, "Training features are scaled")
folds <- stratified_folds(y, 3)
check(all(table(folds, y) == 2), "Inner folds are stratified")

d1 <- data.frame(gene_symbol = c("A", "B", "C", "D"),
                 log2FoldChange = c(2, -1, 0.1, 1), padj = c(0.01, 0.01, 0.8, 0.02))
d2 <- data.frame(gene_symbol = c("A", "B", "C", "D"),
                 log2FoldChange = c(-1, 1, 0.2, 1), padj = c(0.01, 0.01, 0.8, 0.02))
comp <- compare_de(d1, d2, 0.05, 0.5)
check(nrow(comp$shared) == 4 && nrow(comp$overlap) == 3, "Separate common-tested and overlapping-DE scopes")
check(sum(comp$overlap$concordant) == 1, "Compute direction concordance from signs")
check(comp$summary$value[comp$summary$metric == "Pearson_r_shared_DEGs"] ==
        cor(comp$overlap$rnaseq_LFC, comp$overlap$array_logFC), "Correlation uses declared gene scope")

historical <- read.csv("results/overlap_DE_genes.csv")
check(nrow(historical) == 47 && !anyDuplicated(historical$Gene), "Supplied overlap has 47 unique genes")
same <- sign(historical$rnaseq_LFC) == sign(historical$array_logFC)
check(sum(same) == 1 && historical$Gene[same] == "ANGPTL4", "Supplied overlap has one concordant gene")
check(abs(cor(historical$rnaseq_LFC, historical$array_logFC) + 0.8722673) < 1e-6,
      "Supplied 47-row correlation is about -0.872, not full-gene -0.438")
message("Completed ", passed, " checks. No biological/ML end-to-end execution is implied.")

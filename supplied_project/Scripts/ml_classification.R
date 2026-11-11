# =============================================================================
# ml_classification.R
# Project: Cross-Platform Transcriptomic Analysis of Alzheimer's Disease
# Dataset: GSE104704 (RNA-seq, n = 22)
# Author:  Suraj Rathi | SynBioAI Lab, Atria University
#
# Description:
#   Builds and evaluates five machine learning classifiers to distinguish
#   Alzheimer's Disease from Control samples using DEG-derived features.
#   Uses Leave-One-Out Cross-Validation (LOOCV) due to small sample size (n=22).
#   Models: LASSO Logistic Regression, Random Forest, SVM, KNN, Naive Bayes.
#
# Prerequisites:
#   Run deseq2_pipeline.R first (requires RDS/dds.rds, results/rnaseq_deseq2_sig.csv)
#
# Output:
#   - figures/rnaseq_ml_comparison.png       Accuracy + AUC bar chart
#   - figures/rnaseq_rf_importance.png       Random Forest top-20 gene importance
#   - results/rnaseq_ml_metrics.csv          Per-model metrics table
#   - RDS/rnaseq_ml_results.rds             Full results list (for cross-platform)
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(DESeq2)
  library(caret)
  library(randomForest)
  library(glmnet)
  library(e1071)
  library(class)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ── 1. Configuration ──────────────────────────────────────────────────────
TOP_FEATURES <- 100   # number of top DEGs to use as features
OUT_DIR      <- "figures"
RES_DIR      <- "results"
RDS_DIR      <- "RDS"

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(RES_DIR, showWarnings = FALSE)
dir.create(RDS_DIR, showWarnings = FALSE)

set.seed(42)

# ── 2. Prepare Feature Matrix ─────────────────────────────────────────────
# Load DESeq2 object and significant DEGs
dds    <- readRDS(file.path(RDS_DIR, "dds.rds"))
sig_df <- read.csv(file.path(RES_DIR, "rnaseq_deseq2_sig.csv"),
                   stringsAsFactors = FALSE)

# Select top features by adjusted p-value
top_genes <- sig_df %>%
  arrange(padj) %>%
  slice_head(n = TOP_FEATURES) %>%
  pull(gene_id)

cat("Using", length(top_genes), "features for ML\n")

# Extract variance-stabilised counts (VST) as features
vsd <- vst(dds, blind = FALSE)
expr_mat <- t(assay(vsd)[top_genes, ])   # samples × genes

# Build labels
labels <- factor(dds$condition, levels = c("Control", "AD"))

# Combine into data frame
ml_df <- as.data.frame(expr_mat)
ml_df$label <- labels

cat("Feature matrix dimensions:", nrow(ml_df), "samples ×",
    ncol(ml_df) - 1, "features\n")

# ── 3. LOOCV Helper Function ──────────────────────────────────────────────
# Returns list: predictions, actual labels, probabilities
run_loocv <- function(X, y, model_fn) {
  n      <- nrow(X)
  preds  <- character(n)
  probs  <- numeric(n)

  for (i in seq_len(n)) {
    X_train <- X[-i, , drop = FALSE]
    y_train <- y[-i]
    X_test  <- X[i, , drop = FALSE]

    result  <- model_fn(X_train, y_train, X_test)
    preds[i] <- result$pred
    probs[i] <- result$prob
  }

  list(
    preds  = factor(preds, levels = levels(y)),
    actual = y,
    probs  = probs
  )
}

# ── 4. Compute Metrics ────────────────────────────────────────────────────
compute_metrics <- function(loocv_result, model_name) {
  cm       <- confusionMatrix(loocv_result$preds, loocv_result$actual,
                              positive = "AD")
  roc_obj  <- pROC::roc(as.numeric(loocv_result$actual == "AD"),
                        loocv_result$probs, quiet = TRUE)
  data.frame(
    Model    = model_name,
    Accuracy = round(cm$overall["Accuracy"], 3),
    Sensitivity = round(cm$byClass["Sensitivity"], 3),
    Specificity = round(cm$byClass["Specificity"], 3),
    AUC      = round(as.numeric(pROC::auc(roc_obj)), 3)
  )
}

# install pROC if missing
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

X <- as.matrix(ml_df[, -ncol(ml_df)])
y <- ml_df$label

# ── 5. Model Definitions ──────────────────────────────────────────────────

# --- 5a. LASSO Logistic Regression ---
lasso_fn <- function(X_train, y_train, X_test) {
  cv_fit <- cv.glmnet(X_train, as.numeric(y_train) - 1,
                      family = "binomial", alpha = 1, nfolds = 5)
  prob   <- predict(cv_fit, newx = X_test, s = "lambda.min",
                    type = "response")[1]
  pred   <- ifelse(prob > 0.5, "AD", "Control")
  list(pred = pred, prob = prob)
}

# --- 5b. Random Forest ---
rf_fn <- function(X_train, y_train, X_test) {
  rf   <- randomForest(x = X_train, y = y_train, ntree = 500)
  prob <- predict(rf, X_test, type = "prob")[, "AD"]
  pred <- as.character(predict(rf, X_test))
  list(pred = pred, prob = prob)
}

# --- 5c. SVM (RBF kernel) ---
svm_fn <- function(X_train, y_train, X_test) {
  sv   <- svm(X_train, y_train, kernel = "radial",
              probability = TRUE, scale = TRUE)
  prd  <- predict(sv, X_test, probability = TRUE)
  prob <- attr(prd, "probabilities")[, "AD"]
  list(pred = as.character(prd), prob = prob)
}

# --- 5d. KNN (k = 5) ---
knn_fn <- function(X_train, y_train, X_test) {
  pred <- knn(X_train, X_test, y_train, k = 5, prob = TRUE)
  prob <- attr(pred, "prob")
  # knn returns prob of winning class; adjust if predicted is "Control"
  if (as.character(pred) == "Control") prob <- 1 - prob
  list(pred = as.character(pred), prob = prob)
}

# --- 5e. Naive Bayes ---
nb_fn <- function(X_train, y_train, X_test) {
  nb   <- naiveBayes(X_train, y_train)
  prob <- predict(nb, X_test, type = "raw")[, "AD"]
  pred <- predict(nb, X_test)
  list(pred = as.character(pred), prob = prob)
}

# ── 6. Run LOOCV for All Models ───────────────────────────────────────────
cat("\nRunning LOOCV — this may take a few minutes...\n")

models <- list(
  "LASSO"         = lasso_fn,
  "Random Forest" = rf_fn,
  "SVM"           = svm_fn,
  "KNN"           = knn_fn,
  "Naive Bayes"   = nb_fn
)

all_results <- list()
metrics_list <- list()

for (model_name in names(models)) {
  cat("  Running:", model_name, "...\n")
  res <- run_loocv(X, y, models[[model_name]])
  all_results[[model_name]] <- res
  metrics_list[[model_name]] <- compute_metrics(res, model_name)
}

metrics_df <- bind_rows(metrics_list)
metrics_df$Dataset <- "GSE104704 (RNA-seq, LOOCV)"

write.csv(metrics_df,
          file = file.path(RES_DIR, "rnaseq_ml_metrics.csv"),
          row.names = FALSE)

cat("\n=== ML Performance Summary (RNA-seq, LOOCV) ===\n")
print(metrics_df)

saveRDS(list(results = all_results, metrics = metrics_df, genes = top_genes),
        file = file.path(RDS_DIR, "rnaseq_ml_results.rds"))

# ── 7. Bar Chart: Accuracy + AUC ─────────────────────────────────────────
plot_df <- metrics_df %>%
  tidyr::pivot_longer(cols = c(Accuracy, AUC),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(Model = factor(Model, levels = metrics_df$Model[order(metrics_df$AUC)]))

ml_bar <- ggplot(plot_df, aes(x = Model, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.2f", Value)),
            position = position_dodge(0.7),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(Accuracy = "#4E79A7", AUC = "#F28E2B")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  coord_flip() +
  labs(
    title    = "ML Classification Performance — GSE104704 (RNA-seq)",
    subtitle = "Leave-One-Out Cross-Validation (LOOCV) | n = 22",
    x        = NULL,
    y        = "Score",
    fill     = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(OUT_DIR, "rnaseq_ml_comparison.png"),
       plot = ml_bar, width = 8, height = 5, dpi = 300)

# ── 8. Random Forest Feature Importance ───────────────────────────────────
cat("\nFitting full Random Forest for feature importance...\n")
rf_full <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)

imp_df <- as.data.frame(importance(rf_full)) %>%
  tibble::rownames_to_column("Gene") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n = 20)

rf_imp_plot <- ggplot(imp_df,
                      aes(x = reorder(Gene, MeanDecreaseGini),
                          y = MeanDecreaseGini)) +
  geom_col(fill = "#2C7BB6", alpha = 0.85) +
  coord_flip() +
  labs(
    title    = "Random Forest — Top 20 Gene Importance",
    subtitle = "GSE104704 RNA-seq | Mean Decrease in Gini Impurity",
    x        = "Gene",
    y        = "Mean Decrease Gini"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "rnaseq_rf_importance.png"),
       plot = rf_imp_plot, width = 7, height = 6, dpi = 300)

cat("\nAll figures and results saved.\n")
cat("Best model by AUC:", metrics_df$Model[which.max(metrics_df$AUC)],
    "— AUC =", max(metrics_df$AUC), "\n")
cat("\nDone. RNA-seq analysis complete.\n")

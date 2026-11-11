# =============================================================================
# ml_classification.R
# Project: Cross-Platform Transcriptomic Analysis of Alzheimer's Disease
# Dataset: GSE132903 (Microarray, n = 195)
# Author:  Suraj Rathi | SynBioAI Lab, Atria University
#
# Description:
#   Machine learning classification of AD vs. Control using microarray DEG
#   features. Uses stratified 80/20 train-test split (viable due to n = 195).
#   Models: LASSO Logistic Regression, Random Forest, SVM, KNN, Naive Bayes.
#
# Prerequisites:
#   Run limma_pipeline.R first (requires RDS/eset_norm.rds, results/microarray_limma_sig.csv)
#
# Output:
#   - figures/microarray_ml_comparison.png
#   - figures/microarray_rf_importance.png
#   - results/microarray_ml_metrics.csv
#   - RDS/microarray_ml_results.rds
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Biobase)
  library(caret)
  library(randomForest)
  library(glmnet)
  library(e1071)
  library(class)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ── 1. Configuration ──────────────────────────────────────────────────────
TOP_FEATURES <- 100    # top DEGs by adj p-value
TRAIN_FRAC   <- 0.80   # 80/20 split
N_FOLDS      <- 5      # for internal LASSO tuning
OUT_DIR      <- "figures"
RES_DIR      <- "results"
RDS_DIR      <- "RDS"

dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(RES_DIR, showWarnings = FALSE)
dir.create(RDS_DIR, showWarnings = FALSE)

set.seed(42)

# ── 2. Prepare Feature Matrix ─────────────────────────────────────────────
eset   <- readRDS(file.path(RDS_DIR, "eset_norm.rds"))
sig_df <- read.csv(file.path(RES_DIR, "microarray_limma_sig.csv"),
                   stringsAsFactors = FALSE)

# Select top features by adj.P.Val
top_probes <- sig_df %>%
  arrange(adj.P.Val) %>%
  slice_head(n = TOP_FEATURES) %>%
  pull(probe_id)

# Keep only probes present in the ExpressionSet
top_probes <- intersect(top_probes, rownames(eset))
cat("Using", length(top_probes), "probe features for ML\n")

# Feature matrix: samples × probes
expr_mat <- t(exprs(eset)[top_probes, ])
labels   <- factor(eset$condition, levels = c("Control", "AD"))

ml_df       <- as.data.frame(expr_mat)
ml_df$label <- labels

cat("Dataset size:", nrow(ml_df), "samples ×", length(top_probes), "features\n")

# ── 3. Stratified 80/20 Train-Test Split ──────────────────────────────────
train_idx  <- createDataPartition(labels, p = TRAIN_FRAC, list = FALSE)
X_train    <- as.matrix(ml_df[train_idx,  -ncol(ml_df)])
X_test     <- as.matrix(ml_df[-train_idx, -ncol(ml_df)])
y_train    <- ml_df$label[train_idx]
y_test     <- ml_df$label[-train_idx]

cat("\nTrain set:", nrow(X_train), "samples | Test set:", nrow(X_test), "samples\n")
cat("Train class distribution:\n"); print(table(y_train))
cat("Test class distribution:\n");  print(table(y_test))

# ── 4. Model Training & Evaluation ───────────────────────────────────────
evaluate_model <- function(model_name, pred_class, pred_prob) {
  pred_class <- factor(pred_class, levels = levels(y_test))
  cm         <- confusionMatrix(pred_class, y_test, positive = "AD")
  roc_obj    <- roc(as.numeric(y_test == "AD"), pred_prob, quiet = TRUE)
  data.frame(
    Model       = model_name,
    Accuracy    = round(cm$overall["Accuracy"], 3),
    Sensitivity = round(cm$byClass["Sensitivity"], 3),
    Specificity = round(cm$byClass["Specificity"], 3),
    AUC         = round(as.numeric(auc(roc_obj)), 3)
  )
}

metrics_list <- list()

# --- LASSO Logistic Regression ---
cat("\nTraining LASSO...\n")
cv_lasso    <- cv.glmnet(X_train, as.numeric(y_train) - 1,
                         family = "binomial", alpha = 1, nfolds = N_FOLDS)
lasso_prob  <- predict(cv_lasso, newx = X_test, s = "lambda.min",
                       type = "response")[, 1]
lasso_pred  <- ifelse(lasso_prob > 0.5, "AD", "Control")
metrics_list[["LASSO"]] <- evaluate_model("LASSO", lasso_pred, lasso_prob)

# --- Random Forest ---
cat("Training Random Forest...\n")
rf_model   <- randomForest(x = X_train, y = y_train, ntree = 500)
rf_prob    <- predict(rf_model, X_test, type = "prob")[, "AD"]
rf_pred    <- predict(rf_model, X_test)
metrics_list[["Random Forest"]] <- evaluate_model("Random Forest",
                                                   as.character(rf_pred), rf_prob)

# --- SVM (RBF kernel) ---
cat("Training SVM...\n")
svm_model  <- svm(X_train, y_train, kernel = "radial",
                  probability = TRUE, scale = TRUE)
svm_pred   <- predict(svm_model, X_test, probability = TRUE)
svm_prob   <- attr(svm_pred, "probabilities")[, "AD"]
metrics_list[["SVM"]] <- evaluate_model("SVM",
                                        as.character(svm_pred), svm_prob)

# --- KNN (k = 7 — slightly larger k better for n=195) ---
cat("Training KNN...\n")
knn_pred <- knn(X_train, X_test, y_train, k = 7, prob = TRUE)
knn_prob <- attr(knn_pred, "prob")
knn_prob[as.character(knn_pred) == "Control"] <-
  1 - knn_prob[as.character(knn_pred) == "Control"]
metrics_list[["KNN"]] <- evaluate_model("KNN",
                                        as.character(knn_pred), knn_prob)

# --- Naive Bayes ---
cat("Training Naive Bayes...\n")
nb_model  <- naiveBayes(X_train, y_train)
nb_prob   <- predict(nb_model, X_test, type = "raw")[, "AD"]
nb_pred   <- predict(nb_model, X_test)
metrics_list[["Naive Bayes"]] <- evaluate_model("Naive Bayes",
                                                 as.character(nb_pred), nb_prob)

# ── 5. Compile & Save Metrics ─────────────────────────────────────────────
metrics_df           <- bind_rows(metrics_list)
metrics_df$Dataset   <- "GSE132903 (Microarray, 80/20 split)"

write.csv(metrics_df,
          file = file.path(RES_DIR, "microarray_ml_metrics.csv"),
          row.names = FALSE)

cat("\n=== ML Performance Summary (Microarray, 80/20 Split) ===\n")
print(metrics_df)

saveRDS(list(rf_model = rf_model, metrics = metrics_df, probes = top_probes),
        file = file.path(RDS_DIR, "microarray_ml_results.rds"))

# ── 6. Bar Chart: Accuracy + AUC ─────────────────────────────────────────
plot_df <- metrics_df %>%
  pivot_longer(cols = c(Accuracy, AUC),
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
    title    = "ML Classification Performance — GSE132903 (Microarray)",
    subtitle = "Stratified 80/20 Train-Test Split | n = 195",
    x        = NULL, y = "Score", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(file.path(OUT_DIR, "microarray_ml_comparison.png"),
       plot = ml_bar, width = 8, height = 5, dpi = 300)

# ── 7. Random Forest Feature Importance ───────────────────────────────────
cat("\nFitting full Random Forest for feature importance...\n")
rf_full <- randomForest(x = as.matrix(ml_df[, -ncol(ml_df)]),
                        y = ml_df$label,
                        ntree = 1000, importance = TRUE)

imp_df <- as.data.frame(importance(rf_full)) %>%
  tibble::rownames_to_column("Probe") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n = 20)

# Replace probe IDs with gene symbols where available
symbol_lookup <- sig_df %>%
  select(probe_id, gene_symbol) %>%
  distinct() %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

imp_df <- imp_df %>%
  left_join(symbol_lookup, by = c("Probe" = "probe_id")) %>%
  mutate(Label = ifelse(is.na(gene_symbol), Probe, gene_symbol))

rf_imp_plot <- ggplot(imp_df,
                      aes(x = reorder(Label, MeanDecreaseGini),
                          y = MeanDecreaseGini)) +
  geom_col(fill = "#D73027", alpha = 0.85) +
  coord_flip() +
  labs(
    title    = "Random Forest — Top 20 Probe Importance",
    subtitle = "GSE132903 Microarray | Mean Decrease in Gini Impurity",
    x        = "Gene / Probe",
    y        = "Mean Decrease Gini"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "microarray_rf_importance.png"),
       plot = rf_imp_plot, width = 7, height = 6, dpi = 300)

cat("\nAll figures and results saved.\n")
cat("Best model by AUC:", metrics_df$Model[which.max(metrics_df$AUC)],
    "— AUC =", max(metrics_df$AUC), "\n")
cat("\nDone. Microarray analysis complete.\n")

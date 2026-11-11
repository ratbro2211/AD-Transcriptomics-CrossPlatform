# No full-cohort DEG list is accepted by this module.
# Inputs: all eligible expression features and diagnosis labels.
# Feature ranking, variance filtering, centering and scaling are training-only.

prepare_features <- function(train, y, test, top_n = 100L) {
  train <- as.matrix(train)
  test <- as.matrix(test)
  if (!identical(colnames(train), colnames(test))) stop("Feature columns differ.")
  if (any(!is.finite(train)) || any(!is.finite(test))) stop("Non-finite ML input.")
  if (length(y) != nrow(train) || any(table(y) < 2)) stop("Insufficient training groups.")
  means <- colMeans(train)
  spread <- apply(train, 2, sd)
  keep <- is.finite(spread) & spread > 0
  if (sum(keep) < 2) stop("Fewer than two variable training features.")
  train <- train[, keep, drop = FALSE]
  test <- test[, keep, drop = FALSE]
  means <- means[keep]
  spread <- spread[keep]
  # Welch-style training-only t-statistic for ranking, not a DE significance test.
  a <- train[y == "AD", , drop = FALSE]
  b <- train[y == "Control", , drop = FALSE]
  variance <- apply(a, 2, var) / nrow(a) + apply(b, 2, var) / nrow(b)
  score <- abs(colMeans(a) - colMeans(b)) / sqrt(pmax(variance, 1e-12))
  idx <- order(-score, colnames(train))[seq_len(min(top_n, ncol(train)))]
  chosen <- colnames(train)[idx]
  center <- means[idx]
  scale <- spread[idx]
  transform <- function(x) sweep(sweep(x[, chosen, drop = FALSE], 2, center, "-"), 2, scale, "/")
  list(train = transform(train), test = transform(test), genes = chosen,
       center = center, scale = scale)
}

stratified_folds <- function(y, k = 5L) {
  k <- min(k, min(table(y)))
  if (k < 2) stop("Too few samples for nested tuning.")
  out <- integer(length(y))
  for (level in levels(y)) {
    ids <- which(y == level)
    out[ids] <- sample(rep(seq_len(k), length.out = length(ids)))
  }
  out
}

tune_lasso <- function(x, y, top_n) {
  # Feature ranking and scaling are repeated within EVERY inner training fold.
  # cv.glmnet on an already-selected full outer-training matrix would not
  # evaluate that feature-selection operation inside its inner validation.
  grid <- 10^seq(-3, 0, length.out = 12)
  folds <- stratified_folds(y)
  losses <- matrix(NA_real_, nrow = length(y), ncol = length(grid))
  for (fold in unique(folds)) {
    tr <- folds != fold
    f <- prepare_features(x[tr, , drop = FALSE], y[tr],
                          x[!tr, , drop = FALSE], top_n)
    fit <- glmnet::glmnet(f$train, as.integer(y[tr] == "AD"),
                          family = "binomial", alpha = 1,
                          lambda = sort(grid, decreasing = TRUE), standardize = FALSE)
    p <- predict(fit, f$test, s = grid, type = "response")
    p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
    actual <- as.integer(y[!tr] == "AD")
    losses[!tr, ] <- -(actual * log(p) + (1 - actual) * log(1 - p))
  }
  grid[which.min(colMeans(losses))]
}

model_probabilities <- function(train, y, test, lambda) {
  response <- as.integer(y == "AD")
  lasso <- glmnet::glmnet(train, response, family = "binomial", alpha = 1,
                          lambda = lambda, standardize = FALSE)
  rf <- randomForest::randomForest(x = train, y = y, ntree = 500)
  svm <- e1071::svm(x = train, y = y, kernel = "radial",
                    probability = TRUE, scale = FALSE)
  svm_pred <- predict(svm, test, probability = TRUE)
  knn <- class::knn(train, test, y, k = min(5L, nrow(train)), prob = TRUE)
  knn_p <- attr(knn, "prob")
  knn_p[as.character(knn) != "AD"] <- 1 - knn_p[as.character(knn) != "AD"]
  # Gaussian NB needs positive within-class variance for every feature.
  # Smooth its variance estimates rather than permitting division by zero.
  nb <- e1071::naiveBayes(train, y)
  nb$tables <- lapply(nb$tables, function(z) {
    z[, 2] <- pmax(z[, 2], 1e-6)
    z
  })
  list(
    LASSO = as.numeric(predict(lasso, test, s = lambda, type = "response")),
    `Random Forest` = as.numeric(predict(rf, test, type = "prob")[, "AD"]),
    SVM = as.numeric(attr(svm_pred, "probabilities")[, "AD"]),
    KNN = as.numeric(knn_p),
    `Naive Bayes` = as.numeric(predict(nb, test, type = "raw")[, "AD"])
  )
}

score_predictions <- function(pred) {
  rows <- lapply(split(pred, pred$Model), function(d) {
    if (any(!is.finite(d$probability))) stop("Invalid prediction probability.")
    actual <- d$actual == "AD"
    called <- d$probability >= 0.5
    correct <- sum(actual == called)
    ci <- binom.test(correct, nrow(d))$conf.int
    # Fix direction: a high AD probability must mean AD, not auto-invert AUC.
    roc <- pROC::roc(factor(d$actual, levels = c("Control", "AD")),
                     d$probability, levels = c("Control", "AD"),
                     direction = "<", quiet = TRUE)
    data.frame(Model = d$Model[1], n_evaluated = nrow(d),
               Accuracy = correct / nrow(d),
               Accuracy_CI_low = ci[1], Accuracy_CI_high = ci[2],
               Sensitivity = mean(called[actual]),
               Specificity = mean(!called[!actual]),
               AUC = as.numeric(pROC::auc(roc)))
  })
  do.call(rbind, rows)
}

run_ml <- function(platform, kind) {
  require_packages(c("glmnet", "randomForest", "e1071", "class", "pROC"))
  paths <- output_dirs(platform)
  input <- readRDS(file.path(paths$cache, "ml_input.rds"))
  y <- factor(input$labels, levels = c("Control", "AD"))
  if (anyNA(y)) stop("Unknown diagnosis labels.")
  if (kind == "rnaseq") {
    validate_counts(input$raw_counts)
    counts <- input$raw_counts
    # Per-sample transformation: uses only that sample's library size.
    x <- t(log2(sweep(counts, 2, colSums(counts), "/") * 1e6 + 1))
  } else {
    x <- t(input$expression)
  }
  if (anyDuplicated(colnames(x)) || anyDuplicated(input$sample_id)) stop("Duplicate identifiers.")
  if (!identical(rownames(x), input$sample_id) || nrow(x) != length(y))
    stop("ML labels/sample IDs are not aligned.")
  set.seed(config$seed)
  if (kind == "rnaseq") {
    if (nrow(x) != 22L) stop("RNA-seq LOOCV expects exactly 22 AD/Old samples.")
    splits <- lapply(seq_len(nrow(x)), function(i) list(train = setdiff(seq_len(nrow(x)), i), test = i))
    validation <- "Outer LOOCV; inner 5-fold LASSO tuning"
  } else {
    train <- unlist(lapply(levels(y), function(z) {
      ids <- which(y == z)
      sample(ids, floor(length(ids) * config$train_fraction))
    }), use.names = FALSE)
    splits <- list(list(train = sort(train), test = setdiff(seq_len(nrow(x)), train)))
    validation <- "Stratified 80/20 holdout; inner 5-fold LASSO tuning"
  }
  rows <- list()
  selection <- list()
  for (i in seq_along(splits)) {
    tr <- splits[[i]]$train
    te <- splits[[i]]$test
    message("Outer split ", i, "/", length(splits))
    lambda <- tune_lasso(x[tr, , drop = FALSE], y[tr], config$top_features)
    f <- prepare_features(x[tr, , drop = FALSE], y[tr],
                          x[te, , drop = FALSE], config$top_features)
    predictions <- model_probabilities(f$train, y[tr], f$test, lambda)
    for (model in names(predictions)) {
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = input$sample_id[te], actual = as.character(y[te]),
        Model = model, probability = predictions[[model]], outer_split = i)
    }
    selection[[i]] <- data.frame(outer_split = i, gene_id = f$genes, lambda = lambda)
  }
  pred <- do.call(rbind, rows)
  metrics <- score_predictions(pred)
  metrics$validation <- validation
  metrics$status <- "Internal exploratory validation; not external clinical validation"
  write.csv(pred, file.path(paths$results, "ml_predictions.csv"), row.names = FALSE)
  write.csv(metrics, file.path(paths$results, "ml_metrics.csv"), row.names = FALSE)
  write.csv(do.call(rbind, selection), file.path(paths$results, "ml_fold_features.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(seq_along(splits), function(i) {
    data.frame(outer_split = i, sample_id = input$sample_id,
               role = ifelse(seq_along(y) %in% splits[[i]]$train, "train", "test"))
  })), file.path(paths$results, "ml_splits.csv"), row.names = FALSE)
  save_session(paths$results)
  print(metrics)
}

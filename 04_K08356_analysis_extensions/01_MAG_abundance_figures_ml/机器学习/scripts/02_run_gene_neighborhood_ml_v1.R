suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(vegan)
  library(rpart)
})
has_random_forest <- requireNamespace("randomForest", quietly = TRUE)
if (has_random_forest) {
  library(randomForest)
}

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm/机器学习"
model_dir <- file.path(base_dir, "01_gene_neighborhood_model")
fig_dir <- file.path(base_dir, "figures")
out_dir <- file.path(base_dir, "outputs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

feature_file <- file.path(model_dir, "gene_neighborhood_feature_matrix.csv")
feature <- read_csv(feature_file, show_col_types = FALSE)

branch_order <- c(
  "IdrA-associated",
  "canonical aioA-associated",
  "aioA-like-associated",
  "unknown AioA-like / uncertain DMSOR"
)
feature <- feature %>%
  mutate(branch = factor(branch, levels = branch_order))

binary_cols <- grep("^has_", names(feature), value = TRUE)
count_cols <- grep("^n_(Rieske|cyt|arsenic|transporter|nitrogen|sox|redox|annotated|hypothetical)", names(feature), value = TRUE)
X <- feature %>% select(all_of(binary_cols))
X_mat <- as.matrix(X)
rownames(X_mat) <- feature$target_gene_id
feature_variance <- apply(X_mat, 2, var)
zero_var_features <- names(feature_variance)[is.na(feature_variance) | feature_variance == 0]
if (length(zero_var_features) > 0) {
  write_csv(
    data.frame(feature = zero_var_features),
    file.path(model_dir, "zero_variance_binary_features_excluded_from_models.csv")
  )
  binary_cols <- setdiff(binary_cols, zero_var_features)
  X <- feature %>% select(all_of(binary_cols))
  X_mat <- as.matrix(X)
  rownames(X_mat) <- feature$target_gene_id
}

pretty_feature <- function(x) {
  x %>%
    sub("^has_", "", .) %>%
    gsub("_", " ", .)
}

# Heatmap of binary neighborhood features.
ann <- data.frame(
  branch = feature$branch,
  habitat = feature$habitat,
  row.names = feature$target_gene_id
)
pdf(file.path(fig_dir, "gene_neighborhood_binary_feature_heatmap.pdf"), width = 8.8, height = 8.5)
pheatmap(
  X_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = ann,
  color = c("white", "#2C7FB8"),
  breaks = c(-0.01, 0.5, 1.01),
  labels_col = pretty_feature(colnames(X_mat)),
  fontsize_row = 5.5,
  fontsize_col = 8,
  main = "K08356 gene-neighborhood binary feature heatmap"
)
dev.off()

# PCA.
pca <- prcomp(X_mat, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
  mutate(
    target_gene_id = feature$target_gene_id,
    branch = feature$branch,
    habitat = feature$habitat
  )
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
write_csv(pca_scores, file.path(model_dir, "gene_neighborhood_PCA_scores.csv"))

p_pca <- ggplot(pca_scores, aes(PC1, PC2, color = branch, shape = habitat)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_classic(base_family = "Times", base_size = 11) +
  labs(
    title = "PCA of K08356 gene-neighborhood features",
    x = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * pca_var[2]),
    color = "Branch",
    shape = "Habitat"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(fig_dir, "gene_neighborhood_PCA.pdf"), p_pca, width = 6.8, height = 5.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "gene_neighborhood_PCA.svg"), p_pca, width = 6.8, height = 5.2, device = svg)

# PERMANOVA on Jaccard distance for binary features.
set.seed(20260708)
jaccard <- vegdist(X_mat, method = "jaccard", binary = TRUE)
permanova <- adonis2(jaccard ~ branch, data = feature, permutations = 999)
permanova_tbl <- as.data.frame(permanova) %>%
  tibble::rownames_to_column("term")
write_csv(permanova_tbl, file.path(model_dir, "PERMANOVA_jaccard_binary_features.csv"))

# Leave-one-out nearest centroid classifier.
loocv_rows <- list()
for (i in seq_len(nrow(feature))) {
  train_idx <- setdiff(seq_len(nrow(feature)), i)
  train_x <- X_mat[train_idx, , drop = FALSE]
  test_x <- X_mat[i, , drop = FALSE]
  mu <- colMeans(train_x)
  sigma <- apply(train_x, 2, sd)
  sigma[is.na(sigma) | sigma == 0] <- 1
  train_z <- sweep(sweep(train_x, 2, mu, "-"), 2, sigma, "/")
  test_z <- sweep(sweep(test_x, 2, mu, "-"), 2, sigma, "/")
  train_branch <- droplevels(feature$branch[train_idx])
  centroids <- rowsum(train_z, group = train_branch) / as.vector(table(train_branch))
  d <- apply(centroids, 1, function(center) sqrt(sum((test_z[1, ] - center)^2)))
  pred <- names(which.min(d))
  loocv_rows[[i]] <- data.frame(
    target_gene_id = feature$target_gene_id[i],
    true_branch = as.character(feature$branch[i]),
    predicted_branch = pred,
    correct = pred == as.character(feature$branch[i]),
    min_distance = min(d),
    stringsAsFactors = FALSE
  )
}
loocv <- bind_rows(loocv_rows)
write_csv(loocv, file.path(model_dir, "LOOCV_nearest_centroid_predictions.csv"))

confusion <- loocv %>%
  count(true_branch, predicted_branch, name = "n") %>%
  complete(true_branch = branch_order, predicted_branch = branch_order, fill = list(n = 0))
write_csv(confusion, file.path(model_dir, "LOOCV_nearest_centroid_confusion_matrix_long.csv"))

overall_accuracy <- mean(loocv$correct)
branch_accuracy <- loocv %>%
  group_by(true_branch) %>%
  summarise(n = n(), accuracy = mean(correct), .groups = "drop")
balanced_accuracy <- mean(branch_accuracy$accuracy)
model_summary <- bind_rows(
  data.frame(metric = "overall_accuracy", value = overall_accuracy),
  data.frame(metric = "balanced_accuracy", value = balanced_accuracy),
  data.frame(metric = "n_targets", value = nrow(feature)),
  data.frame(metric = "n_binary_features", value = length(binary_cols))
)
write_csv(model_summary, file.path(model_dir, "gene_neighborhood_model_summary.csv"))
write_csv(branch_accuracy, file.path(model_dir, "LOOCV_branch_accuracy.csv"))

p_conf <- ggplot(confusion, aes(predicted_branch, true_branch, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), family = "Times", size = 3.5) +
  scale_fill_gradient(low = "white", high = "#2C7FB8") +
  theme_classic(base_family = "Times", base_size = 10) +
  labs(
    title = "LOOCV nearest-centroid confusion matrix",
    x = "Predicted branch",
    y = "True branch",
    fill = "Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
ggsave(file.path(fig_dir, "LOOCV_nearest_centroid_confusion_matrix.pdf"), p_conf, width = 6.5, height = 5.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "LOOCV_nearest_centroid_confusion_matrix.svg"), p_conf, width = 6.5, height = 5.2, device = svg)

# Random Forest validation when the package is available.
if (has_random_forest) {
  set.seed(20260708)
  rf_df <- feature %>% select(branch, all_of(binary_cols))
  rf_model <- randomForest(
    branch ~ .,
    data = rf_df,
    ntree = 2000,
    importance = TRUE,
    proximity = TRUE
  )
  rf_oob_conf <- as.data.frame.matrix(rf_model$confusion[, setdiff(colnames(rf_model$confusion), "class.error"), drop = FALSE])
  rf_oob_conf <- rf_oob_conf %>%
    tibble::rownames_to_column("true_branch")
  write_csv(rf_oob_conf, file.path(model_dir, "RF_OOB_confusion_matrix.csv"))

  rf_summary <- data.frame(
    metric = c("OOB_error", "OOB_accuracy", "ntree", "mtry", "n_targets", "n_binary_features"),
    value = c(
      rf_model$err.rate[nrow(rf_model$err.rate), "OOB"],
      1 - rf_model$err.rate[nrow(rf_model$err.rate), "OOB"],
      rf_model$ntree,
      rf_model$mtry,
      nrow(feature),
      length(binary_cols)
    )
  )
  write_csv(rf_summary, file.path(model_dir, "RF_model_summary.csv"))

  rf_imp <- as.data.frame(importance(rf_model)) %>%
    tibble::rownames_to_column("feature") %>%
    mutate(feature_label = pretty_feature(feature)) %>%
    arrange(desc(MeanDecreaseGini))
  write_csv(rf_imp, file.path(model_dir, "RF_feature_importance.csv"))

  top_n_rf <- min(12, nrow(rf_imp))
  p_rf_imp <- rf_imp %>%
    slice_head(n = top_n_rf) %>%
    mutate(feature_label = factor(feature_label, levels = rev(feature_label))) %>%
    ggplot(aes(feature_label, MeanDecreaseGini)) +
    geom_col(width = 0.75, fill = "#2C7FB8") +
    coord_flip() +
    theme_classic(base_family = "Times", base_size = 10) +
    labs(
      title = "Random Forest feature importance",
      x = NULL,
      y = "Mean decrease Gini"
    ) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(fig_dir, "RF_feature_importance.pdf"), p_rf_imp, width = 6.8, height = 4.8, device = cairo_pdf)
  ggsave(file.path(fig_dir, "RF_feature_importance.svg"), p_rf_imp, width = 6.8, height = 4.8, device = svg)

  rf_loocv_rows <- list()
  for (i in seq_len(nrow(rf_df))) {
    train <- rf_df[-i, , drop = FALSE]
    test <- rf_df[i, , drop = FALSE]
    set.seed(20260708 + i)
    fit <- randomForest(branch ~ ., data = train, ntree = 1000, importance = FALSE)
    pred <- predict(fit, newdata = test, type = "response")
    probs <- predict(fit, newdata = test, type = "prob")
    rf_loocv_rows[[i]] <- data.frame(
      target_gene_id = feature$target_gene_id[i],
      true_branch = as.character(feature$branch[i]),
      predicted_branch = as.character(pred),
      correct = as.character(pred) == as.character(feature$branch[i]),
      max_probability = max(probs[1, ]),
      stringsAsFactors = FALSE
    )
  }
  rf_loocv <- bind_rows(rf_loocv_rows)
  write_csv(rf_loocv, file.path(model_dir, "RF_LOOCV_predictions.csv"))
  rf_loocv_conf <- rf_loocv %>%
    count(true_branch, predicted_branch, name = "n") %>%
    complete(true_branch = branch_order, predicted_branch = branch_order, fill = list(n = 0))
  write_csv(rf_loocv_conf, file.path(model_dir, "RF_LOOCV_confusion_matrix_long.csv"))
  rf_branch_acc <- rf_loocv %>%
    group_by(true_branch) %>%
    summarise(n = n(), accuracy = mean(correct), .groups = "drop")
  write_csv(rf_branch_acc, file.path(model_dir, "RF_LOOCV_branch_accuracy.csv"))
  rf_loocv_summary <- data.frame(
    metric = c("RF_LOOCV_overall_accuracy", "RF_LOOCV_balanced_accuracy"),
    value = c(mean(rf_loocv$correct), mean(rf_branch_acc$accuracy))
  )
  write_csv(rf_loocv_summary, file.path(model_dir, "RF_LOOCV_summary.csv"))

  p_rf_conf <- ggplot(rf_loocv_conf, aes(predicted_branch, true_branch, fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), family = "Times", size = 3.5) +
    scale_fill_gradient(low = "white", high = "#008B68") +
    theme_classic(base_family = "Times", base_size = 10) +
    labs(
      title = "Random Forest LOOCV confusion matrix",
      x = "Predicted branch",
      y = "True branch",
      fill = "Count"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
  ggsave(file.path(fig_dir, "RF_LOOCV_confusion_matrix.pdf"), p_rf_conf, width = 6.5, height = 5.2, device = cairo_pdf)
  ggsave(file.path(fig_dir, "RF_LOOCV_confusion_matrix.svg"), p_rf_conf, width = 6.5, height = 5.2, device = svg)
}

# Exploratory decision tree and variable importance.
tree_df <- feature %>% select(branch, all_of(binary_cols))
tree <- rpart(branch ~ ., data = tree_df, method = "class", control = rpart.control(cp = 0.02, minsplit = 4))
tree_importance <- data.frame(
  feature = names(tree$variable.importance),
  rpart_importance = as.numeric(tree$variable.importance)
) %>%
  mutate(feature_label = pretty_feature(feature)) %>%
  arrange(desc(rpart_importance))
write_csv(tree_importance, file.path(model_dir, "rpart_exploratory_feature_importance.csv"))

pdf(file.path(fig_dir, "rpart_exploratory_decision_tree.pdf"), width = 9, height = 5.5)
plot(tree, uniform = TRUE, branch = 0.6, margin = 0.08)
text(tree, use.n = TRUE, cex = 0.65)
title("Exploratory decision tree for K08356 neighborhood branch labels")
dev.off()

# Univariate prevalence/importance.
importance_rows <- list()
for (col in binary_cols) {
  x <- feature[[col]]
  p <- tryCatch(fisher.test(table(x, feature$branch), simulate.p.value = TRUE, B = 9999)$p.value, error = function(e) NA_real_)
  prev <- feature %>%
    group_by(branch) %>%
    summarise(prevalence = mean(.data[[col]]), .groups = "drop")
  importance_rows[[col]] <- data.frame(
    feature = col,
    feature_label = pretty_feature(col),
    fisher_p = p,
    prevalence_range = max(prev$prevalence) - min(prev$prevalence),
    max_prevalence_branch = as.character(prev$branch[which.max(prev$prevalence)]),
    stringsAsFactors = FALSE
  )
}
importance <- bind_rows(importance_rows) %>%
  mutate(fisher_fdr = p.adjust(fisher_p, method = "BH")) %>%
  arrange(fisher_fdr, desc(prevalence_range))
write_csv(importance, file.path(model_dir, "univariate_binary_feature_importance.csv"))

p_imp <- importance %>%
  mutate(feature_label = factor(feature_label, levels = rev(feature_label))) %>%
  ggplot(aes(feature_label, -log10(fisher_fdr), fill = max_prevalence_branch)) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_classic(base_family = "Times", base_size = 10) +
  labs(
    title = "Binary neighborhood feature association with K08356 branch",
    x = NULL,
    y = "-log10(FDR)",
    fill = "Highest prevalence branch"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(fig_dir, "univariate_binary_feature_importance.pdf"), p_imp, width = 7.2, height = 4.8, device = cairo_pdf)
ggsave(file.path(fig_dir, "univariate_binary_feature_importance.svg"), p_imp, width = 7.2, height = 4.8, device = svg)

# Prevalence by branch heatmap as a tidy CSV and figure.
prevalence <- read_csv(file.path(model_dir, "gene_neighborhood_feature_prevalence_by_branch.csv"), show_col_types = FALSE)
prevalence_plot <- prevalence %>%
  mutate(
    branch = factor(branch, levels = branch_order),
    function_class = factor(function_class, levels = rev(unique(function_class)))
  ) %>%
  ggplot(aes(branch, function_class, fill = prevalence)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", prevalence)), family = "Times", size = 2.7) +
  scale_fill_gradient(low = "white", high = "#008B68", limits = c(0, 1)) +
  theme_classic(base_family = "Times", base_size = 10) +
  labs(
    title = "Prevalence of neighboring functional classes by K08356 branch",
    x = NULL,
    y = NULL,
    fill = "Prevalence"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
ggsave(file.path(fig_dir, "feature_prevalence_by_branch_heatmap.pdf"), prevalence_plot, width = 7.8, height = 5.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "feature_prevalence_by_branch_heatmap.svg"), prevalence_plot, width = 7.8, height = 5.2, device = svg)

message("Wrote gene-neighborhood ML v1 results to: ", model_dir)
message(sprintf("LOOCV overall accuracy = %.3f; balanced accuracy = %.3f", overall_accuracy, balanced_accuracy))

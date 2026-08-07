suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(ape)
  library(randomForest)
  library(patchwork)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm/机器学习"
model_dir <- file.path(base_dir, "01_gene_neighborhood_model")
fig_dir <- file.path(base_dir, Sys.getenv("K08356_ML_FIG_DIR", "figures_v2"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

branch_order <- c(
  "IdrA-associated",
  "canonical aioA-associated",
  "aioA-like-associated",
  "unknown AioA-like / uncertain DMSOR"
)

branch_colors <- c(
  "IdrA-associated" = "#5B8DB8",
  "canonical aioA-associated" = "#2CA25F",
  "aioA-like-associated" = "#009E73",
  "unknown AioA-like / uncertain DMSOR" = "#CC79A7"
)

feature_order <- c(
  "Rieske/AioB",
  "Cytochrome/peroxidase",
  "Arsenic resistance",
  "Transporter",
  "Sox/sulfur oxidation",
  "Redox/oxidoreductase",
  "Annotated other",
  "Hypothetical/unknown"
)

pretty_branch <- function(x) {
  dplyr::recode(
    as.character(x),
    "IdrA-associated" = "IdrA-associated",
    "canonical aioA-associated" = "canonical aioA",
    "aioA-like-associated" = "aioA-like",
    "unknown AioA-like / uncertain DMSOR" = "unknown/uncertain"
  )
}

pretty_feature <- function(x) {
  x %>%
    gsub("^has_", "", .) %>%
    gsub("_", " ", .) %>%
    gsub("cyt peroxidase", "cyt/peroxidase", .) %>%
    gsub("Rieske AioB", "Rieske/AioB", .) %>%
    gsub("sox sulfur oxidation", "sox/sulfur oxidation", .) %>%
    gsub("redox oxidoreductase", "redox/oxidoreductase", .)
}

feature <- read_csv(file.path(model_dir, "gene_neighborhood_feature_matrix.csv"), show_col_types = FALSE) %>%
  mutate(branch = factor(branch, levels = branch_order))

binary_cols <- grep("^has_", names(feature), value = TRUE)
zero_var_file <- file.path(model_dir, "zero_variance_binary_features_excluded_from_models.csv")
if (file.exists(zero_var_file)) {
  zero_var <- read_csv(zero_var_file, show_col_types = FALSE)$feature
  binary_cols <- setdiff(binary_cols, zero_var)
}
X <- feature %>% select(all_of(binary_cols))
X_mat <- as.matrix(X)
rownames(X_mat) <- feature$target_gene_id

permanova <- read_csv(file.path(model_dir, "PERMANOVA_jaccard_binary_features.csv"), show_col_types = FALSE)
perm_row <- permanova %>% filter(term == "Model") %>% slice(1)
perm_label <- sprintf("PERMANOVA: R² = %.3f, p = %.3f", perm_row$R2, perm_row$`Pr(>F)`)

# RF label permutation test using balanced accuracy.
rf_loocv <- read_csv(file.path(model_dir, "RF_LOOCV_predictions.csv"), show_col_types = FALSE)
obs_branch_acc <- rf_loocv %>%
  group_by(true_branch) %>%
  summarise(acc = mean(correct), .groups = "drop")
obs_bal_acc <- mean(obs_branch_acc$acc)

perm_summary_file <- file.path(model_dir, "RF_label_permutation_test_balanced_accuracy.csv")
perm_dist_file <- file.path(model_dir, "RF_label_permutation_distribution.csv")
if (file.exists(perm_summary_file) && file.exists(perm_dist_file)) {
  perm_tbl_existing <- read_csv(perm_summary_file, show_col_types = FALSE)
  perm_dist_existing <- read_csv(perm_dist_file, show_col_types = FALSE)
  perm_bal_acc <- perm_dist_existing$permutation_balanced_accuracy
  perm_p <- perm_tbl_existing$permutation_p_value[1]
  n_perm <- perm_tbl_existing$n_permutations[1]
} else {
  set.seed(20260708)
  n_perm <- 200
  perm_bal_acc <- numeric(n_perm)
  rf_df <- feature %>% select(branch, all_of(binary_cols))
  for (b in seq_len(n_perm)) {
    perm_branch <- sample(rf_df$branch)
    pred_rows <- vector("list", nrow(rf_df))
    for (i in seq_len(nrow(rf_df))) {
      train <- rf_df[-i, , drop = FALSE]
      train$branch <- perm_branch[-i]
      test <- rf_df[i, , drop = FALSE]
      set.seed(900000 + b * 100 + i)
      fit <- randomForest(branch ~ ., data = train, ntree = 100, importance = FALSE)
      pred <- predict(fit, newdata = test, type = "response")
      pred_rows[[i]] <- data.frame(true_branch = as.character(perm_branch[i]), predicted_branch = as.character(pred))
    }
    pred_df <- bind_rows(pred_rows)
    ba <- pred_df %>%
      group_by(true_branch) %>%
      summarise(acc = mean(true_branch == predicted_branch), .groups = "drop")
    perm_bal_acc[b] <- mean(ba$acc)
  }
  perm_p <- (sum(perm_bal_acc >= obs_bal_acc) + 1) / (n_perm + 1)
  perm_tbl <- data.frame(
    observed_balanced_accuracy = obs_bal_acc,
    permutation_mean_balanced_accuracy = mean(perm_bal_acc),
    permutation_sd_balanced_accuracy = sd(perm_bal_acc),
    permutation_p_value = perm_p,
    n_permutations = n_perm
  )
  write_csv(perm_tbl, perm_summary_file)
  write_csv(data.frame(permutation_balanced_accuracy = perm_bal_acc), perm_dist_file)
}

theme_pub <- function(base_size = 10) {
  theme_classic(base_family = "Times", base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

# A. Feature prevalence heatmap with Fisher FDR stars.
prevalence <- read_csv(file.path(model_dir, "gene_neighborhood_feature_prevalence_by_branch.csv"), show_col_types = FALSE)
tests <- read_csv(file.path(model_dir, "univariate_binary_feature_importance.csv"), show_col_types = FALSE) %>%
  mutate(
    function_class = case_when(
      feature == "has_Rieske_AioB" ~ "Rieske/AioB",
      feature == "has_cyt_peroxidase" ~ "Cytochrome/peroxidase",
      feature == "has_arsenic_resistance" ~ "Arsenic resistance",
      feature == "has_transporter" ~ "Transporter",
      feature == "has_sox_sulfur_oxidation" ~ "Sox/sulfur oxidation",
      feature == "has_redox_oxidoreductase" ~ "Redox/oxidoreductase",
      feature == "has_annotated_other" ~ "Annotated other",
      feature == "has_hypothetical_unknown" ~ "Hypothetical/unknown",
      TRUE ~ NA_character_
    ),
    star = case_when(
      fisher_fdr < 0.001 ~ "***",
      fisher_fdr < 0.01 ~ "**",
      fisher_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  filter(!is.na(function_class)) %>%
  select(function_class, star)

prev_plot_df <- prevalence %>%
  filter(function_class %in% feature_order) %>%
  left_join(tests, by = "function_class") %>%
  mutate(
    branch = factor(branch, levels = branch_order),
    branch_label = factor(pretty_branch(branch), levels = pretty_branch(branch_order)),
    function_class = factor(function_class, levels = rev(feature_order)),
    label = paste0(sprintf("%.0f%%", 100 * prevalence), ifelse(is.na(star), "", star))
  )

p_prev <- ggplot(prev_plot_df, aes(branch_label, function_class, fill = prevalence)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), family = "Times", size = 3.1) +
  scale_fill_gradient(low = "white", high = "#008B68", limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  theme_pub(10) +
  labs(
    title = "A. Gene-neighborhood feature prevalence",
    subtitle = "Asterisks indicate FDR-adjusted Fisher's exact test significance among branches (* FDR < 0.05; ** FDR < 0.01).",
    x = NULL,
    y = NULL,
    fill = "Prevalence"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5)
  )

ggsave(file.path(fig_dir, "A_feature_prevalence_heatmap.pdf"), p_prev, width = 7.8, height = 5.0, device = cairo_pdf)
ggsave(file.path(fig_dir, "A_feature_prevalence_heatmap.svg"), p_prev, width = 7.8, height = 5.0, device = svg)
ggsave(file.path(fig_dir, "A_feature_prevalence_heatmap.png"), p_prev, width = 7.8, height = 5.0, dpi = 450, bg = "white")

# B. Jaccard PCoA.
jaccard <- vegdist(X_mat, method = "jaccard", binary = TRUE)
pcoa <- ape::pcoa(as.dist(jaccard))
pcoa_scores <- as.data.frame(pcoa$vectors[, 1:2]) %>%
  rename(PCoA1 = Axis.1, PCoA2 = Axis.2) %>%
  mutate(
    target_gene_id = feature$target_gene_id,
    branch = feature$branch,
    branch_label = pretty_branch(branch),
    habitat = feature$habitat
  )
write_csv(pcoa_scores, file.path(model_dir, "gene_neighborhood_Jaccard_PCoA_scores.csv"))
pcoa_var <- pcoa$values$Relative_eig[1:2] * 100

p_pcoa <- ggplot(pcoa_scores, aes(PCoA1, PCoA2, color = branch, shape = habitat)) +
  geom_point(size = 3.0, alpha = 0.92) +
  annotate("text", x = -Inf, y = Inf, label = perm_label, hjust = -0.05, vjust = 1.25, family = "Times", size = 3.2) +
  scale_color_manual(values = branch_colors, labels = pretty_branch(branch_order)) +
  theme_pub(10) +
  labs(
    title = "B. Jaccard PCoA of neighborhood features",
    subtitle = "Binary target-centered gene-neighborhood features; target K08356 ORF excluded",
    x = sprintf("PCoA1 (%.1f%%)", pcoa_var[1]),
    y = sprintf("PCoA2 (%.1f%%)", pcoa_var[2]),
    color = "Branch",
    shape = "Habitat"
  ) +
  theme(plot.subtitle = element_text(size = 8.5, hjust = 0.5))
ggsave(file.path(fig_dir, "B_Jaccard_PCoA_PERMANOVA.pdf"), p_pcoa, width = 6.8, height = 5.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "B_Jaccard_PCoA_PERMANOVA.svg"), p_pcoa, width = 6.8, height = 5.2, device = svg)
ggsave(file.path(fig_dir, "B_Jaccard_PCoA_PERMANOVA.png"), p_pcoa, width = 6.8, height = 5.2, dpi = 450, bg = "white")

# C. Random Forest feature importance.
rf_imp <- read_csv(file.path(model_dir, "RF_feature_importance.csv"), show_col_types = FALSE) %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  slice_head(n = 10) %>%
  mutate(
    feature_label = pretty_feature(feature),
    feature_label = factor(feature_label, levels = rev(feature_label))
  )
p_rf <- ggplot(rf_imp, aes(feature_label, MeanDecreaseAccuracy)) +
  geom_col(width = 0.72, fill = "#2C7FB8") +
  coord_flip() +
  theme_pub(10) +
  labs(
    title = "C. Exploratory RF feature ranking",
    x = NULL,
    y = "Mean decrease accuracy"
  ) +
  annotate("text", x = 1, y = max(rf_imp$MeanDecreaseAccuracy, na.rm = TRUE) * 0.55,
           label = "Exploratory Random Forest importance", family = "Times", size = 3)
ggsave(file.path(fig_dir, "C_RF_feature_importance_top10.pdf"), p_rf, width = 6.8, height = 4.8, device = cairo_pdf)
ggsave(file.path(fig_dir, "C_RF_feature_importance_top10.svg"), p_rf, width = 6.8, height = 4.8, device = svg)
ggsave(file.path(fig_dir, "C_RF_feature_importance_top10.png"), p_rf, width = 6.8, height = 4.8, dpi = 450, bg = "white")

# D. Branch-wise RF LOOCV performance.
branch_acc <- read_csv(file.path(model_dir, "RF_LOOCV_branch_accuracy.csv"), show_col_types = FALSE) %>%
  mutate(
    true_branch = factor(true_branch, levels = branch_order),
    branch_label = factor(pretty_branch(true_branch), levels = pretty_branch(branch_order)),
    label = paste0("n=", n, "\n", sprintf("%.0f%%", 100 * accuracy))
  )
p_acc <- ggplot(branch_acc, aes(branch_label, accuracy, fill = true_branch)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_text(aes(label = label), vjust = -0.25, family = "Times", size = 3.1) +
  scale_y_continuous(limits = c(0, 1.12), labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = branch_colors, guide = "none") +
  theme_pub(10) +
  labs(
    title = "D. Branch-wise RF LOOCV performance",
    subtitle = "Canonical aioA-associated performance is unstable because n = 3",
    x = NULL,
    y = "Classification accuracy"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5)
  )
ggsave(file.path(fig_dir, "D_RF_branchwise_LOOCV_accuracy.pdf"), p_acc, width = 6.6, height = 4.8, device = cairo_pdf)
ggsave(file.path(fig_dir, "D_RF_branchwise_LOOCV_accuracy.svg"), p_acc, width = 6.6, height = 4.8, device = svg)
ggsave(file.path(fig_dir, "D_RF_branchwise_LOOCV_accuracy.png"), p_acc, width = 6.6, height = 4.8, dpi = 450, bg = "white")

# E. Permutation distribution.
perm_dist <- data.frame(permutation_balanced_accuracy = perm_bal_acc)
p_perm <- ggplot(perm_dist, aes(permutation_balanced_accuracy)) +
  geom_histogram(bins = 30, fill = "grey80", color = "white") +
  geom_vline(xintercept = obs_bal_acc, color = "#D7301F", linewidth = 0.8) +
  annotate(
    "text",
    x = obs_bal_acc,
    y = Inf,
    label = sprintf("Observed balanced accuracy = %.3f\nPermutation p = %.3f", obs_bal_acc, perm_p),
    hjust = -0.05,
    vjust = 1.2,
    family = "Times",
    size = 3.1
  ) +
  theme_pub(10) +
  labs(
    title = "RF label-permutation test",
    x = "Balanced accuracy under permuted branch labels",
    y = "Frequency"
  )
ggsave(file.path(fig_dir, "E_RF_label_permutation_test.pdf"), p_perm, width = 6.4, height = 4.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "E_RF_label_permutation_test.svg"), p_perm, width = 6.4, height = 4.5, device = svg)
ggsave(file.path(fig_dir, "E_RF_label_permutation_test.png"), p_perm, width = 6.4, height = 4.5, dpi = 450, bg = "white")

# Combined 4-panel figure.
combined <- (p_prev | p_pcoa) / (p_rf | p_acc) +
  plot_annotation(
    title = "Gene-neighborhood signatures of K08356-related DMSOR lineages",
    caption = sprintf(
      "Asterisks indicate FDR-adjusted Fisher's exact test significance among branches (* FDR < 0.05; ** FDR < 0.01). B: Jaccard PCoA from binary features. C: exploratory RF ranking. Target K08356 ORF excluded; RF permutation BA = %.3f, p = %.3f.",
      obs_bal_acc, perm_p
    ),
    theme = theme(
      plot.title = element_text(family = "Times", face = "bold", hjust = 0.5, size = 16),
      plot.caption = element_text(family = "Times", hjust = 0.5, size = 9)
    )
  )
ggsave(file.path(fig_dir, "Figure_gene_neighborhood_ML_four_panel.pdf"), combined, width = 14, height = 10.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "Figure_gene_neighborhood_ML_four_panel.svg"), combined, width = 14, height = 10.2, device = svg)
ggsave(file.path(fig_dir, "Figure_gene_neighborhood_ML_four_panel.png"), combined, width = 14, height = 10.2, dpi = 450, bg = "white")

main_abc <- (p_prev | p_pcoa) / p_rf +
  plot_layout(heights = c(1, 0.82)) +
  plot_annotation(
    title = "Gene-neighborhood signatures of K08356-related DMSOR lineages",
    caption = sprintf(
      "Asterisks indicate FDR-adjusted Fisher's exact test significance among branches (* FDR < 0.05; ** FDR < 0.01). B: Jaccard PCoA from binary features. C: exploratory RF ranking. Target K08356 ORF excluded; RF permutation BA = %.3f, p = %.3f. Branch-wise LOOCV is provided separately.",
      obs_bal_acc, perm_p
    ),
    theme = theme(
      plot.title = element_text(family = "Times", face = "bold", hjust = 0.5, size = 16),
      plot.caption = element_text(family = "Times", hjust = 0.5, size = 9)
    )
  )
ggsave(file.path(fig_dir, "Figure_gene_neighborhood_ML_main_ABC.pdf"), main_abc, width = 14, height = 10.2, device = cairo_pdf)
ggsave(file.path(fig_dir, "Figure_gene_neighborhood_ML_main_ABC.svg"), main_abc, width = 14, height = 10.2, device = svg)
ggsave(file.path(fig_dir, "Figure_gene_neighborhood_ML_main_ABC.png"), main_abc, width = 14, height = 10.2, dpi = 450, bg = "white")

notes <- c(
  "# revised gene-neighborhood ML figure notes",
  "",
  "Generated revised gene-neighborhood ML figures:",
  "- A_feature_prevalence_heatmap.pdf/svg",
  "- B_Jaccard_PCoA_PERMANOVA.pdf/svg",
  "- C_RF_feature_importance_top10.pdf/svg",
  "- D_RF_branchwise_LOOCV_accuracy.pdf/svg",
  "- E_RF_label_permutation_test.pdf/svg",
  "- Figure_gene_neighborhood_ML_four_panel.pdf/svg",
  "- Figure_gene_neighborhood_ML_main_ABC.pdf/svg",
  "",
  sprintf("PERMANOVA: R2 = %.3f, F = %.3f, p = %.3f.", perm_row$R2, perm_row$F, perm_row$`Pr(>F)`),
  sprintf("RF label permutation test: observed balanced accuracy = %.3f; random mean = %.3f; p = %.3f; permutations = %d.",
          obs_bal_acc, mean(perm_bal_acc), perm_p, n_perm),
  "",
  "Interpretation: use Random Forest feature importance as exploratory support, not causal proof."
)
writeLines(notes, file.path(fig_dir, "README_figures_v2.md"))

message("Wrote revised figures to: ", fig_dir)
message(sprintf("RF permutation test: observed BA = %.3f, random mean = %.3f, p = %.3f", obs_bal_acc, mean(perm_bal_acc), perm_p))

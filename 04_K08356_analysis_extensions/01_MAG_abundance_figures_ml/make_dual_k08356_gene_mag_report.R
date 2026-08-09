suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

base_out <- "/Users/catherine/Downloads/MAG-bin-tpm/K08356_gene_and_MAG_dual_evidence_20260708"
gene_out <- file.path(base_out, "01_gene_level_K08356_TPM")
mag_out <- file.path(base_out, "02_MAG_host_TPM")
fig_out <- file.path(base_out, "figures")
dir.create(gene_out, recursive = TRUE, showWarnings = FALSE)
dir.create(mag_out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

habitat_levels <- c("IS", "AS", "ES", "NS")
hab_cols <- c(IS = "#2ca25f", AS = "#8e63c7", ES = "#e39d25", NS = "#3b78c8")
branch_cols <- c(
  "IdrA-associated" = "#5B8DB8",
  "canonical aioA-associated" = "#D55E00",
  "aioA-like-associated" = "#009E73",
  "unknown AioA-like / uncertain DMSOR" = "#CC79A7"
)

sample_to_habitat <- function(sample) {
  case_when(
    grepl("^SY", sample) ~ "IS",
    grepl("^SQ_", sample) ~ "AS",
    grepl("^(S1_|S2_|S4_)", sample) ~ "AS",
    grepl("^(C1_|C2_|C3_|S13_|S14_|S15_|ES_)", sample) ~ "ES",
    grepl("^(S3_|R2111_|NS_)", sample) ~ "NS",
    TRUE ~ NA_character_
  )
}

sample_group <- function(sample) {
  case_when(
    grepl("^SY365", sample) ~ "SY365",
    grepl("^SY366YB", sample) ~ "SY366YB",
    grepl("^SY366YW", sample) ~ "SY366YW",
    grepl("^SY368", sample) ~ "SY368",
    grepl("^SY456", sample) ~ "SY456",
    grepl("^SY457", sample) ~ "SY457",
    grepl("^SY459", sample) ~ "SY459",
    grepl("^SQ_58", sample) ~ "SQ_58",
    grepl("^SQ_81", sample) ~ "SQ_81",
    grepl("^R2111_N500", sample) ~ "R2111_N500",
    grepl("^R2111_S500", sample) ~ "R2111_S500",
    grepl("^R2111_N300", sample) ~ "R2111_N300",
    grepl("^R2111_S300", sample) ~ "R2111_S300",
    grepl("^ES_2", sample) ~ "ES_2",
    grepl("^NS_", sample) ~ "NS",
    TRUE ~ sub("_.*$", "", sample)
  )
}

depth_start <- function(sample) {
  x <- sub("^.*_([0-9]+)-[0-9]+$", "\\1", sample)
  suppressWarnings(as.numeric(ifelse(x == sample, 0, x)))
}

make_sample_order <- function(samples) {
  group_order <- c(
    "SY365", "SY366YB", "SY366YW", "SY368", "SY456", "SY457", "SY459",
    "S1", "S2", "S4", "SQ_58", "SQ_81",
    "C1", "C2", "C3", "ES_2", "S13", "S14", "S15",
    "S3", "NS", "R2111_N500", "R2111_S500", "R2111_N300", "R2111_S300"
  )
  tibble(sample = unique(samples)) %>%
    mutate(
      habitat = factor(sample_to_habitat(sample), levels = habitat_levels),
      group = sample_group(sample),
      group_rank = match(group, group_order),
      group_rank = ifelse(is.na(group_rank), 999, group_rank),
      depth = depth_start(sample)
    ) %>%
    arrange(habitat, group_rank, depth, sample) %>%
    pull(sample)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 3))))
}

kw_table <- function(df, feature_col = "Feature", value_col = "logTPM", group_col = "Habitat") {
  df %>%
    group_by(.data[[feature_col]]) %>%
    group_modify(~ {
      if (nrow(.x) < 3 || length(unique(.x[[group_col]])) < 2) {
        return(tibble(statistic = NA_real_, df = NA_real_, p_value = NA_real_, n = nrow(.x)))
      }
      kt <- kruskal.test(.x[[value_col]] ~ .x[[group_col]])
      tibble(statistic = unname(kt$statistic), df = unname(kt$parameter), p_value = kt$p.value, n = nrow(.x))
    }) %>%
    ungroup()
}

pairwise_table <- function(df, feature_col = "Feature", value_col = "logTPM", group_col = "Habitat") {
  out <- list()
  for (ft in unique(df[[feature_col]])) {
    dat <- df %>% filter(.data[[feature_col]] == ft)
    if (nrow(dat) >= 3 && length(unique(dat[[group_col]])) >= 2) {
      pw <- pairwise.wilcox.test(dat[[value_col]], dat[[group_col]], p.adjust.method = "BH", exact = FALSE)
      mat <- pw$p.value
      if (!is.null(mat)) {
        for (r in rownames(mat)) {
          for (c in colnames(mat)) {
            p <- mat[r, c]
            if (!is.na(p)) {
              out[[length(out) + 1]] <- tibble(Feature = ft, group1 = c, group2 = r, p_adj_BH = p)
            }
          }
        }
      }
    }
  }
  bind_rows(out)
}

theme_pub <- theme_classic(base_family = "Times", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.text = element_text(color = "grey15"),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    plot.margin = margin(10, 14, 10, 10)
  )

## 1. Gene/contig-level K08356 TPM, matching previous contig workflow
gene_file <- "/Users/catherine/Downloads/王镨蒂师兄交作业系列 ddl/师兄给的 contig 功能基因相对丰度数据/all_koid_tpm.xlsx"
gene_raw <- read_excel(gene_file, sheet = "all_koid_tpm")
gene_total <- gene_raw %>%
  filter(as.character(KO_ID) == "K08356") %>%
  pivot_longer(-KO_ID, names_to = "Sample", values_to = "TPM") %>%
  mutate(
    TPM = as.numeric(TPM),
    TPM = ifelse(is.na(TPM), 0, TPM),
    logTPM = log10(TPM + 1),
    Habitat = factor(sample_to_habitat(Sample), levels = habitat_levels),
    Sample_ordered = factor(Sample, levels = make_sample_order(Sample)),
    Feature = "K08356 gene/contig-level TPM"
  ) %>%
  filter(!is.na(Habitat))

write_csv(gene_total, file.path(gene_out, "K08356_gene_total_sample_long.csv"))
gene_summary <- gene_total %>%
  group_by(Habitat) %>%
  summarise(
    n_samples = n(),
    n_detected_TPM_gt_0 = sum(TPM > 0),
    detection_rate = n_detected_TPM_gt_0 / n_samples,
    mean_TPM = mean(TPM),
    median_TPM = median(TPM),
    max_TPM = max(TPM),
    mean_log10_TPM_plus_1 = mean(logTPM),
    median_log10_TPM_plus_1 = median(logTPM),
    .groups = "drop"
  )
write_csv(gene_summary, file.path(gene_out, "K08356_gene_total_habitat_summary.csv"))
gene_kw <- kw_table(gene_total)
write_csv(gene_kw, file.path(gene_out, "K08356_gene_total_Kruskal.csv"))
gene_pw <- pairwise_table(gene_total)
write_csv(gene_pw, file.path(gene_out, "K08356_gene_total_pairwise_Wilcoxon_BH.csv"))

p_gene_box <- ggplot(gene_total, aes(Habitat, logTPM, fill = Habitat)) +
  geom_boxplot(width = 0.58, alpha = 0.72, color = "grey25", outlier.shape = NA, linewidth = 0.45) +
  geom_jitter(aes(color = Habitat), width = 0.12, size = 2.2, alpha = 0.9, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.0, fill = "white", color = "black") +
  scale_fill_manual(values = hab_cols, drop = FALSE) +
  scale_color_manual(values = hab_cols, drop = FALSE) +
  labs(
    title = "Contig/KO-level total K08356 TPM",
    subtitle = paste0("Kruskal-Wallis on log10(TPM+1): ", fmt_p(gene_kw$p_value[1])),
    x = "Habitat",
    y = "K08356 gene abundance, log10(TPM + 1)"
  ) +
  theme_pub +
  theme(legend.position = "none")
ggsave(file.path(fig_out, "01_gene_total_K08356_TPM_boxplot.png"), p_gene_box, width = 6.4, height = 6.0, dpi = 450)
ggsave(file.path(fig_out, "01_gene_total_K08356_TPM_boxplot.pdf"), p_gene_box, width = 6.4, height = 6.0)

p_gene_bar <- ggplot(gene_total, aes(Sample_ordered, TPM, fill = Habitat)) +
  geom_col(width = 0.82, color = "grey35", linewidth = 0.12) +
  scale_fill_manual(values = hab_cols, drop = FALSE) +
  labs(
    title = "Contig/KO-level K08356 TPM by sample order",
    subtitle = "One K08356 gene/KO TPM value per sample",
    x = "Sample ordered by habitat/station",
    y = "K08356 gene TPM"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
ggsave(file.path(fig_out, "02_gene_total_K08356_TPM_sample_order_bar.png"), p_gene_bar, width = 12, height = 6.5, dpi = 450)
ggsave(file.path(fig_out, "02_gene_total_K08356_TPM_sample_order_bar.pdf"), p_gene_bar, width = 12, height = 6.5)

## 2. MAG host TPM for K08356-bearing derepMAGs
mag_matrix_file <- "/Users/catherine/Downloads/MAG-bin-tpm/merged_tpm_matrix.txt"
branch_map_file <- "/Users/catherine/Downloads/MAG-bin-tpm/derepMAG_abundance_stats/k08356_branch_feature_map.tsv"
mag_matrix <- read.delim(mag_matrix_file, check.names = FALSE)
colnames(mag_matrix)[1] <- "MAG"
branch_map <- read.delim(branch_map_file, check.names = FALSE) %>%
  select(MAG, Feature, Original_MAG_ID) %>%
  distinct() %>%
  group_by(MAG) %>%
  mutate(
    K08356_copy_rows_in_MAG = n(),
    copy_fraction = 1 / K08356_copy_rows_in_MAG
  ) %>%
  ungroup()

branch_denominators <- branch_map %>%
  group_by(Feature) %>%
  summarise(
    n_K08356_copy_rows = n(),
    n_unique_MAGs = n_distinct(MAG),
    sum_copy_fraction = sum(copy_fraction),
    .groups = "drop"
  )
write_csv(branch_denominators, file.path(mag_out, "K08356_branch_MAG_denominators.csv"))

mag_long <- mag_matrix %>%
  inner_join(branch_map, by = "MAG") %>%
  pivot_longer(-c(MAG, Feature, Original_MAG_ID, K08356_copy_rows_in_MAG, copy_fraction), names_to = "Sample", values_to = "raw_MAG_TPM") %>%
  mutate(
    raw_MAG_TPM = as.numeric(raw_MAG_TPM),
    raw_MAG_TPM = ifelse(is.na(raw_MAG_TPM), 0, raw_MAG_TPM),
    TPM = raw_MAG_TPM * copy_fraction
  )

mag_branch_sample <- mag_long %>%
  group_by(Feature, Sample) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    logTPM = log10(TPM + 1),
    Habitat = factor(sample_to_habitat(Sample), levels = habitat_levels),
    Sample_ordered = factor(Sample, levels = make_sample_order(Sample))
  ) %>%
  filter(!is.na(Habitat))

unique_mags <- branch_map %>% distinct(MAG)
mag_unique_total <- mag_matrix %>%
  inner_join(unique_mags, by = "MAG") %>%
  pivot_longer(-MAG, names_to = "Sample", values_to = "TPM") %>%
  group_by(Sample) %>%
  summarise(TPM = sum(as.numeric(TPM), na.rm = TRUE), .groups = "drop") %>%
  mutate(
    logTPM = log10(TPM + 1),
    Habitat = factor(sample_to_habitat(Sample), levels = habitat_levels),
    Sample_ordered = factor(Sample, levels = make_sample_order(Sample)),
    Feature = "K08356-bearing host MAG total TPM"
  ) %>%
  filter(!is.na(Habitat))

mag_branch_mean_per_copy <- mag_branch_sample %>%
  left_join(branch_denominators, by = "Feature") %>%
  mutate(
    TPM = TPM / n_K08356_copy_rows,
    logTPM = log10(TPM + 1),
    Metric = "copy-number-adjusted mean TPM per branch-assigned K08356 copy row"
  )

mag_branch_prop <- mag_branch_sample %>%
  select(Feature, Sample, Branch_TPM = TPM, Habitat, Sample_ordered) %>%
  left_join(mag_unique_total %>% select(Sample, Total_K08356_host_MAG_TPM = TPM), by = "Sample") %>%
  mutate(
    TPM = ifelse(Total_K08356_host_MAG_TPM > 0, Branch_TPM / Total_K08356_host_MAG_TPM, 0),
    logTPM = TPM,
    Metric = "branch proportion within unique K08356-bearing host MAG total"
  )

write_csv(mag_branch_sample, file.path(mag_out, "K08356_branch_host_MAG_TPM_by_sample.csv"))
write_csv(mag_unique_total, file.path(mag_out, "K08356_unique_host_MAG_total_TPM_by_sample.csv"))
write_csv(mag_branch_mean_per_copy, file.path(mag_out, "K08356_branch_host_MAG_mean_per_copy_by_sample.csv"))
write_csv(mag_branch_prop, file.path(mag_out, "K08356_branch_host_MAG_proportion_by_sample.csv"))

mag_summary <- mag_branch_sample %>%
  group_by(Feature, Habitat) %>%
  summarise(
    n_samples = n(),
    n_detected_TPM_gt_0 = sum(TPM > 0),
    detection_rate = n_detected_TPM_gt_0 / n_samples,
    mean_TPM = mean(TPM),
    median_TPM = median(TPM),
    max_TPM = max(TPM),
    mean_log10_TPM_plus_1 = mean(logTPM),
    median_log10_TPM_plus_1 = median(logTPM),
    .groups = "drop"
  )
write_csv(mag_summary, file.path(mag_out, "K08356_branch_host_MAG_habitat_summary.csv"))

mag_mean_summary <- mag_branch_mean_per_copy %>%
  group_by(Feature, Habitat) %>%
  summarise(
    n_samples = n(),
    n_detected_TPM_gt_0 = sum(TPM > 0),
    detection_rate = n_detected_TPM_gt_0 / n_samples,
    mean_TPM = mean(TPM),
    median_TPM = median(TPM),
    max_TPM = max(TPM),
    mean_log10_TPM_plus_1 = mean(logTPM),
    median_log10_TPM_plus_1 = median(logTPM),
    .groups = "drop"
  )
write_csv(mag_mean_summary, file.path(mag_out, "K08356_branch_host_MAG_mean_per_copy_habitat_summary.csv"))

mag_prop_summary <- mag_branch_prop %>%
  group_by(Feature, Habitat) %>%
  summarise(
    n_samples = n(),
    n_detected_prop_gt_0 = sum(TPM > 0),
    detection_rate = n_detected_prop_gt_0 / n_samples,
    mean_proportion = mean(TPM),
    median_proportion = median(TPM),
    max_proportion = max(TPM),
    .groups = "drop"
  )
write_csv(mag_prop_summary, file.path(mag_out, "K08356_branch_host_MAG_proportion_habitat_summary.csv"))

mag_kw <- kw_table(mag_branch_sample)
write_csv(mag_kw, file.path(mag_out, "K08356_branch_host_MAG_Kruskal.csv"))
mag_pw <- pairwise_table(mag_branch_sample)
write_csv(mag_pw, file.path(mag_out, "K08356_branch_host_MAG_pairwise_Wilcoxon_BH.csv"))
mag_mean_kw <- kw_table(mag_branch_mean_per_copy)
write_csv(mag_mean_kw, file.path(mag_out, "K08356_branch_host_MAG_mean_per_copy_Kruskal.csv"))
mag_prop_kw <- kw_table(mag_branch_prop)
write_csv(mag_prop_kw, file.path(mag_out, "K08356_branch_host_MAG_proportion_Kruskal.csv"))
mag_prop_pw <- pairwise_table(mag_branch_prop)
write_csv(mag_prop_pw, file.path(mag_out, "K08356_branch_host_MAG_proportion_pairwise_Wilcoxon_BH.csv"))
mag_total_kw <- kw_table(mag_unique_total)
write_csv(mag_total_kw, file.path(mag_out, "K08356_unique_host_MAG_total_Kruskal.csv"))
mag_total_pw <- pairwise_table(mag_unique_total)
write_csv(mag_total_pw, file.path(mag_out, "K08356_unique_host_MAG_total_pairwise_Wilcoxon_BH.csv"))

kw_labels <- mag_kw %>%
  mutate(label = paste0("KW ", fmt_p(p_value)))

p_mag_branch_box <- ggplot(mag_branch_sample, aes(Habitat, logTPM, fill = Habitat)) +
  geom_boxplot(width = 0.58, alpha = 0.72, color = "grey25", outlier.shape = NA, linewidth = 0.45) +
  geom_jitter(aes(color = Habitat), width = 0.12, size = 1.7, alpha = 0.88, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white", color = "black") +
  geom_text(data = kw_labels, aes(x = 2.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.35, size = 3.1) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = hab_cols, drop = FALSE) +
  scale_color_manual(values = hab_cols, drop = FALSE) +
  labs(
    title = "K08356-bearing host MAG TPM by branch",
    subtitle = "MAG TPM matrix is column-normalized; this measures host MAG abundance, not K08356 gene TPM",
    x = "Habitat",
    y = "Host MAG abundance, log10(TPM + 1)"
  ) +
  theme_pub
ggsave(file.path(fig_out, "03_MAG_host_TPM_by_branch_boxplots.png"), p_mag_branch_box, width = 8, height = 8, dpi = 450)
ggsave(file.path(fig_out, "03_MAG_host_TPM_by_branch_boxplots.pdf"), p_mag_branch_box, width = 8, height = 8)

p_mag_total_box <- ggplot(mag_unique_total, aes(Habitat, logTPM, fill = Habitat)) +
  geom_boxplot(width = 0.58, alpha = 0.72, color = "grey25", outlier.shape = NA, linewidth = 0.45) +
  geom_jitter(aes(color = Habitat), width = 0.12, size = 2.1, alpha = 0.88, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.0, fill = "white", color = "black") +
  scale_fill_manual(values = hab_cols, drop = FALSE) +
  scale_color_manual(values = hab_cols, drop = FALSE) +
  labs(
    title = "Total K08356-bearing host MAG TPM",
    subtitle = paste0("Kruskal-Wallis on log10(TPM+1): ", fmt_p(mag_total_kw$p_value[1])),
    x = "Habitat",
    y = "Host MAG abundance, log10(TPM + 1)"
  ) +
  theme_pub +
  theme(legend.position = "none")
ggsave(file.path(fig_out, "04_MAG_host_total_TPM_boxplot.png"), p_mag_total_box, width = 6.4, height = 6.0, dpi = 450)
ggsave(file.path(fig_out, "04_MAG_host_total_TPM_boxplot.pdf"), p_mag_total_box, width = 6.4, height = 6.0)

p_mag_stack <- mag_branch_sample %>%
  mutate(Feature = factor(Feature, levels = names(branch_cols))) %>%
  ggplot(aes(Sample_ordered, TPM, fill = Feature)) +
  geom_col(width = 0.82, color = "grey35", linewidth = 0.08) +
  scale_fill_manual(values = branch_cols, drop = FALSE) +
  labs(
    title = "K08356-bearing host MAG TPM by sample order",
    subtitle = "Stacked by K08356-associated MAG branch; values are host MAG TPM",
    x = "Sample ordered by habitat/station",
    y = "K08356-bearing host MAG TPM",
    fill = "MAG branch"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
ggsave(file.path(fig_out, "05_MAG_host_branch_TPM_sample_order_stacked_bar.png"), p_mag_stack, width = 12, height = 6.8, dpi = 450)
ggsave(file.path(fig_out, "05_MAG_host_branch_TPM_sample_order_stacked_bar.pdf"), p_mag_stack, width = 12, height = 6.8)

method_note <- tibble(
  analysis = c("Gene/contig-level K08356 TPM", "MAG host TPM"),
  input_file = c(gene_file, mag_matrix_file),
  abundance_meaning = c(
    "Functional gene/KO abundance: one K08356 TPM value per sample from all_koid_tpm.xlsx",
    "Host genome abundance: sum of TPM values for derepMAGs carrying each K08356-associated branch"
  ),
  normalization = c(
    "TPM normalized upstream in contig/KO workflow; no habitat-sum normalization before tests",
    "MAG matrix columns are already normalized to ~1,000,000 TPM; no second normalization before tests; multi-clade/copy MAGs are split equally among their K08356 copy rows"
  ),
  recommended_use = c(
    "Use for claims that K08356 functional gene itself is enriched in IS or other habitats",
    "Use for claims about abundance/spatial distribution of K08356-bearing host MAGs; total, proportion, and mean-per-copy metrics are all reported"
  )
)
write_csv(method_note, file.path(base_out, "README_method_and_interpretation.csv"))

cat(base_out, "\n")

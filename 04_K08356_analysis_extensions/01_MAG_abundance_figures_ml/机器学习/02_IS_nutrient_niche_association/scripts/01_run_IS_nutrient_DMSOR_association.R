suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(scales)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "机器学习", "02_IS_nutrient_niche_association")
input_dir <- file.path(out_dir, "inputs")
result_dir <- file.path(out_dir, "results")
fig_dir <- file.path(out_dir, "figures")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

nutrient_file <- "/Users/catherine/Downloads/营养盐热图绘制/varechem.csv"
mag_file <- file.path(base_dir, "K08356_gene_and_MAG_dual_evidence_20260708", "02_MAG_host_TPM", "K08356_branch_host_MAG_TPM_by_sample.csv")
gene_file <- file.path(base_dir, "K08356_gene_and_MAG_dual_evidence_20260708", "01_gene_level_K08356_TPM", "K08356_gene_total_sample_long.csv")

theme_pub <- function(base_size = 10) {
  theme_classic(base_family = "Times", base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

label_env <- c(
  NH4_mgkg = "NH4+",
  NO3_mgkg = "NO3-",
  DIN = "DIN",
  NH4_NO3_ratio = "NH4+/NO3-",
  TOP_mgkg = "TOP",
  TOC_mgkg = "TOC",
  TOC_TOP_ratio = "TOC/TOP"
)

label_k <- c(
  log_total_K08356_gene_TPM = "Total K08356 gene TPM",
  log_total_K08356_bearing_MAG_TPM = "Total K08356-bearing MAG TPM",
  log_IdrA_MAG_TPM = "IdrA MAG TPM",
  log_canonical_aioA_MAG_TPM = "canonical aioA MAG TPM",
  log_aioA_like_MAG_TPM = "aioA-like MAG TPM",
  log_unknown_DMSOR_MAG_TPM = "unknown/uncertain DMSOR MAG TPM"
)

clean_feature <- function(x) {
  dplyr::recode(
    x,
    "IdrA-associated" = "IdrA_MAG_TPM",
    "canonical aioA-associated" = "canonical_aioA_MAG_TPM",
    "aioA-like-associated" = "aioA_like_MAG_TPM",
    "unknown AioA-like / uncertain DMSOR" = "unknown_DMSOR_MAG_TPM",
    .default = gsub("[^A-Za-z0-9]+", "_", x)
  )
}

extract_layer <- function(sample) {
  sub("^.*?(\\d+-\\d+)$", "\\1", sample)
}

extract_core <- function(sample) {
  sub("-(\\d+-\\d+)$", "", sample)
}

nutrient_raw <- read_csv(nutrient_file, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
nutrient <- nutrient_raw %>%
  transmute(
    Sample = ID,
    Depth_m = `depth（m)`,
    NH4_mgkg = `NH4+(mg/kg)`,
    NO3_mgkg = `NO3- (mg/kg)`,
    TOP_mgkg = `TOP (mg/kg)`,
    TOC_mgkg = `TOC (mg/kg)`
  ) %>%
  mutate(
    Core = extract_core(Sample),
    Layer = factor(extract_layer(Sample), levels = c("0-4", "4-8", "8-12")),
    Habitat = "IS",
    DIN = NH4_mgkg + NO3_mgkg,
    NH4_NO3_ratio = ifelse(NO3_mgkg > 0, NH4_mgkg / NO3_mgkg, NA_real_),
    TOC_TOP_ratio = ifelse(TOP_mgkg > 0, TOC_mgkg / TOP_mgkg, NA_real_)
  ) %>%
  select(Sample, Habitat, Core, Layer, Depth_m, NH4_mgkg, NO3_mgkg, DIN, NH4_NO3_ratio, TOP_mgkg, TOC_mgkg, TOC_TOP_ratio)

write_csv(nutrient, file.path(input_dir, "IS_nutrient_cleaned.csv"))

mag_long <- read_csv(mag_file, show_col_types = FALSE) %>%
  filter(Habitat == "IS") %>%
  mutate(metric = clean_feature(Feature))

mag_wide <- mag_long %>%
  select(Sample, metric, TPM) %>%
  group_by(Sample, metric) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = metric, values_from = TPM, values_fill = 0)

expected_mag_cols <- c("IdrA_MAG_TPM", "canonical_aioA_MAG_TPM", "aioA_like_MAG_TPM", "unknown_DMSOR_MAG_TPM")
for (cc in expected_mag_cols) {
  if (!cc %in% names(mag_wide)) mag_wide[[cc]] <- 0
}
mag_wide <- mag_wide %>%
  mutate(total_K08356_bearing_MAG_TPM = IdrA_MAG_TPM + canonical_aioA_MAG_TPM + aioA_like_MAG_TPM + unknown_DMSOR_MAG_TPM) %>%
  select(Sample, total_K08356_bearing_MAG_TPM, all_of(expected_mag_cols))

gene_total <- read_csv(gene_file, show_col_types = FALSE) %>%
  filter(Habitat == "IS") %>%
  group_by(Sample) %>%
  summarise(total_K08356_gene_TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

merged <- nutrient %>%
  left_join(gene_total, by = "Sample") %>%
  left_join(mag_wide, by = "Sample") %>%
  mutate(across(c(total_K08356_gene_TPM, total_K08356_bearing_MAG_TPM, all_of(expected_mag_cols)), ~replace_na(.x, 0))) %>%
  mutate(
    log_total_K08356_gene_TPM = log10(total_K08356_gene_TPM + 1),
    log_total_K08356_bearing_MAG_TPM = log10(total_K08356_bearing_MAG_TPM + 1),
    log_IdrA_MAG_TPM = log10(IdrA_MAG_TPM + 1),
    log_canonical_aioA_MAG_TPM = log10(canonical_aioA_MAG_TPM + 1),
    log_aioA_like_MAG_TPM = log10(aioA_like_MAG_TPM + 1),
    log_unknown_DMSOR_MAG_TPM = log10(unknown_DMSOR_MAG_TPM + 1)
  )

match_check <- data.frame(
  source = c("nutrient_samples", "matched_gene_TPM_samples", "matched_MAG_TPM_samples"),
  n = c(nrow(nutrient), sum(nutrient$Sample %in% gene_total$Sample), sum(nutrient$Sample %in% mag_wide$Sample))
)
write_csv(match_check, file.path(result_dir, "sample_matching_check.csv"))
write_csv(merged, file.path(result_dir, "IS_K08356_nutrient_merged_matrix.csv"))

env_vars <- names(label_env)
k_vars <- names(label_k)

rho <- matrix(NA_real_, nrow = length(env_vars), ncol = length(k_vars), dimnames = list(env_vars, k_vars))
pval <- rho
for (e in env_vars) {
  for (k in k_vars) {
    ct <- suppressWarnings(cor.test(merged[[e]], merged[[k]], method = "spearman", exact = FALSE))
    rho[e, k] <- unname(ct$estimate)
    pval[e, k] <- ct$p.value
  }
}
fdr_vals <- p.adjust(as.vector(pval), method = "BH")
fdr <- matrix(fdr_vals, nrow = nrow(pval), ncol = ncol(pval), dimnames = dimnames(pval))

write_csv(as.data.frame(rho) %>% tibble::rownames_to_column("Environmental_variable"), file.path(result_dir, "Spearman_rho_matrix.csv"))
write_csv(as.data.frame(pval) %>% tibble::rownames_to_column("Environmental_variable"), file.path(result_dir, "Spearman_p_matrix.csv"))
write_csv(as.data.frame(fdr) %>% tibble::rownames_to_column("Environmental_variable"), file.path(result_dir, "Spearman_FDR_matrix.csv"))

cor_long <- expand.grid(Environmental_variable = env_vars, K08356_metric = k_vars, stringsAsFactors = FALSE) %>%
  mutate(
    rho = mapply(function(e, k) rho[e, k], Environmental_variable, K08356_metric),
    p_value = mapply(function(e, k) pval[e, k], Environmental_variable, K08356_metric),
    FDR = mapply(function(e, k) fdr[e, k], Environmental_variable, K08356_metric),
    abs_rho = abs(rho),
    env_label = label_env[Environmental_variable],
    metric_label = label_k[K08356_metric],
    star = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      p_value < 0.05 ~ "†",
      TRUE ~ ""
    )
  ) %>%
  arrange(FDR, desc(abs_rho))
write_csv(cor_long, file.path(result_dir, "Spearman_correlation_long.csv"))

priority_metrics <- c("log_aioA_like_MAG_TPM", "log_unknown_DMSOR_MAG_TPM", "log_total_K08356_gene_TPM", "log_total_K08356_bearing_MAG_TPM")
selected_pairs <- cor_long %>%
  mutate(priority = ifelse(K08356_metric %in% priority_metrics, 0, 1)) %>%
  arrange(priority, FDR, p_value, desc(abs_rho)) %>%
  group_by(K08356_metric) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(priority, FDR, p_value, desc(abs_rho)) %>%
  slice_head(n = 4)
if (nrow(selected_pairs) < 4) {
  selected_pairs <- bind_rows(
    selected_pairs,
    anti_join(cor_long, selected_pairs, by = c("Environmental_variable", "K08356_metric")) %>%
      arrange(FDR, p_value, desc(abs_rho)) %>%
      slice_head(n = 4 - nrow(selected_pairs))
  )
}
write_csv(selected_pairs, file.path(result_dir, "selected_scatterplot_pairs.csv"))

# Nutrient boxplots by layer.
nut_long <- nutrient %>%
  select(Sample, Core, Layer, NH4_mgkg, NO3_mgkg, DIN, NH4_NO3_ratio, TOP_mgkg, TOC_mgkg, TOC_TOP_ratio) %>%
  pivot_longer(cols = all_of(env_vars), names_to = "variable", values_to = "value") %>%
  mutate(variable_label = factor(label_env[variable], levels = label_env[env_vars]))

p_box <- ggplot(nut_long, aes(Layer, value, fill = Layer)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.75, linewidth = 0.25) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.85) +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("0-4" = "#74ADD1", "4-8" = "#ABDDA4", "8-12" = "#FDAE61")) +
  theme_pub(10) +
  labs(
    title = "Nutrient and organic matter gradients within IS samples",
    x = "Sediment layer (cm)",
    y = NULL,
    fill = "Layer"
  ) +
  theme(strip.background = element_rect(fill = "grey95", color = "grey70"), strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "nutrient_layer_boxplots.pdf"), p_box, width = 10, height = 6, device = cairo_pdf)
ggsave(file.path(fig_dir, "nutrient_layer_boxplots.svg"), p_box, width = 10, height = 6, device = svg)
ggsave(file.path(fig_dir, "nutrient_layer_boxplots.png"), p_box, width = 10, height = 6, dpi = 450, bg = "white")

# Spearman heatmap.
heat_df <- cor_long %>%
  mutate(
    env_label = factor(env_label, levels = rev(label_env[env_vars])),
    metric_label = factor(metric_label, levels = label_k[k_vars]),
    label = paste0(sprintf("%.2f", rho), star)
  )
p_heat <- ggplot(heat_df, aes(metric_label, env_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), family = "Times", size = 3.0) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0, limits = c(-1, 1)) +
  theme_pub(10) +
  labs(
    title = "IS nutrient-DMSOR Spearman correlations",
    subtitle = "Values are Spearman rho; * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001, † nominal p < 0.05",
    x = NULL,
    y = NULL,
    fill = "rho"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.subtitle = element_text(size = 8.5, hjust = 0.5))
ggsave(file.path(fig_dir, "nutrient_K08356_spearman_heatmap.pdf"), p_heat, width = 9.2, height = 5.8, device = cairo_pdf)
ggsave(file.path(fig_dir, "nutrient_K08356_spearman_heatmap.svg"), p_heat, width = 9.2, height = 5.8, device = svg)
ggsave(file.path(fig_dir, "nutrient_K08356_spearman_heatmap.png"), p_heat, width = 9.2, height = 5.8, dpi = 450, bg = "white")

# Selected scatterplots.
scatter_plots <- list()
for (i in seq_len(nrow(selected_pairs))) {
  pair <- selected_pairs[i, ]
  e <- pair$Environmental_variable
  k <- pair$K08356_metric
  label <- sprintf("rho = %.2f, p = %.3g, FDR = %.3g", pair$rho, pair$p_value, pair$FDR)
  scatter_plots[[i]] <- ggplot(merged, aes(.data[[e]], .data[[k]], color = Layer)) +
    geom_point(size = 2.4, alpha = 0.92) +
    geom_smooth(method = "lm", se = TRUE, color = "grey35", linewidth = 0.45, linetype = "dashed") +
    scale_color_manual(values = c("0-4" = "#74ADD1", "4-8" = "#66BD63", "8-12" = "#F46D43")) +
    theme_pub(9.5) +
    labs(
      title = paste0(label_k[k], " ~ ", label_env[e]),
      subtitle = label,
      x = label_env[e],
      y = paste0("log10(", sub("^log_", "", k), " + 1)"),
      color = "Layer"
    ) +
    theme(plot.subtitle = element_text(size = 8, hjust = 0.5))
}
p_scatter <- wrap_plots(scatter_plots, ncol = 2) +
  plot_annotation(
    title = "Selected nutrient-DMSOR associations within IS samples",
    caption = "Pairs were selected automatically from Spearman results by prioritizing high |rho| and lower p/FDR among non-canonical and total K08356 metrics.",
    theme = theme(
      plot.title = element_text(family = "Times", face = "bold", hjust = 0.5, size = 15),
      plot.caption = element_text(family = "Times", hjust = 0.5, size = 8.5)
    )
  )
ggsave(file.path(fig_dir, "selected_scatterplots.pdf"), p_scatter, width = 9.5, height = 7.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "selected_scatterplots.svg"), p_scatter, width = 9.5, height = 7.5, device = svg)
ggsave(file.path(fig_dir, "selected_scatterplots.png"), p_scatter, width = 9.5, height = 7.5, dpi = 450, bg = "white")

# Exploratory RDA with original nutrient variables only.
rda_env <- merged %>% select(NH4_mgkg, NO3_mgkg, TOP_mgkg, TOC_mgkg) %>% mutate(across(everything(), scale))
branch_mat <- merged %>%
  select(log_IdrA_MAG_TPM, log_canonical_aioA_MAG_TPM, log_aioA_like_MAG_TPM, log_unknown_DMSOR_MAG_TPM)
rda_fit <- rda(branch_mat ~ NH4_mgkg + NO3_mgkg + TOP_mgkg + TOC_mgkg, data = as.data.frame(rda_env))
rda_perm <- anova.cca(rda_fit, permutations = 999)
rda_axis <- anova.cca(rda_fit, by = "axis", permutations = 999)
rda_terms <- anova.cca(rda_fit, by = "term", permutations = 999)
write_csv(as.data.frame(rda_perm) %>% tibble::rownames_to_column("term"), file.path(result_dir, "RDA_overall_permutation.csv"))
write_csv(as.data.frame(rda_axis) %>% tibble::rownames_to_column("axis"), file.path(result_dir, "RDA_axis_permutation.csv"))
write_csv(as.data.frame(rda_terms) %>% tibble::rownames_to_column("term"), file.path(result_dir, "RDA_term_permutation.csv"))

site_scores <- as.data.frame(scores(rda_fit, display = "sites", choices = 1:2)) %>%
  mutate(Sample = merged$Sample, Layer = merged$Layer, Core = merged$Core)
species_scores <- as.data.frame(scores(rda_fit, display = "species", choices = 1:2)) %>%
  tibble::rownames_to_column("metric")
biplot_scores <- as.data.frame(scores(rda_fit, display = "bp", choices = 1:2)) %>%
  tibble::rownames_to_column("env")
write_csv(site_scores, file.path(result_dir, "RDA_site_scores.csv"))
write_csv(species_scores, file.path(result_dir, "RDA_branch_metric_scores.csv"))
write_csv(biplot_scores, file.path(result_dir, "RDA_environment_biplot_scores.csv"))

eig <- summary(rda_fit)$cont$importance
rda_label <- sprintf("Exploratory RDA: overall p = %.3g", as.data.frame(rda_perm)$`Pr(>F)`[1])
arrow_mult <- 1.0
p_rda <- ggplot(site_scores, aes(RDA1, RDA2, color = Layer)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.25) +
  geom_point(size = 2.6, alpha = 0.92) +
  geom_segment(
    data = biplot_scores,
    aes(x = 0, y = 0, xend = RDA1 * arrow_mult, yend = RDA2 * arrow_mult),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.16, "cm")),
    color = "grey25",
    linewidth = 0.45
  ) +
  geom_text(
    data = biplot_scores,
    aes(x = RDA1 * arrow_mult * 1.10, y = RDA2 * arrow_mult * 1.10, label = env),
    inherit.aes = FALSE,
    family = "Times",
    size = 3.2
  ) +
  scale_color_manual(values = c("0-4" = "#74ADD1", "4-8" = "#66BD63", "8-12" = "#F46D43")) +
  theme_pub(10) +
  labs(
    title = "Exploratory RDA of K08356 branch abundance constrained by nutrients",
    subtitle = rda_label,
    x = sprintf("RDA1 (%.1f%%)", eig["Proportion Explained", "RDA1"] * 100),
    y = sprintf("RDA2 (%.1f%%)", eig["Proportion Explained", "RDA2"] * 100),
    color = "Layer"
  ) +
  theme(plot.subtitle = element_text(size = 8.5, hjust = 0.5))
ggsave(file.path(fig_dir, "exploratory_RDA_nutrient_branch_abundance.pdf"), p_rda, width = 6.8, height = 5.4, device = cairo_pdf)
ggsave(file.path(fig_dir, "exploratory_RDA_nutrient_branch_abundance.svg"), p_rda, width = 6.8, height = 5.4, device = svg)
ggsave(file.path(fig_dir, "exploratory_RDA_nutrient_branch_abundance.png"), p_rda, width = 6.8, height = 5.4, dpi = 450, bg = "white")

readme <- c(
  "# IS internal nutrient-DMSOR niche association analysis",
  "",
  "This analysis links nutrient/organic matter gradients with K08356-related DMSOR gene- and MAG-level abundance within IS samples only.",
  "",
  "## Boundary",
  "- This is an IS-internal microenvironment association analysis.",
  "- It must not be interpreted as a four-habitat environmental driver analysis because AS/ES/NS nutrient data are not available here.",
  "",
  "## Data",
  "- Nutrient table: `/Users/catherine/Downloads/营养盐热图绘制/varechem.csv`.",
  "- MAG branch TPM: `K08356_branch_host_MAG_TPM_by_sample.csv`.",
  "- Gene-level K08356 TPM: `K08356_gene_total_sample_long.csv`.",
  sprintf("- Sample matching: %d nutrient samples; %d matched to gene TPM; %d matched to MAG TPM.", match_check$n[1], match_check$n[2], match_check$n[3]),
  "",
  "## Derived nutrient variables",
  "- DIN = NH4_mgkg + NO3_mgkg.",
  "- NH4_NO3_ratio = NH4_mgkg / NO3_mgkg.",
  "- TOC_TOP_ratio = TOC_mgkg / TOP_mgkg.",
  "",
  "## Abundance transformation",
  "- All K08356 abundance metrics were transformed as log10(TPM + 1) for correlation plots.",
  "",
  "## Figure interpretation",
  "- `nutrient_layer_boxplots`: IS nutrient and organic matter variation among sediment layers.",
  "- `nutrient_K08356_spearman_heatmap`: Spearman rho between nutrient variables and log-transformed K08356 abundance metrics.",
  "- `selected_scatterplots`: automatically selected nutrient-DMSOR pairs with relatively high |rho| and lower p/FDR, prioritizing non-canonical and total K08356 metrics.",
  "- `exploratory_RDA_nutrient_branch_abundance`: exploratory constrained ordination using NH4, NO3, TOP and TOC only.",
  "",
  "## Recommended wording",
  "These results should be described as nutrient-DMSOR niche association within IS samples, supporting micro-scale environmental heterogeneity potentially related to non-canonical DMSOR-bearing MAG distributions. They should not be described as causal evidence or as drivers of four-habitat differences."
)
writeLines(readme, file.path(out_dir, "README_interpretation.md"))

message("Wrote IS nutrient-DMSOR association outputs to: ", out_dir)
message("Matched samples: ", match_check$n[2], " gene TPM; ", match_check$n[3], " MAG TPM out of ", match_check$n[1], " nutrient samples.")

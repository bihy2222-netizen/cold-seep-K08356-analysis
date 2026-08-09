suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "Figure4_predefined_module_reanalysis_17_vs_28")
assignment_file <- file.path(base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected",
                             "Rhodobacteraceae_MAG_group_assignments_corrected.csv")
coverage_file <- file.path(out_dir, "predefined_module_MAG_coverage_45x17.csv")
ranked_file <- file.path(out_dir, "predefined_module_statistics_ranked.csv")
negative_17v23_file <- file.path(base_dir, "K08356_negative_control_module_analysis_17_5_23",
                                 "module_statistics_Clade4_17_vs_K08356negative_23.csv")
primary_is_17v18_file <- file.path(base_dir, "four_layer_IS_priority_module_analysis_17_5_18_5",
                                   "primary_IS_Clade4_17_vs_IS_K08356negative_18.csv")
abundance_file <- file.path(base_dir, "mag在各个样品中的丰度.xls")
sample_info_file <- "/Users/catherine/Downloads/Rstudio/sample_info.csv"

assignments <- read.csv(assignment_file, check.names = FALSE, stringsAsFactors = FALSE)
coverage <- read.csv(coverage_file, check.names = FALSE, stringsAsFactors = FALSE)
ranked <- read.csv(ranked_file, check.names = FALSE, stringsAsFactors = FALSE)
negative_17v23 <- read.csv(negative_17v23_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  select(-module) %>% rename(module = module_definition_name)
primary_is_17v18 <- read.csv(primary_is_17v18_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  select(-module) %>% rename(module = module_definition_name)

selected_modules <- ranked %>%
  filter(main_significant) %>%
  arrange(desc(robust_significant), main_min_fdr, desc(abs(cliffs_delta))) %>%
  pull(module)
if (length(selected_modules) != 5) stop("Figure 4 requires the five main-analysis significant modules.")

fixed_clade4 <- c(
  "SY365BB-8-12_bin8", "SY457BB-8-12_bin4", "SY459WG-0-4_bin37",
  "SY368YW-8-12_bin33", "SY456YB-8-12_bin34", "SY459WG-8-12_bin2",
  "SY365BB-8-12_bin10", "SY366YW-8-12_bin2", "SY456YB-0-4_bin8",
  "SY457BB-8-12_bin15", "SY457BB-4-8_bin29", "SY457BB-8-12_bin25",
  "SY456YB-8-12_bin27", "SY457BB-0-4_bin35", "SY457BB-0-4_bin6",
  "SY365BB-8-12_bin26", "SY457BB-4-8_bin33"
)
mag_order <- c(
  fixed_clade4,
  assignments %>% filter(comparison_group == "non-Clade 4") %>% arrange(genus, MAG) %>% pull(MAG)
)
if (length(mag_order) != 45 || anyDuplicated(mag_order) || !identical(mag_order[1:17], fixed_clade4)) {
  stop("45-MAG axis QC failed.")
}

module_labels <- c(
  "B12 biosynthesis" = "B12 biosynthesis",
  "TCA/reductive TCA central-cycle coverage" = "TCA / reductive TCA central cycle",
  "Dissimilatory sulfate/sulfite reduction" = "Dsr/Sat-associated sulfur-redox module",
  "Molybdopterin cofactor biosynthesis" = "Molybdopterin cofactor biosynthesis",
  "Denitrification" = "Denitrification"
)
module_order <- selected_modules
module_display_order <- unname(module_labels[module_order])

group_colors <- c("Clade 4 carrier" = "#7651A8", "non-Clade 4" = "#2E8B57")
category_colors <- c("Cofactor" = "#009E73", "Carbon" = "#4C78A8",
                     "Sulfur" = "#D99B2B", "Nitrogen" = "#7A5195")
habitat_colors <- c("IS" = "#E69F00", "AS" = "#56B4E9", "ES" = "#009E73", "NS" = "#CC79A7")

plot_coverage <- coverage %>%
  filter(module %in% selected_modules) %>%
  mutate(
    MAG = factor(MAG, levels = mag_order),
    module_label = unname(module_labels[module]),
    module_label = factor(module_label, levels = rev(module_display_order))
  )

annotation_data <- assignments %>%
  filter(MAG %in% mag_order) %>%
  mutate(MAG = factor(MAG, levels = mag_order))

p_annotation <- ggplot(annotation_data, aes(MAG, y = 1, fill = comparison_group)) +
  geom_tile(width = 0.98, height = 0.75, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 1.8) +
  annotate("text", x = 8.5, y = 1.75, label = "Clade 4 (n=17)",
           fontface = "bold", size = 3.6, color = group_colors[["Clade 4 carrier"]]) +
  annotate("text", x = 32, y = 1.75, label = "non-Clade 4 (n=28)",
           fontface = "bold", size = 3.6, color = group_colors[["non-Clade 4"]]) +
  scale_fill_manual(values = group_colors) +
  coord_cartesian(ylim = c(0.55, 2.0), clip = "off") +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(14, 74, 0, 238))

robust_lookup <- ranked %>%
  filter(module %in% selected_modules) %>%
  transmute(module_label = unname(module_labels[module]), robust_significant,
            status = if_else(robust_significant, "Robust in IS-only", "Main analysis only"))

p_heatmap <- ggplot(plot_coverage, aes(MAG, module_label, fill = module_coverage)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 2.1) +
  scale_fill_gradientn(colors = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                       limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                       name = "Module coverage (%)") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, title = "A  Predefined-module coverage across 45 Rhodobacteraceae MAGs") +
  theme_minimal(base_size = 9.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.8, color = "black"),
    axis.text.y = element_text(size = 8.6, color = "black"),
    plot.title = element_text(face = "bold", size = 11.5),
    legend.position = "right", legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.5), plot.margin = margin(3, 6, 5, 6)
  )

panel_a <- p_annotation / p_heatmap + plot_layout(heights = c(0.17, 0.83))

set.seed(20260805)
bootstrap_difference <- function(left, right, iterations = 5000) {
  estimates <- replicate(iterations,
                         mean(sample(left, length(left), replace = TRUE)) -
                           mean(sample(right, length(right), replace = TRUE)))
  tibble(ci_low = unname(quantile(estimates, 0.025)),
         ci_high = unname(quantile(estimates, 0.975)))
}

stability_wide <- ranked %>%
  filter(module %in% selected_modules) %>%
  select(module, category, difference_17v28 = coverage_difference,
         fdr_17v28 = coverage_fdr) %>%
  left_join(negative_17v23 %>%
              select(module, difference_17v23 = coverage_difference,
                     fdr_17v23 = coverage_BH_FDR), by = "module") %>%
  left_join(primary_is_17v18 %>%
              select(module, difference_IS17v18 = coverage_difference,
                     fdr_IS17v18 = coverage_BH_FDR), by = "module") %>%
  mutate(
    direction_consistent = sign(difference_17v28) == sign(difference_17v23) &
      sign(difference_17v28) == sign(difference_IS17v18),
    robust_all_three = direction_consistent & fdr_17v28 < 0.05 & fdr_17v23 < 0.05 &
      fdr_IS17v18 < 0.05
  )

coverage_ci <- primary_is_17v18 %>%
  filter(module %in% selected_modules) %>%
  transmute(
    module, category,
    coverage_difference,
    ci_low = coverage_bootstrap_95CI_low,
    ci_high = coverage_bootstrap_95CI_high,
    coverage_fdr = coverage_BH_FDR
  ) %>%
  left_join(stability_wide %>% select(module, robust_all_three), by = "module") %>%
  mutate(module_label = factor(unname(module_labels[module]), levels = rev(module_display_order)),
         status = if_else(robust_all_three, "Stable across all three", "Not stable in all three"),
         fdr_label = paste0("FDR=", formatC(coverage_fdr, format = "g", digits = 2)))

p_forest <- ggplot(coverage_ci, aes(coverage_difference, module_label, color = category)) +
  geom_vline(xintercept = 0, color = "#555555", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.16, linewidth = 0.8,
                orientation = "y") +
  geom_point(aes(shape = status), size = 3.1, stroke = 0.8) +
  geom_text(aes(x = ci_high + 1.2, label = fdr_label), hjust = 0, size = 2.65, color = "#333333") +
  scale_color_manual(values = category_colors, guide = "none") +
  scale_shape_manual(values = c("Stable across all three" = 16, "Not stable in all three" = 1), name = NULL) +
  coord_cartesian(clip = "off") +
  labs(x = "Coverage difference (IS Clade 4 − IS K08356-negative, pp)", y = NULL,
       title = "B  Core IS-only coverage differences and bootstrap 95% CI") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 7.5, color = "black"),
        plot.title = element_text(face = "bold", size = 10.5), legend.position = "bottom",
        plot.margin = margin(5, 42, 5, 5))

prevalence_ci <- primary_is_17v18 %>%
  filter(module %in% selected_modules) %>%
  transmute(
    module, category,
    prevalence_difference_pp,
    ci_low = prevalence_bootstrap_95CI_low_pp,
    ci_high = prevalence_bootstrap_95CI_high_pp,
    fisher_fdr = fisher_BH_FDR
  ) %>%
  mutate(module_label = factor(unname(module_labels[module]), levels = rev(module_display_order)),
         fdr_label = paste0("Fisher FDR=", formatC(fisher_fdr, format = "g", digits = 2)))

p_prevalence <- ggplot(prevalence_ci, aes(prevalence_difference_pp, module_label, color = category)) +
  geom_vline(xintercept = 0, color = "#555555", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.16, linewidth = 0.8,
                orientation = "y") +
  geom_point(size = 3.0) +
  geom_text(aes(x = ci_high + 1.5, label = fdr_label), hjust = 0, size = 2.55, color = "#333333") +
  scale_color_manual(values = category_colors, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = "Prevalence difference at ≥75% coverage (pp)", y = NULL,
       title = "C  Core IS-only module-prevalence differences",
       subtitle = "Only TCA/rTCA passed Fisher BH-FDR < 0.05") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 7.5, color = "black"),
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 7.3, color = "#555555"),
        plot.margin = margin(5, 48, 5, 5))

sensitivity <- stability_wide %>%
  select(module, category,
         `17 vs 28` = difference_17v28,
         `17 vs 23` = difference_17v23,
         `IS 17 vs 18` = difference_IS17v18,
         fdr_17v28, fdr_17v23, fdr_IS17v18) %>%
  pivot_longer(c(`17 vs 28`, `17 vs 23`, `IS 17 vs 18`),
               names_to = "analysis", values_to = "coverage_difference") %>%
  mutate(
         coverage_fdr = case_when(
           analysis == "17 vs 28" ~ fdr_17v28,
           analysis == "17 vs 23" ~ fdr_17v23,
           TRUE ~ fdr_IS17v18
         ),
         coverage_significant = coverage_fdr < 0.05,
         analysis = factor(analysis, levels = c("17 vs 28", "17 vs 23", "IS 17 vs 18")),
         module_label = unname(module_labels[module]),
         label_offset = case_when(
           module == "TCA/reductive TCA central-cycle coverage" ~ 1.2,
           module == "Denitrification" ~ -1.0,
           module == "Molybdopterin cofactor biosynthesis" ~ 1.0,
           module == "Dissimilatory sulfate/sulfite reduction" ~ -0.8,
           TRUE ~ -0.4
         ),
         label_y = coverage_difference + label_offset)

p_sensitivity <- ggplot(sensitivity,
                         aes(analysis, coverage_difference, group = module_label, color = category)) +
  geom_hline(yintercept = 0, color = "#666666", linewidth = 0.4) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(shape = coverage_significant), size = 2.8, stroke = 0.8) +
  geom_text(data = sensitivity %>% filter(analysis == "IS 17 vs 18"),
            aes(y = label_y, label = module_label), hjust = -0.08, size = 2.35, show.legend = FALSE) +
  scale_color_manual(values = category_colors, name = "Module category") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "Coverage FDR < 0.05", `FALSE` = "FDR ≥ 0.05"),
                     name = "Per-comparison status") +
  coord_cartesian(xlim = c(0.85, 3.75), clip = "off") +
  labs(x = NULL, y = "Coverage difference (pp)", title = "D  Three-level sensitivity") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.5),
        legend.position = "bottom", legend.box = "vertical", plot.margin = margin(5, 96, 5, 5))

abundance <- read_excel(abundance_file, sheet = "merged_tpm_matrix", .name_repair = "unique")
names(abundance)[1] <- "MAG"
sample_columns <- names(abundance)[-1]
if (length(sample_columns) != 56) stop("TPM workbook must contain 56 sample columns.")

sample_metadata_raw <- read.csv(sample_info_file, check.names = FALSE, stringsAsFactors = FALSE)
names(sample_metadata_raw)[names(sample_metadata_raw) == ""] <- "row_id"
sample_metadata <- sample_metadata_raw %>%
  transmute(Sample = Sample, habitat = Group,
            Sample = str_replace(Sample, "^SY457YW", "SY457BB"),
            Sample = str_replace(Sample, "^SY459YW", "SY459WG"),
            Sample = str_replace(Sample, "^SQ_58_-8-12$", "SQ_58_8-12"),
            Sample = str_replace(Sample, "^SQ_81_-8-12$", "SQ_81_8-12"))
if (!setequal(sample_columns, sample_metadata$Sample)) stop("56-sample metadata matching failed.")

abundance_long <- abundance %>%
  filter(MAG %in% mag_order) %>%
  pivot_longer(all_of(sample_columns), names_to = "Sample", values_to = "TPM")

weighted_potential <- abundance_long %>%
  inner_join(coverage %>% filter(module %in% selected_modules) %>%
               select(MAG, module, comparison_group, module_coverage), by = "MAG",
             relationship = "many-to-many") %>%
  mutate(module_weighted_TPM = TPM * module_coverage / 100) %>%
  group_by(Sample, module, comparison_group) %>%
  summarise(total_module_weighted_TPM = sum(module_weighted_TPM), .groups = "drop") %>%
  left_join(sample_metadata, by = "Sample") %>%
  mutate(module_label = factor(unname(module_labels[module]), levels = module_display_order),
         habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS")),
         comparison_group = factor(comparison_group, levels = c("Clade 4 carrier", "non-Clade 4")))
if (nrow(weighted_potential) != 56 * 5 * 2 || any(is.na(weighted_potential$habitat))) {
  stop("56-sample weighted-potential QC failed.")
}
write.csv(weighted_potential, file.path(out_dir, "Figure4_final_v2_abundance_weighted_module_potential_56samples.csv"), row.names = FALSE)

p_weighted <- ggplot(weighted_potential,
                     aes(habitat, log10(total_module_weighted_TPM + 1), fill = comparison_group)) +
  geom_boxplot(position = position_dodge(width = 0.72), width = 0.62, outlier.shape = NA,
               alpha = 0.78, linewidth = 0.4) +
  geom_point(aes(color = comparison_group),
             position = position_jitterdodge(jitter.width = 0.13, dodge.width = 0.72),
             size = 0.75, alpha = 0.50, show.legend = FALSE) +
  facet_wrap(~module_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors, name = "MAG group") +
  scale_color_manual(values = group_colors) +
  labs(x = "Sample habitat",
       y = expression(log[10](1 + Sigma~TPM %*% (coverage/100))),
       title = "E  Abundance-weighted genomic module potential across 56 metagenomes") +
  theme_minimal(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 7.5),
        plot.title = element_text(face = "bold", size = 10.8), legend.position = "bottom",
        axis.text.x = element_text(face = "bold"))

right_column <- p_forest / p_prevalence / p_sensitivity + plot_layout(heights = c(1, 1, 1.08))
top_row <- (panel_a | right_column) + plot_layout(widths = c(1.42, 1))
figure4 <- top_row / p_weighted +
  plot_layout(heights = c(1.95, 0.72)) +
  plot_annotation(
    title = "Predefined metabolic-module profiles within Rhodobacteraceae",
    subtitle = "The core IS-only 17 vs 18 analysis identifies four coverage-level differences that remain stable across three nested comparisons",
    caption = str_wrap(paste0(
      "Panels B and C show the core IS-only comparison of 17 Clade 4 carriers versus 18 K08356-negative MAGs. ",
      "Panel D compares 17 vs 28, 17 vs 23 and IS-only 17 vs 18; filled points denote coverage BH-FDR < 0.05 in each comparison. ",
      "Four modules remain directionally consistent and coverage-FDR significant across all three comparisons, whereas denitrification loses significance in IS-only analysis. ",
      "Prevalence is defined a priori as ≥75% module coverage; only TCA/rTCA passes Fisher BH-FDR < 0.05 in the core comparison. ",
      "Abundance-weighted potential = Σ[MAG TPM × (module coverage/100)], so coverage enters as a 0–1 fraction. ",
      "B12, TCA/rTCA, Moco and Dsr/Sat-associated results indicate higher coverage, not enrichment of complete active pathways. ",
      "The five other K08356-positive MAGs are a descriptive third group and show no pairwise FDR-significant module differences. ",
      "Results may still reflect MAG completeness and taxonomic structure and should not be interpreted as pathway activity without additional validation."
    ), width = 190),
    theme = theme(plot.title = element_text(face = "bold", size = 17),
                  plot.subtitle = element_text(size = 10.5, color = "#404040"),
                  plot.caption = element_text(size = 8.1, color = "#444444", hjust = 0))
  )

save_plot <- function(extension, device, dpi = 450) {
  args <- list(filename = file.path(out_dir, paste0("Figure4_predefined_modules_final_v2.", extension)),
               plot = figure4, width = 19, height = 15.5, units = "in", bg = "white", device = device)
  if (extension == "png") args$dpi <- dpi
  do.call(ggsave, args)
}
save_plot("pdf", cairo_pdf)
save_plot("svg", grDevices::svg)
save_plot("png", "png")

ci_output <- coverage_ci %>%
  transmute(module, category, coverage_difference, bootstrap_95CI_low = ci_low,
            bootstrap_95CI_high = ci_high, coverage_fdr, robust_all_three)
write.csv(ci_output, file.path(out_dir, "Figure4_final_v2_core_IS_coverage_bootstrap_CI.csv"), row.names = FALSE)

figure_qc <- c(
  "Figure 4 predefined-module QC",
  "Status: PASSED",
  sprintf("Figure A MAG columns: %d", n_distinct(plot_coverage$MAG)),
  sprintf("Figure A modules: %d", n_distinct(plot_coverage$module)),
  sprintf("Modules stable across 17v28, 17v23 and IS17v18: %d", sum(stability_wide$robust_all_three)),
  sprintf("Modules not stable across all three: %d", sum(!stability_wide$robust_all_three)),
  "Panel B/C core comparison: IS Clade 4 n=17 vs IS K08356-negative n=18.",
  "Panel D comparison columns: 17 vs 28; 17 vs 23; IS-only 17 vs 18.",
  sprintf("Bootstrap iterations per module: %d", 5000),
  sprintf("Weighted-potential rows: %d (required 560)", nrow(weighted_potential)),
  sprintf("Weighted-potential samples: %d (required 56)", n_distinct(weighted_potential$Sample)),
  "Sulfur-module display name: Dsr/Sat-associated sulfur-redox module.",
  "E-panel formula: sum(MAG TPM x module_coverage/100); coverage is used as a 0-1 fraction.",
  sprintf("Core IS-only modules with Fisher BH-FDR < 0.05: %d",
          sum(primary_is_17v18$fisher_BH_FDR < 0.05)),
  "K08355/K08356 do not contribute to any displayed module score.",
  "QC PASSED refers only to data integrity, grouping, and computation.",
  "MAG completeness and taxonomic structure are not fully adjusted.",
  "No claim of measured pathway activity or complete-module enrichment is made."
)
writeLines(figure_qc, file.path(out_dir, "Figure4_predefined_modules_final_v2_QC.txt"))
cat(paste(figure_qc, collapse = "\n"), "\n")

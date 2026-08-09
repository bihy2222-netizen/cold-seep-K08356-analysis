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
four_layer_dir <- file.path(base_dir, "four_layer_IS_priority_module_analysis_17_5_18_5")
negative_dir <- file.path(base_dir, "K08356_negative_control_module_analysis_17_5_23")

assignments <- read.csv(
  file.path(base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected",
            "Rhodobacteraceae_MAG_group_assignments_corrected.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  mutate(
    IS_three_group = case_when(
      comparison_group == "Clade 4 carrier" & source_habitat == "IS" ~ "Clade 4",
      comparison_group == "non-Clade 4" & K08356_status == "K08356-positive" &
        source_habitat == "IS" ~ "other K08356-positive",
      K08356_status == "K08356-negative" & source_habitat == "IS" ~ "K08356-negative",
      TRUE ~ NA_character_
    )
  )
coverage <- read.csv(file.path(out_dir, "predefined_module_MAG_coverage_45x17.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE) %>%
  select(-any_of(c("comparison_group", "source_habitat", "genus"))) %>%
  left_join(assignments %>% select(MAG, IS_three_group, source_habitat, genus), by = "MAG")
ranked <- read.csv(file.path(out_dir, "predefined_module_statistics_ranked.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
negative_17v23 <- read.csv(
  file.path(negative_dir, "module_statistics_Clade4_17_vs_K08356negative_23.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
) %>% select(-module) %>% rename(module = module_definition_name)
primary_17v18 <- read.csv(
  file.path(four_layer_dir, "primary_IS_Clade4_17_vs_IS_K08356negative_18.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
) %>% select(-module) %>% rename(module = module_definition_name)

robust_modules <- c(
  "B12 biosynthesis",
  "TCA/reductive TCA central-cycle coverage",
  "Dissimilatory sulfate/sulfite reduction",
  "Molybdopterin cofactor biosynthesis"
)
stability_modules <- c(robust_modules, "Denitrification")
module_labels <- c(
  "B12 biosynthesis" = "B12 biosynthesis",
  "TCA/reductive TCA central-cycle coverage" = "TCA / reductive TCA central cycle",
  "Dissimilatory sulfate/sulfite reduction" = "Dsr/Sat-associated sulfur-redox module",
  "Molybdopterin cofactor biosynthesis" = "Molybdopterin cofactor biosynthesis",
  "Denitrification" = "Denitrification"
)
module_order <- robust_modules
module_display_order <- unname(module_labels[module_order])

group_levels <- c("Clade 4", "other K08356-positive", "K08356-negative")
group_colors <- c("Clade 4" = "#7651A8", "other K08356-positive" = "#E69F00",
                  "K08356-negative" = "#2E8B57")
category_colors <- c("Cofactor" = "#009E73", "Carbon" = "#4C78A8",
                     "Sulfur" = "#D99B2B", "Nitrogen" = "#8C6BB1")

fixed_clade4 <- c(
  "SY365BB-8-12_bin8", "SY457BB-8-12_bin4", "SY459WG-0-4_bin37",
  "SY368YW-8-12_bin33", "SY456YB-8-12_bin34", "SY459WG-8-12_bin2",
  "SY365BB-8-12_bin10", "SY366YW-8-12_bin2", "SY456YB-0-4_bin8",
  "SY457BB-8-12_bin15", "SY457BB-4-8_bin29", "SY457BB-8-12_bin25",
  "SY456YB-8-12_bin27", "SY457BB-0-4_bin35", "SY457BB-0-4_bin6",
  "SY365BB-8-12_bin26", "SY457BB-4-8_bin33"
)
other_order <- assignments %>%
  filter(IS_three_group == "other K08356-positive") %>% arrange(genus, MAG) %>% pull(MAG)
negative_order <- assignments %>%
  filter(IS_three_group == "K08356-negative") %>% arrange(genus, MAG) %>% pull(MAG)
mag_order <- c(fixed_clade4, other_order, negative_order)

group_counts <- table(assignments$IS_three_group, useNA = "no")
mandatory_qc <- c(
  length(mag_order) == 40,
  length(unique(mag_order)) == 40,
  identical(mag_order[1:17], fixed_clade4),
  group_counts[["Clade 4"]] == 17,
  group_counts[["other K08356-positive"]] == 5,
  group_counts[["K08356-negative"]] == 18,
  all(robust_modules %in% coverage$module),
  all(primary_17v18$coverage_BH_FDR[match(robust_modules, primary_17v18$module)] < 0.05)
)
if (!all(mandatory_qc)) stop("Final v3 core grouping/module QC failed.")

plot_coverage <- coverage %>%
  filter(MAG %in% mag_order, module %in% robust_modules) %>%
  mutate(
    MAG = factor(MAG, levels = mag_order),
    IS_three_group = factor(IS_three_group, levels = group_levels),
    module_label = factor(unname(module_labels[module]), levels = rev(module_display_order))
  )
annotation_data <- assignments %>%
  filter(MAG %in% mag_order) %>%
  mutate(MAG = factor(MAG, levels = mag_order),
         IS_three_group = factor(IS_three_group, levels = group_levels))

p_annotation <- ggplot(annotation_data, aes(MAG, y = 1, fill = IS_three_group)) +
  geom_tile(width = 0.98, height = 0.72, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(17.5, 22.5), color = "white", linewidth = 1.8) +
  annotate("text", x = 8.8, y = 1.72, label = "Clade 4 (n=17)",
           fontface = "bold", size = 3.6, color = group_colors[["Clade 4"]]) +
  annotate("text", x = 20, y = 1.72, label = "other K+ (n=5)",
           fontface = "bold", size = 3.2, color = group_colors[["other K08356-positive"]]) +
  annotate("text", x = 31.5, y = 1.72, label = "K08356-negative (n=18)",
           fontface = "bold", size = 3.5, color = group_colors[["K08356-negative"]]) +
  scale_fill_manual(values = group_colors) +
  coord_cartesian(ylim = c(0.55, 1.95), clip = "off") +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(13, 55, 0, 230))

p_heatmap <- ggplot(plot_coverage, aes(MAG, module_label, fill = module_coverage)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_vline(xintercept = c(17.5, 22.5), color = "white", linewidth = 2.0) +
  scale_fill_gradientn(colors = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                       limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                       name = "Module coverage (%)") +
  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, title = "A  Four reproducible module-coverage profiles in IS MAGs") +
  theme_minimal(base_size = 9.2) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.7, color = "black"),
        axis.text.y = element_text(size = 8.4, color = "black"),
        plot.title = element_text(face = "bold", size = 11.2),
        legend.position = "right", legend.title = element_text(face = "bold", size = 8.3),
        legend.text = element_text(size = 7.3), plot.margin = margin(3, 5, 5, 5))
panel_a <- p_annotation / p_heatmap + plot_layout(heights = c(0.18, 0.82))

primary_labels <- primary_17v18 %>%
  filter(module %in% robust_modules) %>%
  transmute(module_label = unname(module_labels[module]),
            label = paste0("FDR=", formatC(coverage_BH_FDR, format = "g", digits = 2)))

p_distribution <- ggplot(plot_coverage,
                          aes(IS_three_group, module_coverage, fill = IS_three_group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.72, linewidth = 0.45) +
  geom_jitter(aes(color = IS_three_group), width = 0.12, size = 1.25, alpha = 0.72,
              show.legend = FALSE) +
  geom_text(data = primary_labels, aes(x = 2, y = 104, label = label),
            inherit.aes = FALSE, size = 2.5, color = "#333333") +
  facet_wrap(
    ~module_label,
    ncol = 2,
    labeller = as_labeller(c(
      "B12 biosynthesis" = "B12 biosynthesis",
      "TCA / reductive TCA central cycle" = "TCA/rTCA central cycle",
      "Dsr/Sat-associated sulfur-redox module" = "Dsr/Sat sulfur-redox",
      "Molybdopterin cofactor biosynthesis" = "Moco biosynthesis"
    ))
  ) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = c("Clade 4" = "Clade 4\n(n=17)",
                              "other K08356-positive" = "other K+\n(n=5)",
                              "K08356-negative" = "K-negative\n(n=18)")) +
  coord_cartesian(ylim = c(0, 108), clip = "off") +
  labs(x = NULL, y = "Module coverage (%)",
       title = "B  IS three-group distributions",
       subtitle = "FDR refers to the core 17 vs 18 contrast") +
  theme_minimal(base_size = 8.4) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 7.3),
        axis.text.x = element_text(size = 6.8), plot.title = element_text(face = "bold", size = 10.2),
        plot.subtitle = element_text(size = 7.7, color = "#444444"),
        legend.position = "bottom")

forest <- primary_17v18 %>%
  filter(module %in% robust_modules) %>%
  transmute(module, category, coverage_difference,
            ci_low = coverage_bootstrap_95CI_low,
            ci_high = coverage_bootstrap_95CI_high,
            coverage_fdr = coverage_BH_FDR,
            module_label = factor(unname(module_labels[module]), levels = rev(module_display_order)),
            fdr_label = paste0("FDR=", formatC(coverage_fdr, format = "g", digits = 2)))

p_forest <- ggplot(forest, aes(coverage_difference, module_label, color = category)) +
  geom_vline(xintercept = 0, color = "#555555", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.16, linewidth = 0.8,
                orientation = "y") +
  geom_point(size = 3.0) +
  geom_text(aes(x = ci_high + 1.0, label = fdr_label), hjust = 0, size = 2.55, color = "#333333") +
  scale_color_manual(values = category_colors, guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = "Coverage difference (IS Clade 4 − IS K-negative, pp)", y = NULL,
       title = "C  Core IS-only effects and bootstrap 95% CI") +
  theme_minimal(base_size = 8.7) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 7.3, color = "black"),
        plot.title = element_text(face = "bold", size = 10.2), plot.margin = margin(5, 40, 5, 5))

stability <- ranked %>%
  filter(module %in% stability_modules) %>%
  select(module, category, difference_17v28 = coverage_difference, fdr_17v28 = coverage_fdr) %>%
  left_join(negative_17v23 %>%
              select(module, difference_17v23 = coverage_difference, fdr_17v23 = coverage_BH_FDR),
            by = "module") %>%
  left_join(primary_17v18 %>%
              select(module, difference_IS17v18 = coverage_difference, fdr_IS17v18 = coverage_BH_FDR),
            by = "module") %>%
  pivot_longer(c(difference_17v28, difference_17v23, difference_IS17v18),
               names_to = "analysis", values_to = "coverage_difference") %>%
  mutate(
    coverage_fdr = case_when(
      analysis == "difference_17v28" ~ fdr_17v28,
      analysis == "difference_17v23" ~ fdr_17v23,
      TRUE ~ fdr_IS17v18
    ),
    significant = coverage_fdr < 0.05,
    analysis = factor(analysis,
                      levels = c("difference_17v28", "difference_17v23", "difference_IS17v18"),
                      labels = c("17 vs 28", "17 vs 23", "IS 17 vs 18")),
    module_label = unname(module_labels[module]),
    line_color = if_else(module == "Denitrification", "#9E9E9E", unname(category_colors[category])),
    endpoint_label = case_when(
      module == "TCA/reductive TCA central-cycle coverage" ~ "TCA/rTCA",
      module == "Denitrification" ~ "Denitrification (not robust)",
      module == "Molybdopterin cofactor biosynthesis" ~ "Moco",
      module == "Dissimilatory sulfate/sulfite reduction" ~ "Dsr/Sat sulfur-redox",
      TRUE ~ "B12 biosynthesis"
    ),
    label_y = case_when(
      module == "TCA/reductive TCA central-cycle coverage" ~ 26.2,
      module == "Denitrification" ~ 19.1,
      module == "Molybdopterin cofactor biosynthesis" ~ 16.7,
      module == "Dissimilatory sulfate/sulfite reduction" ~ 14.8,
      TRUE ~ 12.6
    )
  )

p_stability <- ggplot(stability,
                      aes(analysis, coverage_difference, group = module_label, color = line_color)) +
  geom_hline(yintercept = 0, color = "#666666", linewidth = 0.4) +
  geom_line(linewidth = 0.85, alpha = 0.88) +
  geom_point(aes(shape = significant), size = 2.8, stroke = 0.8) +
  geom_text(data = stability %>% filter(analysis == "IS 17 vs 18"),
            aes(y = label_y, label = endpoint_label), hjust = -0.08, size = 2.3,
            show.legend = FALSE) +
  scale_color_identity() +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "Coverage FDR < 0.05", `FALSE` = "FDR ≥ 0.05"),
                     name = NULL) +
  coord_cartesian(xlim = c(0.85, 3.75), ylim = c(0, 34), clip = "off") +
  labs(x = NULL, y = "Coverage difference (pp)",
       title = "D  Reproducibility across three control definitions") +
  theme_minimal(base_size = 8.7) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.2),
        legend.position = "bottom", plot.margin = margin(5, 95, 5, 5))

abundance_file <- file.path(base_dir, "mag在各个样品中的丰度.xls")
sample_info_file <- "/Users/catherine/Downloads/Rstudio/sample_info.csv"
abundance <- read_excel(abundance_file, sheet = "merged_tpm_matrix", .name_repair = "unique")
names(abundance)[1] <- "MAG"
sample_columns <- names(abundance)[-1]
if (length(sample_columns) != 56) stop("TPM workbook must contain 56 samples.")
sample_metadata_raw <- read.csv(sample_info_file, check.names = FALSE, stringsAsFactors = FALSE)
names(sample_metadata_raw)[names(sample_metadata_raw) == ""] <- "row_id"
sample_metadata <- sample_metadata_raw %>%
  transmute(Sample = Sample, habitat = Group,
            Sample = str_replace(Sample, "^SY457YW", "SY457BB"),
            Sample = str_replace(Sample, "^SY459YW", "SY459WG"),
            Sample = str_replace(Sample, "^SQ_58_-8-12$", "SQ_58_8-12"),
            Sample = str_replace(Sample, "^SQ_81_-8-12$", "SQ_81_8-12"))
if (!setequal(sample_columns, sample_metadata$Sample)) stop("56-sample metadata matching failed.")

weighted_potential <- abundance %>%
  filter(MAG %in% mag_order) %>%
  pivot_longer(all_of(sample_columns), names_to = "Sample", values_to = "TPM") %>%
  inner_join(coverage %>% filter(MAG %in% mag_order, module %in% robust_modules) %>%
               select(MAG, module, IS_three_group, module_coverage),
             by = "MAG", relationship = "many-to-many") %>%
  mutate(module_weighted_TPM = TPM * module_coverage / 100) %>%
  group_by(Sample, module, IS_three_group) %>%
  summarise(total_module_weighted_TPM = sum(module_weighted_TPM), .groups = "drop") %>%
  left_join(sample_metadata, by = "Sample") %>%
  mutate(
    module_label = factor(unname(module_labels[module]), levels = module_display_order),
    IS_three_group = factor(IS_three_group, levels = group_levels),
    habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS"))
  )
if (nrow(weighted_potential) != 56 * 4 * 3 || any(is.na(weighted_potential$habitat))) {
  stop("Final v3 abundance-weighted potential QC failed.")
}
write.csv(weighted_potential,
          file.path(out_dir, "Figure4_final_v3_abundance_weighted_potential_56samples.csv"),
          row.names = FALSE)

p_weighted <- ggplot(weighted_potential,
                     aes(habitat, log10(total_module_weighted_TPM + 1), fill = IS_three_group)) +
  geom_boxplot(position = position_dodge(width = 0.76), width = 0.66, outlier.shape = NA,
               alpha = 0.75, linewidth = 0.4) +
  geom_point(aes(color = IS_three_group),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.76),
             size = 0.68, alpha = 0.46, show.legend = FALSE) +
  facet_wrap(~module_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors, name = "IS-source MAG group") +
  scale_color_manual(values = group_colors) +
  labs(x = "Metagenome habitat",
       y = expression(log[10](1 + Sigma~TPM %*% (coverage/100))),
       title = "E  Abundance-weighted genomic potential of the four modules across 56 metagenomes") +
  theme_minimal(base_size = 8.3) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 7.4),
        plot.title = element_text(face = "bold", size = 10.6), legend.position = "bottom",
        axis.text.x = element_text(face = "bold"))

right_column <- p_distribution / p_forest / p_stability +
  plot_layout(heights = c(1.28, 0.82, 0.95))
top_row <- (panel_a | right_column) + plot_layout(widths = c(1.33, 1))
figure4_v3 <- top_row / p_weighted +
  plot_layout(heights = c(1.95, 0.70)) +
  plot_annotation(
    title = "Clade 4 carriers show reproducibly higher coverage of four metabolic modules",
    subtitle = "The B12, TCA/rTCA, Dsr/Sat-associated sulfur-redox and molybdopterin-cofactor modules remain significant across three nested control definitions",
    caption = str_wrap(paste0(
      "The core inference uses IS Clade 4 carriers (n=17) versus IS K08356-negative Rhodobacteraceae (n=18). ",
      "The five other K08356-positive IS MAGs are shown as a descriptive third group only. ",
      "Coverage differences were tested by Mann–Whitney tests with BH correction across 17 predefined modules; confidence intervals are based on 5,000 bootstrap resamples. ",
      "Filled points in panel D denote coverage BH-FDR < 0.05. Denitrification is retained only as a grey non-robust reference because it loses significance in IS-only analysis. ",
      "Higher coverage indicates that more predefined module steps are encoded; it does not mean that the modules are unique to Clade 4 or prove pathway activity. ",
      "Abundance-weighted potential = Σ[MAG TPM × (module coverage/100)]. MAG completeness and taxonomic structure remain potential limitations."
    ), width = 190),
    theme = theme(plot.title = element_text(face = "bold", size = 17),
                  plot.subtitle = element_text(size = 10.4, color = "#404040"),
                  plot.caption = element_text(size = 8.0, color = "#444444", hjust = 0))
  )

save_figure <- function(extension, device, dpi = 450) {
  args <- list(filename = file.path(out_dir, paste0("Figure4_module_coverage_story_final_v3.", extension)),
               plot = figure4_v3, width = 19, height = 15.7, units = "in", bg = "white",
               device = device)
  if (extension == "png") args$dpi <- dpi
  do.call(ggsave, args)
}
save_figure("pdf", cairo_pdf)
save_figure("svg", grDevices::svg)
save_figure("png", "png")

figure_qc <- c(
  "Figure 4 module-coverage story final v3 QC",
  "Status: PASSED",
  sprintf("Panel A IS MAGs: %d (17 + 5 + 18)", n_distinct(plot_coverage$MAG)),
  sprintf("Robust modules in A/B/C/E: %d", n_distinct(plot_coverage$module)),
  sprintf("Panel D stability modules: %d (four robust + denitrification reference)",
          n_distinct(stability$module)),
  sprintf("Four robust modules coverage-FDR significant in all three comparisons: %s",
          all(stability$significant[stability$module %in% robust_modules])),
  sprintf("Denitrification significant in IS 17 vs 18: %s",
          stability$significant[stability$module == "Denitrification" &
                                  stability$analysis == "IS 17 vs 18"]),
  sprintf("Weighted-potential rows: %d (required 672)", nrow(weighted_potential)),
  "Coverage is used as a 0-1 fraction in the TPM-weighted formula.",
  "No prevalence-enrichment or pathway-activity claim is made.",
  "Existing Figure 4 versions and original statistics were not overwritten."
)
writeLines(figure_qc, file.path(out_dir, "Figure4_module_coverage_story_final_v3_QC.txt"))
cat(paste(figure_qc, collapse = "\n"), "\n")

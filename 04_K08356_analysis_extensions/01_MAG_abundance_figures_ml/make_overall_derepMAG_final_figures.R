suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

base_out <- "/Users/catherine/Downloads/MAG-bin-tpm/K08356_gene_and_MAG_dual_evidence_20260708"
mag_out <- file.path(base_out, "02_MAG_host_TPM")
fig_out <- file.path(base_out, "overall_derepMAG_final_figures")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

habitat_levels <- c("IS", "AS", "ES", "NS")
hab_cols <- c(IS = "#2CA25F", AS = "#8E63C7", ES = "#E39D25", NS = "#3B78C8")
feature_levels <- c(
  "K08356-bearing host MAG total TPM",
  "aioA-like-associated",
  "unknown AioA-like / uncertain DMSOR",
  "IdrA-associated",
  "canonical aioA-associated"
)
feature_labels <- c(
  "K08356-bearing host MAG total TPM" = "Total K08356-bearing\nhost MAGs",
  "aioA-like-associated" = "aioA-like-associated",
  "unknown AioA-like / uncertain DMSOR" = "unknown AioA-like /\nuncertain DMSOR",
  "IdrA-associated" = "IdrA-associated",
  "canonical aioA-associated" = "canonical\naioA-associated"
)

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "KW p = NA",
    ifelse(p < 0.001, paste0("KW p = ", formatC(p, format = "e", digits = 2)),
           paste0("KW p = ", signif(p, 3)))
  )
}

theme_mag <- theme_classic(base_family = "Times", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5),
    axis.text = element_text(color = "grey15"),
    axis.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "grey35", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "top",
    legend.title = element_blank(),
    panel.spacing = unit(1.0, "lines"),
    plot.margin = margin(10, 12, 10, 10)
  )

total <- read_csv(
  file.path(mag_out, "K08356_unique_host_MAG_total_TPM_by_sample.csv"),
  show_col_types = FALSE
) %>%
  select(Feature, Sample, TPM, logTPM, Habitat, Sample_ordered)

branch <- read_csv(
  file.path(mag_out, "K08356_branch_host_MAG_TPM_by_sample.csv"),
  show_col_types = FALSE
) %>%
  select(Feature, Sample, TPM, logTPM, Habitat, Sample_ordered)

plot_dat <- bind_rows(total, branch) %>%
  mutate(
    Habitat = factor(Habitat, levels = habitat_levels),
    Feature = factor(Feature, levels = feature_levels)
  ) %>%
  filter(!is.na(Feature), !is.na(Habitat))

kw_total <- read_csv(
  file.path(mag_out, "K08356_unique_host_MAG_total_Kruskal.csv"),
  show_col_types = FALSE
) %>%
  select(Feature, p_value)

kw_branch <- read_csv(
  file.path(mag_out, "K08356_branch_host_MAG_Kruskal.csv"),
  show_col_types = FALSE
) %>%
  select(Feature, p_value)

kw_labels <- bind_rows(kw_total, kw_branch) %>%
  mutate(
    Feature = factor(Feature, levels = feature_levels),
    label = fmt_p(p_value),
    x = 2.5
  ) %>%
  filter(!is.na(Feature))

main_fig <- ggplot(plot_dat, aes(Habitat, logTPM, fill = Habitat)) +
  geom_boxplot(width = 0.58, alpha = 0.78, color = "grey20", outlier.shape = NA, linewidth = 0.42) +
  geom_jitter(aes(color = Habitat), width = 0.11, size = 1.75, alpha = 0.82, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white", color = "black") +
  geom_text(
    data = kw_labels,
    aes(x = x, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.35,
    size = 3.35,
    family = "Times"
  ) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 2, labeller = as_labeller(feature_labels)) +
  scale_fill_manual(values = hab_cols, drop = FALSE) +
  scale_color_manual(values = hab_cols, drop = FALSE) +
  labs(
    title = "Overall dRep MAG abundance of K08356-bearing host genomes",
    subtitle = "MAG TPM is sample-normalized; branch panels use copy-number adjusted MAG TPM; tests use sample-level log10(TPM + 1)",
    x = "Habitat",
    y = "Host MAG abundance, log10(TPM + 1)"
  ) +
  theme_mag

ggsave(file.path(fig_out, "overall_derepMAG_main_figure.png"), main_fig, width = 8.2, height = 8.2, dpi = 450)
ggsave(file.path(fig_out, "overall_derepMAG_main_figure.pdf"), main_fig, width = 8.2, height = 8.2)

total_label <- kw_labels %>%
  filter(Feature == "K08356-bearing host MAG total TPM") %>%
  pull(label)

total_fig <- ggplot(total %>% mutate(Habitat = factor(Habitat, levels = habitat_levels)),
                    aes(Habitat, logTPM, fill = Habitat)) +
  geom_boxplot(width = 0.58, alpha = 0.78, color = "grey20", outlier.shape = NA, linewidth = 0.45) +
  geom_jitter(aes(color = Habitat), width = 0.11, size = 2.1, alpha = 0.85, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.0, fill = "white", color = "black") +
  annotate("text", x = 2.5, y = Inf, label = total_label, vjust = 1.45, family = "Times", size = 4.2) +
  scale_fill_manual(values = hab_cols, drop = FALSE) +
  scale_color_manual(values = hab_cols, drop = FALSE) +
  labs(
    title = "Overall dRep MAG abundance of K08356-bearing host genomes",
    subtitle = "Sample-normalized host MAG TPM; Kruskal-Wallis test on log10(TPM + 1)",
    x = "Habitat",
    y = "Host MAG abundance, log10(TPM + 1)"
  ) +
  theme_mag +
  theme(legend.position = "none")

ggsave(file.path(fig_out, "overall_derepMAG_total_boxplot_square.png"), total_fig, width = 6.2, height = 6.2, dpi = 450)
ggsave(file.path(fig_out, "overall_derepMAG_total_boxplot_square.pdf"), total_fig, width = 6.2, height = 6.2)

branch_fig <- ggplot(
  branch %>%
    mutate(
      Habitat = factor(Habitat, levels = habitat_levels),
      Feature = factor(Feature, levels = feature_levels[-1])
    ),
  aes(Habitat, logTPM, fill = Habitat)
) +
  geom_boxplot(width = 0.58, alpha = 0.78, color = "grey20", outlier.shape = NA, linewidth = 0.42) +
  geom_jitter(aes(color = Habitat), width = 0.11, size = 1.85, alpha = 0.82, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white", color = "black") +
  geom_text(
    data = kw_labels %>% filter(Feature != "K08356-bearing host MAG total TPM"),
    aes(x = x, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.35,
    size = 3.4,
    family = "Times"
  ) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 2, labeller = as_labeller(feature_labels)) +
  scale_fill_manual(values = hab_cols, drop = FALSE) +
  scale_color_manual(values = hab_cols, drop = FALSE) +
  labs(
    title = "K08356 branch-resolved host MAG abundance",
    subtitle = "Copy-number adjusted MAG TPM; Kruskal-Wallis tests on sample-level log10(TPM + 1)",
    x = "Habitat",
    y = "Host MAG abundance, log10(TPM + 1)"
  ) +
  theme_mag

ggsave(file.path(fig_out, "overall_derepMAG_branch_boxplot_square.png"), branch_fig, width = 7.8, height = 7.8, dpi = 450)
ggsave(file.path(fig_out, "overall_derepMAG_branch_boxplot_square.pdf"), branch_fig, width = 7.8, height = 7.8)

svg_reference_sample_order <- c(
  "SY365BB-0-4", "SY365BB-4-8", "SY365BB-8-12",
  "SY366YB-0-4", "SY366YB-4-8", "SY366YB-8-12",
  "SY366YW-0-4", "SY366YW-4-8", "SY366YW-8-12",
  "SY368YW-0-4", "SY368YW-4-8", "SY368YW-8-12",
  "SY456YB-0-4", "SY456YB-4-8", "SY456YB-8-12",
  "SY457BB-0-4", "SY457BB-4-8", "SY457BB-8-12",
  "SY459WG-0-4", "SY459WG-4-8", "SY459WG-8-12",
  "S4_12-15", "S4_9-12",
  "SQ_58_0-4", "SQ_58_4-8", "SQ_58_8-12",
  "SQ_81_0-4", "SQ_81_4-8", "SQ_81_8-12",
  "S2_0-3", "S2_3-6", "S2_12-15",
  "S1_0-3", "S1_6-9", "S1_9-12",
  "C2_0-6", "C2_6-12", "C2_12-18",
  "C1_0-6", "C1_6-12", "C1_12-18",
  "C3_0-6", "C3_6-12", "C3_12-18",
  "S13_0-2", "S14_4-6", "S15_8-10", "ES_2_0-6",
  "S3_0-3", "S3_6-9", "S3_9-12",
  "R2111_N300_0-10", "R2111_N500_0-10",
  "R2111_S300_0-10", "R2111_S500_0-10",
  "NS_0-6"
)

missing_from_reference <- setdiff(unique(branch$Sample), svg_reference_sample_order)
missing_from_data <- setdiff(svg_reference_sample_order, unique(branch$Sample))
if (length(missing_from_reference) > 0 || length(missing_from_data) > 0) {
  warning(
    "Sample order mismatch. Missing from SVG reference order: ",
    paste(missing_from_reference, collapse = ", "),
    "; missing from MAG data: ",
    paste(missing_from_data, collapse = ", ")
  )
}

stack_dat <- branch %>%
  mutate(
    Habitat = factor(Habitat, levels = habitat_levels),
    Feature = factor(Feature, levels = feature_levels[-1]),
    Sample_ordered = factor(Sample, levels = svg_reference_sample_order)
  )

branch_cols <- c(
  "aioA-like-associated" = "#009E73",
  "unknown AioA-like / uncertain DMSOR" = "#CC79A7",
  "IdrA-associated" = "#5B8DB8",
  "canonical aioA-associated" = "#D55E00"
)

stack_fig <- ggplot(stack_dat, aes(Sample_ordered, TPM, fill = Feature)) +
  geom_col(width = 0.82, color = "grey35", linewidth = 0.08) +
  facet_grid(. ~ Habitat, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = branch_cols, drop = FALSE, labels = feature_labels[names(branch_cols)]) +
  labs(
    title = "Overall dRep MAG branch composition across samples",
    subtitle = "Stacked bars show copy-number adjusted host MAG TPM for each K08356 branch",
    x = "Sample order matched to class-level grouped figure",
    y = "Host MAG TPM",
    fill = "K08356 branch"
  ) +
  theme_mag +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 6.4),
    strip.text = element_text(size = 11),
    legend.position = "bottom"
  )

ggsave(file.path(fig_out, "overall_derepMAG_branch_sample_order_stacked_bar.png"), stack_fig, width = 10.2, height = 8.2, dpi = 450)
ggsave(file.path(fig_out, "overall_derepMAG_branch_sample_order_stacked_bar.pdf"), stack_fig, width = 10.2, height = 8.2)

pairwise <- read_csv(
  file.path(mag_out, "K08356_branch_host_MAG_pairwise_Wilcoxon_BH.csv"),
  show_col_types = FALSE
)

readme <- c(
  "# Overall dRepMAG figure notes",
  "",
  "This folder contains the finalized overall dRepMAG figures for K08356-bearing host MAG abundance.",
  "",
  "Main statistical unit: sample.",
  "Input abundance: `merged_tpm_matrix.txt` MAG TPM. Each sample column is already TPM-normalized, with column sums close to 1,000,000.",
  "Biological meaning: MAG TPM describes the abundance of host derepMAGs carrying K08356-related genes; it is not K08356 gene/contig TPM.",
  "Branch panels: if one derepMAG contains multiple K08356 copies assigned to different branches, that MAG TPM is split equally among its K08356 copy rows before branch-level summation.",
  "Significance: Kruskal-Wallis tests are run on sample-level log10(TPM + 1). Pairwise Wilcoxon tests use BH adjustment and are retained in the CSV/Excel tables.",
  "",
  "Key results for wording:",
  "- Total K08356-bearing host MAG TPM differs significantly among habitats: Kruskal-Wallis p = 1.32e-6.",
  "- aioA-like-associated host MAG TPM differs significantly among habitats: Kruskal-Wallis p = 2.20e-7; IS is significantly higher than AS, ES, and NS in pairwise BH-adjusted tests.",
  "- unknown AioA-like / uncertain DMSOR host MAG TPM differs significantly among habitats: Kruskal-Wallis p = 4.01e-4; IS is higher than AS and ES, but IS vs NS is not significant after BH adjustment.",
  "- IdrA-associated host MAG TPM is not significant among habitats after copy-number adjustment: Kruskal-Wallis p = 0.121.",
  "- canonical aioA-associated host MAG TPM is not significant among habitats: Kruskal-Wallis p = 0.877.",
  "",
  "Recommended placement:",
  "- `overall_derepMAG_total_boxplot_square` and `overall_derepMAG_branch_boxplot_square` are the cleanest main Figure A/B files.",
  "- `overall_derepMAG_main_figure` is a one-file combined overview, but it has one blank facet slot because it contains five panels.",
  "- `overall_derepMAG_branch_sample_order_stacked_bar` is better used as a supplementary/sample-order distribution figure; its x-axis order matches `/Users/catherine/Downloads/metaphlan 作图/class34.csv 作图/最终版图片/class_level_grouped_final_v2.svg` and the previous grouped-analysis order.",
  "- Grouped-derep results should be added later as an appendix/sensitivity analysis, not as a replacement for the overall dRepMAG abundance figure."
)
writeLines(readme, file.path(fig_out, "README_overall_derepMAG_figure.md"))

write_csv(kw_labels %>% select(Feature, p_value, label), file.path(fig_out, "overall_derepMAG_Kruskal_labels_used.csv"))
write_csv(pairwise, file.path(fig_out, "overall_derepMAG_branch_pairwise_Wilcoxon_BH.csv"))

message("Wrote overall dRepMAG final figures to: ", fig_out)

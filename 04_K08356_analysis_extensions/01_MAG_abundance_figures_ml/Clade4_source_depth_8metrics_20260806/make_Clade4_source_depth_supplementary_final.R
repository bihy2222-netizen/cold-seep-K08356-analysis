suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "Clade4_source_depth_8metrics_20260806")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clade4 <- readLines(file.path(
  base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected", "Clade4_carrier_17_MAGs.txt"
))
depth_levels <- c("0-4", "4-8", "8-12")
depth_colors <- c("0-4" = "#8DD3C7", "4-8" = "#2CA25F", "8-12" = "#225EA8")
status_colors <- c("Clade 4" = "#7B4FA3", "non-Clade 4" = "#B8B8B8")

if (length(clade4) != 17 || anyDuplicated(clade4)) stop("Expected 17 unique Clade 4 MAGs.")

quality_report <- read_tsv(
  "/Users/catherine/Downloads/Contig 热图柱状图 new/quality_report.tsv",
  show_col_types = FALSE
) %>%
  filter(Name %in% clade4) %>%
  transmute(
    MAG = Name,
    source_depth = str_match(Name, "-(0-4|4-8|8-12)_bin")[, 2],
    Genome_size_Mbp = Genome_Size / 1e6,
    Coding_Density_fraction = Coding_Density,
    Completeness_pct = Completeness,
    Contamination_pct = Contamination,
    Contig_N50_kbp = Contig_N50 / 1e3,
    Total_CDS_k = Total_Coding_Sequences / 1e3
  )

high_precision_gc <- read_excel(
  "/Users/catherine/Downloads/MAGs 建树/bin_tax.xlsx",
  sheet = "gtdbtk.ar53.bac120.summary"
) %>%
  filter(user_genome %in% clade4) %>%
  transmute(MAG = user_genome, GC_Content_pct = 100 * GC)

quality <- quality_report %>%
  left_join(high_precision_gc, by = "MAG") %>%
  mutate(source_depth = factor(source_depth, levels = depth_levels))

if (nrow(quality) != 17 || n_distinct(quality$MAG) != 17 || anyNA(quality)) {
  stop("All 17 MAGs require complete quality and high-precision GC data.")
}

assignments <- read.csv(
  file.path(base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected",
            "Rhodobacteraceae_MAG_group_assignments_corrected.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
) %>%
  filter(source_habitat == "IS") %>%
  mutate(
    source_depth = str_match(MAG, "-(0-4|4-8|8-12)_bin")[, 2],
    carrier_status = if_else(comparison_group == "Clade 4 carrier", "Clade 4", "non-Clade 4"),
    source_depth = factor(source_depth, levels = depth_levels),
    carrier_status = factor(carrier_status, levels = c("Clade 4", "non-Clade 4"))
  )

depth_counts <- assignments %>% count(source_depth, carrier_status, name = "MAG_count")
depth_table <- xtabs(MAG_count ~ source_depth + carrier_status, data = depth_counts)
depth_fisher <- fisher.test(depth_table)

expected_depth_table <- matrix(c(4, 7, 2, 6, 11, 10), nrow = 3, byrow = TRUE,
                               dimnames = list(depth_levels, c("Clade 4", "non-Clade 4")))
if (!isTRUE(all.equal(as.numeric(depth_table), as.numeric(expected_depth_table)))) {
  stop("IS Rhodobacteraceae depth contingency table does not match required 4/7, 2/6, 11/10 counts.")
}

metric_info <- tibble(
  metric = c(
    "GC_Content_pct", "Genome_size_Mbp", "Coding_Density_fraction", "Completeness_pct",
    "Contamination_pct", "Contig_N50_kbp", "Total_CDS_k"
  ),
  panel = LETTERS[2:8],
  label = c(
    "GC content (%)", "Genome size (Mb)", "Coding density (fraction)", "Completeness (%)",
    "Contamination (%)", "Contig N50 (kb)", "Coding sequences (×10³)"
  )
)

plot_df <- quality %>%
  pivot_longer(cols = all_of(metric_info$metric), names_to = "metric", values_to = "value") %>%
  left_join(metric_info, by = "metric") %>%
  mutate(
    facet_label = factor(
      paste0(panel, "\n", label),
      levels = paste0(metric_info$panel, "\n", metric_info$label)
    )
  )

p_to_star <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

pair_rows <- list()
for (current_metric in metric_info$metric) {
  current_data <- plot_df %>% filter(metric == current_metric)
  for (depth_pair in combn(depth_levels, 2, simplify = FALSE)) {
    left_values <- current_data %>% filter(source_depth == depth_pair[1]) %>% pull(value)
    right_values <- current_data %>% filter(source_depth == depth_pair[2]) %>% pull(value)
    test_p <- if (length(left_values) < 2 || length(right_values) < 2) {
      NA_real_
    } else {
      suppressWarnings(wilcox.test(left_values, right_values, exact = FALSE)$p.value)
    }
    pair_rows[[length(pair_rows) + 1]] <- tibble(
      metric = current_metric, group1 = depth_pair[1], group2 = depth_pair[2],
      group1_n = length(left_values), group2_n = length(right_values), wilcox_p = test_p
    )
  }
}

pairwise <- bind_rows(pair_rows) %>%
  group_by(metric) %>%
  mutate(BH_FDR = p.adjust(wilcox_p, method = "BH"), star = p_to_star(BH_FDR)) %>%
  ungroup() %>%
  left_join(metric_info, by = "metric")

global_tests <- plot_df %>%
  group_by(metric, panel, label, facet_label) %>%
  group_modify(~ {
    result <- kruskal.test(value ~ source_depth, data = .x)
    tibble(KW_statistic = unname(result$statistic), KW_p = result$p.value)
  }) %>%
  ungroup() %>%
  mutate(
    panel_note = paste0(
      "Kruskal–Wallis P=", formatC(KW_p, digits = 3, format = "f"),
      "\nAll pairwise BH-FDR >0.05"
    )
  )

if (abs(global_tests$KW_p[global_tests$metric == "GC_Content_pct"] - 0.6472639) > 1e-5) {
  stop("GC P-value provenance QC failed; expected high-precision bin_tax result P=0.647264.")
}
if (any(pairwise$BH_FDR < 0.05, na.rm = TRUE)) {
  stop("Unexpected significant pairwise result; inspect before plotting.")
}

summary_table <- plot_df %>%
  group_by(metric, panel, label, source_depth) %>%
  summarise(
    n = n(), mean = mean(value), sd = sd(value), median = median(value),
    q1 = quantile(value, 0.25), q3 = quantile(value, 0.75),
    min = min(value), max = max(value), .groups = "drop"
  )

plot_ranges <- plot_df %>%
  group_by(metric, facet_label) %>%
  summarise(ymin = min(value), ymax = max(value), .groups = "drop") %>%
  mutate(yrange = pmax(ymax - ymin, abs(ymax) * 0.08, 0.1), note_y = ymax + 0.16 * yrange)
panel_notes <- global_tests %>% left_join(plot_ranges, by = c("metric", "facet_label"))

outlier_label <- quality %>%
  filter(Contamination_pct > 20) %>%
  transmute(
    MAG, source_depth, metric = "Contamination_pct", value = Contamination_pct,
    facet_label = factor("F\nContamination (%)", levels = levels(plot_df$facet_label)),
    label = MAG
  )
if (nrow(outlier_label) != 1 || outlier_label$MAG != "SY365BB-8-12_bin8") {
  stop("Expected one >20% contamination outlier: SY365BB-8-12_bin8.")
}

p_counts <- ggplot(depth_counts, aes(source_depth, MAG_count, fill = carrier_status)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.66,
           color = "#303030", linewidth = 0.40) +
  geom_text(aes(label = MAG_count), position = position_dodge(width = 0.74),
            vjust = -0.28, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = status_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 12, 2), expand = c(0, 0)) +
  labs(
    x = "Source depth (cm)", y = "Number of IS Rhodobacteraceae MAGs",
    title = "A\nSource-depth recovery counts",
    subtitle = paste0("Clade 4 vs non-Clade 4: Fisher P=",
                      formatC(depth_fisher$p.value, digits = 3, format = "f"), " (ns)")
  ) +
  theme_bw(base_size = 9.5, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", size = 11.2, hjust = 0.5),
    plot.subtitle = element_text(size = 8.0, hjust = 0.5),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, face = "bold"),
    legend.position = "bottom"
  )

p_metrics <- ggplot(plot_df, aes(source_depth, value)) +
  geom_boxplot(aes(color = source_depth), width = 0.30, outlier.shape = NA,
               fill = "white", linewidth = 0.55, show.legend = FALSE) +
  geom_point(aes(fill = source_depth), shape = 21, color = "#303030", stroke = 0.25,
             size = 2.3, alpha = 0.92,
             position = position_jitter(width = 0.10, height = 0, seed = 20260806),
             show.legend = FALSE) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.8,
               fill = "black", color = "black") +
  geom_text(
    data = panel_notes,
    aes(x = 1, y = note_y, label = panel_note),
    inherit.aes = FALSE, hjust = 0, size = 2.55, lineheight = 0.95, color = "#333333"
  ) +
  geom_text(
    data = outlier_label,
    aes(x = source_depth, y = value, label = label),
    inherit.aes = FALSE, nudge_x = -0.18, nudge_y = 0.70,
    hjust = 1, size = 2.35, color = "#8B1A1A"
  ) +
  facet_wrap(~facet_label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = depth_colors) +
  scale_fill_manual(values = depth_colors) +
  scale_x_discrete(labels = c("0-4" = "0–4", "4-8" = "4–8", "8-12" = "8–12")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.28))) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9.5, base_family = "Times New Roman") +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 10.5, lineheight = 0.95),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, face = "bold", size = 8.5),
    axis.text.y = element_text(size = 8.2)
  )

metric_panels <- lapply(levels(plot_df$facet_label), function(current_facet) {
  current_data <- plot_df %>% filter(facet_label == current_facet)
  current_note <- panel_notes %>% filter(facet_label == current_facet)
  current_outlier <- outlier_label %>% filter(facet_label == current_facet)

  ggplot(current_data, aes(source_depth, value)) +
    geom_boxplot(aes(color = source_depth), width = 0.30, outlier.shape = NA,
                 fill = "white", linewidth = 0.55, show.legend = FALSE) +
    geom_point(aes(fill = source_depth), shape = 21, color = "#303030", stroke = 0.25,
               size = 2.3, alpha = 0.92,
               position = position_jitter(width = 0.10, height = 0, seed = 20260806),
               show.legend = FALSE) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2.8,
                 fill = "black", color = "black") +
    geom_text(
      data = current_note,
      aes(x = 1, y = note_y, label = panel_note),
      inherit.aes = FALSE, hjust = 0, size = 2.55, lineheight = 0.95, color = "#333333"
    ) +
    geom_text(
      data = current_outlier,
      aes(x = source_depth, y = value, label = label),
      inherit.aes = FALSE, nudge_x = -0.18, nudge_y = 0.70,
      hjust = 1, size = 2.35, color = "#8B1A1A"
    ) +
    facet_wrap(~facet_label, scales = "free_y") +
    scale_color_manual(values = depth_colors) +
    scale_fill_manual(values = depth_colors) +
    scale_x_discrete(labels = c("0-4" = "0–4", "4-8" = "4–8", "8-12" = "8–12")) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.28))) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9.5, base_family = "Times New Roman") +
    theme(
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", size = 10.5, lineheight = 0.95),
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 40, hjust = 1, face = "bold", size = 8.5),
      axis.text.y = element_text(size = 8.2)
    )
})

figure <- wrap_plots(c(list(p_counts), metric_panels), ncol = 4) +
  plot_annotation(
    title = "Source-depth recovery and genomic features of 17 Clade 4 / DIRM-like carriers",
    subtitle = "Raw MAG-level points, boxplots and group medians show no detectable depth-associated genomic-feature differences",
    caption = paste0(
      "Black diamonds denote group medians. Panels B–H use 17 MAGs distributed as n=4, 2 and 11 across 0–4, 4–8 and 8–12 cm. ",
      "All genomic-feature pairwise Wilcoxon tests had BH-FDR >0.05; null results should be interpreted cautiously because the 4–8 cm group contains only two MAGs. ",
      "GC uses the higher-precision assembly GC field from bin_tax.xlsx (Kruskal–Wallis P=0.647); the two-decimal quality_report.tsv GC field yields P=0.700 because rounding changes ranks/ties. ",
      "SY365BB-8-12_bin8 (contamination 23.4%) is retained and labeled. Panel A includes the IS Rhodobacteraceae background and shows no source-depth association with Clade 4 status. ",
      "ANI is withheld because a complete 17-genome all-vs-all fastANI matrix (136 non-self pairs; 16 comparisons per MAG) is not currently available."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15.5),
      plot.subtitle = element_text(size = 9.5, color = "#404040"),
      plot.caption = element_text(size = 7.6, color = "#444444", hjust = 0)
    )
  )

save_plot <- function(extension, device, dpi = 450) {
  args <- list(
    filename = file.path(out_dir, paste0("Clade4_source_depth_genomic_features_supplementary_final.", extension)),
    plot = figure, width = 15.5, height = 10.5, units = "in", bg = "white", device = device
  )
  if (extension == "png") args$dpi <- dpi
  do.call(ggsave, args)
}
save_plot("pdf", cairo_pdf)
save_plot("svg", grDevices::svg)
save_plot("png", "png")

write.csv(quality, file.path(out_dir, "Clade4_17_MAG_genomic_features_final.csv"), row.names = FALSE)
write.csv(depth_counts, file.path(out_dir, "Clade4_vs_nonClade4_IS_source_depth_counts.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(depth_table),
          file.path(out_dir, "Clade4_vs_nonClade4_IS_source_depth_contingency.csv"))
write.csv(summary_table, file.path(out_dir, "Clade4_source_depth_genomic_feature_summary_final.csv"), row.names = FALSE)
write.csv(global_tests, file.path(out_dir, "Clade4_source_depth_genomic_feature_KW_final.csv"), row.names = FALSE)
write.csv(pairwise, file.path(out_dir, "Clade4_source_depth_genomic_feature_pairwise_BH_FDR_final.csv"), row.names = FALSE)

qc <- c(
  "Clade 4 source-depth genomic-feature supplementary figure final QC",
  "Status: PASSED",
  paste0("Clade 4 MAGs: ", nrow(quality), " (required 17)"),
  "Source-depth counts Clade 4: 4 / 2 / 11",
  "Source-depth counts non-Clade 4 IS Rhodobacteraceae: 7 / 6 / 10",
  paste0("Depth contingency Fisher P: ", signif(depth_fisher$p.value, 6)),
  paste0("GC Kruskal-Wallis P (high-precision bin_tax field): ",
         signif(global_tests$KW_p[global_tests$metric == "GC_Content_pct"], 6)),
  paste0("Significant pairwise genomic-feature BH-FDR tests: ",
         sum(pairwise$BH_FDR < 0.05, na.rm = TRUE)),
  "Plot geometry: raw points + boxplots + black group-median diamonds; no violins.",
  "Contamination outlier SY365BB-8-12_bin8 is retained and labeled.",
  "Coding density is labeled as a fraction.",
  "ANI panel withheld: complete 17x17 all-vs-all fastANI results are unavailable.",
  "No deep-source enrichment, depth adaptation or redox-interface localization claim is made."
)
writeLines(qc, file.path(out_dir, "Clade4_source_depth_genomic_features_supplementary_final_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

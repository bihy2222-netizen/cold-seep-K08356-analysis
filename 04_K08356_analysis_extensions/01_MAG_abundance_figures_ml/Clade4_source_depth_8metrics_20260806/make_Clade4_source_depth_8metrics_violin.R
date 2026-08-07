suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "Clade4_source_depth_8metrics_20260806")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

quality_file <- "/Users/catherine/Downloads/Contig 热图柱状图 new/quality_report.tsv"
ani_file <- "/Users/catherine/Downloads/Contig 热图柱状图 new/outputs/mag_ani_drep_tables/Ndb.csv"
clade4 <- readLines(file.path(
  base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected", "Clade4_carrier_17_MAGs.txt"
))
depth_levels <- c("0-4", "4-8", "8-12")
depth_colors <- c("0-4" = "#8DD3C7", "4-8" = "#2CA25F", "8-12" = "#225EA8")

if (length(clade4) != 17 || anyDuplicated(clade4)) stop("Expected 17 unique Clade 4 MAGs.")

quality <- read_tsv(quality_file, show_col_types = FALSE) %>%
  filter(Name %in% clade4) %>%
  transmute(
    MAG = Name,
    source_depth = str_match(Name, "-(0-4|4-8|8-12)_bin")[, 2],
    Genome_size_Mbp = Genome_Size / 1e6,
    Coding_Density = Coding_Density,
    GC_Content_pct = 100 * GC_Content,
    Completeness_pct = Completeness,
    Contamination_pct = Contamination,
    Contig_N50_kbp = Contig_N50 / 1e3,
    Total_CDS_k = Total_Coding_Sequences / 1e3
  ) %>%
  mutate(source_depth = factor(source_depth, levels = depth_levels))

if (nrow(quality) != 17 || n_distinct(quality$MAG) != 17 || anyNA(quality$source_depth)) {
  stop("All 17 Clade 4 MAGs must match the quality table and source depths.")
}

ani_per_mag <- read_csv(ani_file, show_col_types = FALSE) %>%
  mutate(
    reference = str_remove(basename(reference), "\\.(fa|fna|fasta)$"),
    querry = str_remove(basename(querry), "\\.(fa|fna|fasta)$")
  ) %>%
  filter(reference %in% clade4, reference != querry, is.finite(ani), ani > 0) %>%
  group_by(reference) %>%
  summarise(
    ANI_median_pct = 100 * median(ani),
    ANI_mean_pct = 100 * mean(ani),
    ANI_neighbor_comparisons = n(),
    .groups = "drop"
  ) %>%
  rename(MAG = reference) %>%
  left_join(quality %>% select(MAG, source_depth), by = "MAG")

metric_info <- tibble(
  metric = c(
    "Genome_size_Mbp", "Coding_Density", "GC_Content_pct", "Completeness_pct",
    "Contamination_pct", "Contig_N50_kbp", "Total_CDS_k", "ANI_median_pct"
  ),
  panel = LETTERS[1:8],
  label = c(
    "Genome size (Mb)", "Coding density", "GC content (%)", "Completeness (%)",
    "Contamination (%)", "Contig N50 (kb)", "Coding sequences (×10³)",
    "Median dRep-neighbor ANI (%)"
  )
)

quality_long <- quality %>%
  pivot_longer(
    cols = all_of(metric_info$metric[1:7]),
    names_to = "metric", values_to = "value"
  )
ani_long <- ani_per_mag %>%
  transmute(MAG, source_depth, metric = "ANI_median_pct", value = ANI_median_pct)

plot_df <- bind_rows(quality_long, ani_long) %>%
  left_join(metric_info, by = "metric") %>%
  mutate(
    source_depth = factor(source_depth, levels = depth_levels),
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
  metric_data <- plot_df %>% filter(metric == current_metric)
  for (depth_pair in combn(depth_levels, 2, simplify = FALSE)) {
    left_values <- metric_data %>% filter(source_depth == depth_pair[1]) %>% pull(value)
    right_values <- metric_data %>% filter(source_depth == depth_pair[2]) %>% pull(value)
    test_p <- if (length(left_values) < 2 || length(right_values) < 2 ||
                  length(unique(c(left_values, right_values))) <= 1) {
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
    valid_groups <- .x %>% count(source_depth) %>% filter(n >= 2)
    if (nrow(valid_groups) < 2) return(tibble(KW_statistic = NA_real_, KW_p = NA_real_))
    result <- kruskal.test(value ~ source_depth, data = .x)
    tibble(KW_statistic = unname(result$statistic), KW_p = result$p.value)
  }) %>%
  ungroup() %>%
  mutate(global_label = if_else(
    is.na(KW_p), "Global test: insufficient n",
    paste0("KW P=", formatC(KW_p, digits = 3, format = "f"), " (", p_to_star(KW_p), ")")
  ))

summary_table <- plot_df %>%
  group_by(metric, panel, label, source_depth) %>%
  summarise(
    n = n(), mean = mean(value), sd = sd(value), median = median(value),
    q1 = quantile(value, 0.25), q3 = quantile(value, 0.75), .groups = "drop"
  )

plot_ranges <- plot_df %>%
  group_by(metric, facet_label) %>%
  summarise(ymin = min(value), ymax = max(value), .groups = "drop") %>%
  mutate(yrange = pmax(ymax - ymin, abs(ymax) * 0.08, 0.1), label_y = ymax + 0.13 * yrange)

significant_brackets <- pairwise %>%
  filter(!is.na(BH_FDR), BH_FDR < 0.05) %>%
  left_join(plot_ranges, by = "metric") %>%
  group_by(metric) %>%
  arrange(BH_FDR, .by_group = TRUE) %>%
  mutate(
    bracket_id = row_number(),
    xmin = match(group1, depth_levels), xmax = match(group2, depth_levels),
    y = ymax + (0.18 + 0.14 * (bracket_id - 1)) * yrange,
    yend = y - 0.04 * yrange
  ) %>%
  ungroup()

global_labels <- global_tests %>%
  left_join(plot_ranges, by = c("metric", "facet_label"))

p <- ggplot(plot_df, aes(source_depth, value, fill = source_depth)) +
  geom_violin(width = 0.88, trim = TRUE, scale = "width", color = "#303030",
              linewidth = 0.38, alpha = 0.38, na.rm = TRUE) +
  geom_boxplot(width = 0.20, outlier.shape = NA, fill = "white", color = "#303030",
               linewidth = 0.38, alpha = 0.86, na.rm = TRUE) +
  geom_point(aes(color = source_depth), shape = 21, stroke = 0.15, size = 1.5,
             alpha = 0.75, position = position_jitter(width = 0.10, height = 0, seed = 20260806),
             show.legend = FALSE, na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.3,
               fill = "black", color = "black", na.rm = TRUE) +
  geom_text(
    data = global_labels,
    aes(x = 1, y = label_y, label = global_label),
    inherit.aes = FALSE, hjust = 0, size = 2.9, color = "#333333"
  ) +
  geom_segment(
    data = significant_brackets,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE, color = "black", linewidth = 0.42
  ) +
  geom_segment(
    data = significant_brackets,
    aes(x = xmin, xend = xmin, y = y, yend = yend),
    inherit.aes = FALSE, color = "black", linewidth = 0.42
  ) +
  geom_segment(
    data = significant_brackets,
    aes(x = xmax, xend = xmax, y = y, yend = yend),
    inherit.aes = FALSE, color = "black", linewidth = 0.42
  ) +
  geom_text(
    data = significant_brackets,
    aes(x = (xmin + xmax) / 2, y = y + 0.03 * yrange, label = star),
    inherit.aes = FALSE, fontface = "bold", size = 4.0
  ) +
  facet_wrap(~facet_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = depth_colors, name = "Source depth (cm)") +
  scale_color_manual(values = depth_colors) +
  scale_x_discrete(labels = c("0-4" = "0–4", "4-8" = "4–8", "8-12" = "8–12")) +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.24))) +
  labs(
    x = NULL, y = NULL,
    title = "Genomic features and dRep-neighbor ANI of Clade 4 carriers across source depths",
    subtitle = "Stars indicate BH-FDR-adjusted pairwise Wilcoxon tests (* <0.05, ** <0.01, *** <0.001); no pairwise comparison passed FDR, so panels are marked ns"
  ) +
  theme_bw(base_size = 10.5, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", size = 15.5),
    plot.subtitle = element_text(size = 9.1),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 11.2, lineheight = 0.95),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, face = "bold", size = 9),
    axis.text.y = element_text(size = 8.5), legend.position = "right",
    plot.caption = element_text(size = 7.8, color = "#444444", hjust = 0)
  ) +
  labs(caption = paste0(
    "Each point is one MAG for panels A–G (n=17: 4, 2 and 11 MAGs at 0–4, 4–8 and 8–12 cm). ",
    "Panel H uses one median ANI value per MAG across available non-self dRep Ndb comparisons (n=9: 2, 1 and 6); eight MAGs lacked non-self ANI entries. ",
    "The 4–8 cm distributions are descriptive because n=2 (and n=1 for ANI). Source depth is parsed from MAG names and does not itself establish depth adaptation."
  ))

save_plot <- function(extension, device, dpi = 450) {
  args <- list(
    filename = file.path(out_dir, paste0("Clade4_source_depth_8metrics_violin.", extension)),
    plot = p, width = 15.2, height = 10.4, units = "in", bg = "white", device = device
  )
  if (extension == "png") args$dpi <- dpi
  do.call(ggsave, args)
}
save_plot("pdf", cairo_pdf)
save_plot("svg", grDevices::svg)
save_plot("png", "png")

write.csv(quality, file.path(out_dir, "Clade4_17_MAG_quality_metrics.csv"), row.names = FALSE)
write.csv(ani_per_mag, file.path(out_dir, "Clade4_MAG_median_dRep_neighbor_ANI.csv"), row.names = FALSE)
write.csv(plot_df, file.path(out_dir, "Clade4_source_depth_8metrics_plot_data.csv"), row.names = FALSE)
write.csv(summary_table, file.path(out_dir, "Clade4_source_depth_8metrics_summary.csv"), row.names = FALSE)
write.csv(global_tests, file.path(out_dir, "Clade4_source_depth_8metrics_global_KW.csv"), row.names = FALSE)
write.csv(pairwise, file.path(out_dir, "Clade4_source_depth_8metrics_pairwise_Wilcoxon_BH_FDR.csv"), row.names = FALSE)

qc <- c(
  "Clade 4 source-depth eight-metric violin figure QC",
  "Status: PASSED",
  paste0("Quality-metric MAGs: ", nrow(quality), " (required 17)"),
  paste0("Per-MAG ANI summaries: ", nrow(ani_per_mag), " (9 available; 8 without non-self Ndb entries)"),
  paste0("Metrics: ", n_distinct(plot_df$metric), " (required 8)"),
  paste0("Significant pairwise comparisons at BH-FDR < 0.05: ", sum(pairwise$BH_FDR < 0.05, na.rm = TRUE)),
  paste0("Significance brackets/stars drawn: ", nrow(significant_brackets)),
  "No significance star is fabricated; ns is shown when the global comparison is not significant.",
  "ANI is summarized once per MAG to reduce pairwise pseudoreplication.",
  "The 4-8 cm group is descriptive because of small n."
)
writeLines(qc, file.path(out_dir, "Clade4_source_depth_8metrics_violin_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

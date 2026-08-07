suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "Clade4_source_depth_GC_genome_20260806")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clade4 <- readLines(file.path(
  base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected", "Clade4_carrier_17_MAGs.txt"
))
quality_file <- "/Users/catherine/Downloads/MAGs 建树/bin_tax.xlsx"

if (length(clade4) != 17 || anyDuplicated(clade4)) stop("Clade 4 list must contain 17 unique MAGs.")

quality <- read_excel(quality_file, sheet = "gtdbtk.ar53.bac120.summary") %>%
  filter(user_genome %in% clade4) %>%
  transmute(
    MAG = user_genome,
    source_depth = str_match(user_genome, "-(0-4|4-8|8-12)_bin")[, 2],
    GC_percent = 100 * GC,
    genome_size_Mb = size / 1e6,
    completeness,
    contamination,
    N50
  ) %>%
  mutate(
    source_depth = factor(source_depth, levels = c("0-4", "4-8", "8-12")),
    MAG = factor(MAG, levels = clade4)
  )

if (nrow(quality) != 17 || n_distinct(quality$MAG) != 17 || anyNA(quality$source_depth) ||
    anyNA(quality$GC_percent) || anyNA(quality$genome_size_Mb)) {
  stop("All 17 Clade 4 MAGs must have valid source depth, GC and genome-size data.")
}

depth_counts <- quality %>%
  count(source_depth, name = "MAG_count") %>%
  mutate(
    percentage = 100 * MAG_count / sum(MAG_count),
    count_label = paste0("n=", MAG_count, "\n", formatC(percentage, digits = 1, format = "f"), "%")
  )

kw_gc <- kruskal.test(GC_percent ~ source_depth, data = quality)
kw_size <- kruskal.test(genome_size_Mb ~ source_depth, data = quality)
kw_completeness <- kruskal.test(completeness ~ source_depth, data = quality)

pairwise_gc <- pairwise.wilcox.test(
  quality$GC_percent, quality$source_depth, p.adjust.method = "BH", exact = FALSE
)
pairwise_size <- pairwise.wilcox.test(
  quality$genome_size_Mb, quality$source_depth, p.adjust.method = "BH", exact = FALSE
)

summary_table <- quality %>%
  group_by(source_depth) %>%
  summarise(
    n = n(),
    GC_mean = mean(GC_percent), GC_median = median(GC_percent), GC_sd = sd(GC_percent),
    genome_size_mean_Mb = mean(genome_size_Mb),
    genome_size_median_Mb = median(genome_size_Mb), genome_size_sd_Mb = sd(genome_size_Mb),
    completeness_mean = mean(completeness), completeness_median = median(completeness),
    .groups = "drop"
  )

statistics <- tibble(
  metric = c("GC content (%)", "Genome size (Mb)", "Completeness (%)"),
  test = "Kruskal-Wallis",
  statistic = c(unname(kw_gc$statistic), unname(kw_size$statistic), unname(kw_completeness$statistic)),
  df = c(unname(kw_gc$parameter), unname(kw_size$parameter), unname(kw_completeness$parameter)),
  p_value = c(kw_gc$p.value, kw_size$p.value, kw_completeness$p.value)
)

write.csv(quality %>% mutate(MAG = as.character(MAG)),
          file.path(out_dir, "Clade4_17_MAG_source_depth_GC_genome_quality.csv"), row.names = FALSE)
write.csv(depth_counts, file.path(out_dir, "Clade4_source_depth_counts.csv"), row.names = FALSE)
write.csv(summary_table, file.path(out_dir, "Clade4_source_depth_GC_genome_summary.csv"), row.names = FALSE)
write.csv(statistics, file.path(out_dir, "Clade4_source_depth_GC_genome_statistics.csv"), row.names = FALSE)

depth_colors <- c("0-4" = "#8DD3C7", "4-8" = "#2CA25F", "8-12" = "#225EA8")
depth_fills <- c("0-4" = "#CDECE7", "4-8" = "#A8DDB5", "8-12" = "#9ECAE1")

p_depth <- ggplot(depth_counts, aes(source_depth, MAG_count, fill = source_depth)) +
  geom_col(width = 0.66, color = "#303030", linewidth = 0.45) +
  geom_text(aes(label = count_label), vjust = -0.25, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = depth_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 12, 2), expand = c(0, 0)) +
  labs(
    x = "Source depth of Clade 4 MAG (cm)", y = "Number of MAGs",
    title = "A  Source-depth distribution of the 17 Clade 4 carriers",
    subtitle = "Depth was parsed directly from each MAG name"
  ) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.8, color = "#444444"),
        axis.text.x = element_text(face = "bold"))

violin_layer <- function() {
  list(
    geom_violin(aes(fill = source_depth), width = 0.82, trim = FALSE,
                alpha = 0.35, color = "#404040", linewidth = 0.5),
    geom_boxplot(aes(fill = source_depth), width = 0.18, outlier.shape = NA,
                 alpha = 0.72, color = "#303030", linewidth = 0.45),
    geom_jitter(aes(color = source_depth), width = 0.10, height = 0,
                size = 2.5, alpha = 0.95, show.legend = FALSE),
    scale_fill_manual(values = depth_fills, guide = "none"),
    scale_color_manual(values = depth_colors, guide = "none")
  )
}

p_gc <- ggplot(quality, aes(source_depth, GC_percent)) +
  violin_layer() +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.1,
               fill = "white", color = "black") +
  labs(
    x = "Source depth (cm)", y = "Genome GC content (%)",
    title = "B  Genome GC content across source depths",
    subtitle = paste0("Kruskal–Wallis P=", formatC(kw_gc$p.value, digits = 3, format = "f"))
  ) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.8, color = "#444444"),
        axis.text.x = element_text(face = "bold"))

p_size <- ggplot(quality, aes(source_depth, genome_size_Mb)) +
  violin_layer() +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3.1,
               fill = "white", color = "black") +
  labs(
    x = "Source depth (cm)", y = "Genome size (Mb)",
    title = "C  Genome size across source depths",
    subtitle = paste0("Kruskal–Wallis P=", formatC(kw_size$p.value, digits = 3, format = "f"))
  ) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.8, color = "#444444"),
        axis.text.x = element_text(face = "bold"))

name_strip <- quality %>%
  arrange(source_depth, MAG) %>%
  mutate(
    MAG_label = str_replace(as.character(MAG), "^SY", "SY"),
    row_id = row_number()
  )

p_names <- ggplot(name_strip, aes(row_id, y = 1, fill = source_depth)) +
  geom_tile(width = 0.92, height = 0.60, color = "white", linewidth = 0.35) +
  scale_fill_manual(values = depth_fills, name = "Source depth (cm)") +
  scale_x_continuous(breaks = name_strip$row_id, labels = name_strip$MAG_label,
                     expand = expansion(add = 0.55)) +
  coord_cartesian(ylim = c(0.65, 1.35), clip = "off") +
  labs(x = NULL, y = NULL, title = "D  Individual Clade 4 MAGs grouped by source depth") +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.2, color = "black"),
        axis.ticks.x = element_blank(), plot.title = element_text(face = "bold", size = 11),
        legend.position = "right", plot.margin = margin(5, 25, 50, 15))

figure <- (p_depth | p_gc | p_size) / p_names +
  plot_layout(heights = c(1.0, 0.42)) +
  plot_annotation(
    title = "Source-depth distribution and genome properties of 17 Clade 4 / DIRM-like carriers",
    subtitle = "Most recovered carriers originated from 8–12 cm, whereas GC content and genome size did not differ significantly among source-depth groups",
    caption = paste0(
      "Points represent individual MAGs; diamonds indicate group medians. The 4–8 cm group contains only two MAGs, so its violin density is descriptive and should not be over-interpreted. ",
      "The depth count reflects MAG recovery/source names, not environmental abundance. GC and genome-size comparisons used Kruskal–Wallis tests without claiming depth adaptation."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10, color = "#404040"),
      plot.caption = element_text(size = 8.1, color = "#444444", hjust = 0)
    )
  )

save_plot <- function(extension, device, dpi = 450) {
  args <- list(
    filename = file.path(out_dir, paste0("Clade4_source_depth_GC_genome_violin.", extension)),
    plot = figure, width = 16.5, height = 10.8, units = "in", bg = "white", device = device
  )
  if (extension == "png") args$dpi <- dpi
  do.call(ggsave, args)
}
save_plot("pdf", cairo_pdf)
save_plot("svg", grDevices::svg)
save_plot("png", "png")

qc <- c(
  "Clade 4 source-depth / GC / genome-size figure QC",
  "Status: PASSED",
  paste0("Clade 4 MAGs matched to quality table: ", nrow(quality), " (required 17)"),
  paste0("Depth counts 0-4 / 4-8 / 8-12 cm: ", paste(depth_counts$MAG_count, collapse = " / ")),
  paste0("GC Kruskal-Wallis P: ", signif(kw_gc$p.value, 6)),
  paste0("Genome-size Kruskal-Wallis P: ", signif(kw_size$p.value, 6)),
  paste0("Completeness Kruskal-Wallis P: ", signif(kw_completeness$p.value, 6)),
  paste0("Any pairwise GC BH-FDR < 0.05: ", any(c(pairwise_gc$p.value) < 0.05, na.rm = TRUE)),
  paste0("Any pairwise genome-size BH-FDR < 0.05: ", any(c(pairwise_size$p.value) < 0.05, na.rm = TRUE)),
  "The 4-8 cm violin is descriptive because n=2.",
  "Source-depth recovery counts are not interpreted as environmental abundance."
)
writeLines(qc, file.path(out_dir, "Clade4_source_depth_GC_genome_violin_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

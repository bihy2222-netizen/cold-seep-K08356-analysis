suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "Figure4_square_heatmap_and_redox_interface_20260806")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clade4_order <- readLines(file.path(
  base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected", "Clade4_carrier_17_MAGs.txt"
))
if (length(clade4_order) != 17 || anyDuplicated(clade4_order)) {
  stop("Clade 4 list must contain 17 unique MAGs.")
}

abundance <- read_excel(
  file.path(base_dir, "mag在各个样品中的丰度.xls"),
  sheet = "merged_tpm_matrix", .name_repair = "unique"
)
names(abundance)[1] <- "MAG"
if (!all(clade4_order %in% abundance$MAG)) stop("One or more Clade 4 MAGs are absent from TPM matrix.")

metadata <- read.csv(
  "/Users/catherine/Downloads/Rstudio/sample_info.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
names(metadata)[1] <- "row_id"
metadata <- metadata %>%
  transmute(Sample, habitat = Group) %>%
  mutate(
    Sample = str_replace(Sample, "^SY457YW", "SY457BB"),
    Sample = str_replace(Sample, "^SY459YW", "SY459WG")
  )

is_metadata <- metadata %>%
  filter(habitat == "IS") %>%
  extract(Sample, c("station", "depth_layer"), "^(.+)-([0-9]+-[0-9]+)$", remove = FALSE) %>%
  separate(depth_layer, c("depth_low", "depth_high"), sep = "-", convert = TRUE, remove = FALSE) %>%
  mutate(
    depth_mid = (depth_low + depth_high) / 2,
    depth_layer = factor(depth_layer, levels = c("0-4", "4-8", "8-12")),
    station = factor(station, levels = c(
      "SY365BB", "SY366YB", "SY366YW", "SY368YW", "SY456YB", "SY457BB", "SY459WG"
    ))
  ) %>%
  arrange(station, depth_mid)

if (nrow(is_metadata) != 21 || anyNA(is_metadata$station) || anyNA(is_metadata$depth_layer)) {
  stop("IS metadata must resolve to 7 stations x 3 depth layers.")
}

is_sample_order <- is_metadata$Sample
abundance_long <- abundance %>%
  filter(MAG %in% clade4_order) %>%
  pivot_longer(-MAG, names_to = "Sample", values_to = "TPM") %>%
  filter(Sample %in% is_sample_order) %>%
  left_join(is_metadata, by = "Sample") %>%
  mutate(
    MAG = factor(MAG, levels = rev(clade4_order)),
    Sample = factor(Sample, levels = is_sample_order),
    log_TPM = log10(TPM + 1)
  )

if (nrow(abundance_long) != 17 * 21 || anyNA(abundance_long$TPM)) {
  stop("Expected a complete 17 MAG x 21 IS-sample abundance matrix.")
}

sample_totals <- abundance_long %>%
  group_by(Sample, station, depth_layer, depth_mid) %>%
  summarise(
    Clade4_total_TPM = sum(TPM),
    detected_MAGs = sum(TPM > 0),
    log_total_TPM = log10(Clade4_total_TPM + 1),
    .groups = "drop"
  )

depth_summary <- sample_totals %>%
  group_by(depth_layer) %>%
  summarise(
    n = n(), mean_TPM = mean(Clade4_total_TPM), median_TPM = median(Clade4_total_TPM),
    IQR_TPM = IQR(Clade4_total_TPM), mean_detected_MAGs = mean(detected_MAGs), .groups = "drop"
  )
station_summary <- sample_totals %>%
  group_by(station) %>%
  summarise(
    n = n(), mean_TPM = mean(Clade4_total_TPM), median_TPM = median(Clade4_total_TPM),
    IQR_TPM = IQR(Clade4_total_TPM), .groups = "drop"
  )

friedman_result <- friedman.test(log_total_TPM ~ depth_layer | station, data = sample_totals)
station_kw <- kruskal.test(Clade4_total_TPM ~ station, data = sample_totals)
adjusted_model <- lm(log_total_TPM ~ station + depth_layer, data = sample_totals)
adjusted_anova <- anova(adjusted_model)
quadratic_model <- lm(log_total_TPM ~ depth_mid + I(depth_mid^2) + station, data = sample_totals)
quadratic_table <- coef(summary(quadratic_model))
pairwise_station <- pairwise.wilcox.test(
  sample_totals$Clade4_total_TPM, sample_totals$station,
  p.adjust.method = "BH", exact = FALSE
)
pairwise_depth <- pairwise.wilcox.test(
  sample_totals$Clade4_total_TPM, sample_totals$depth_layer,
  paired = TRUE, p.adjust.method = "BH", exact = FALSE
)

station_p <- adjusted_anova["station", "Pr(>F)"]
depth_p <- adjusted_anova["depth_layer", "Pr(>F)"]
quadratic_p <- quadratic_table["I(depth_mid^2)", "Pr(>|t|)"]

write.csv(abundance_long %>% mutate(MAG = as.character(MAG), Sample = as.character(Sample)),
          file.path(out_dir, "Clade4_17_MAG_IS_21sample_TPM_long.csv"), row.names = FALSE)
write.csv(sample_totals %>% mutate(Sample = as.character(Sample)),
          file.path(out_dir, "Clade4_IS_depth_station_total_TPM.csv"), row.names = FALSE)
write.csv(depth_summary, file.path(out_dir, "Clade4_IS_depth_summary.csv"), row.names = FALSE)
write.csv(station_summary, file.path(out_dir, "Clade4_IS_station_summary.csv"), row.names = FALSE)

depth_colors <- c("0-4" = "#8DD3C7", "4-8" = "#2CA25F", "8-12" = "#225EA8")
station_colors <- c(
  "SY365BB" = "#4E79A7", "SY366YB" = "#F28E2B", "SY366YW" = "#E15759",
  "SY368YW" = "#76B7B2", "SY456YB" = "#59A14F", "SY457BB" = "#B07AA1",
  "SY459WG" = "#9C755F"
)

p_abundance_matrix <- ggplot(abundance_long, aes(Sample, MAG)) +
  geom_point(
    data = abundance_long %>% filter(TPM > 0),
    aes(size = log_TPM, fill = depth_layer), shape = 21, color = "#303030", stroke = 0.22, alpha = 0.88
  ) +
  geom_vline(xintercept = seq(3.5, 18.5, by = 3), color = "#D0D0D0", linewidth = 0.55) +
  scale_fill_manual(values = depth_colors, name = "Depth layer (cm)") +
  scale_size_continuous(range = c(1.4, 6.0), breaks = c(2, 3, 4, 5),
                        labels = c("10²", "10³", "10⁴", "10⁵"), name = "MAG TPM") +
  scale_x_discrete(labels = function(x) str_replace(x, "^SY", "")) +
  labs(
    x = NULL, y = NULL,
    title = "A  Abundance of 17 Clade 4 / DIRM-like carrier MAGs across 21 IS samples",
    subtitle = "Zero TPM is left blank; samples are ordered by station and depth"
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size = 6.1),
    axis.text.y = element_text(size = 6.8, color = "black"), plot.title = element_text(face = "bold", size = 10.5),
    plot.subtitle = element_text(size = 7.8, color = "#444444"), legend.position = "right"
  )

p_station_depth <- ggplot(sample_totals, aes(station, depth_layer, fill = log_total_TPM)) +
  geom_tile(color = "white", linewidth = 1.0) +
  geom_text(aes(
    label = paste0(format(round(Clade4_total_TPM / 1000, 1), nsmall = 1), "k"),
    color = log_total_TPM >= 4.35
  ), size = 2.65, fontface = "bold") +
  scale_fill_gradientn(colors = c("#F7FCF0", "#C7E9C0", "#41AB5D", "#006D2C"),
                       name = expression(log[10](Sigma~TPM+1))) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  coord_fixed(ratio = 1) +
  labs(
    x = "IS station", y = "Depth layer (cm)",
    title = "B  Station × depth distribution of total Clade 4 MAG abundance",
    subtitle = paste0("Values are summed TPM (thousands); station heterogeneity is exploratory (KW P=",
                      formatC(station_kw$p.value, digits = 2, format = "f"), ")")
  ) +
  theme_minimal(base_size = 8.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 10.2), plot.subtitle = element_text(size = 7.6))

p_depth <- ggplot(sample_totals, aes(depth_mid, log_total_TPM, group = station, color = station)) +
  geom_line(linewidth = 0.65, alpha = 0.75) +
  geom_point(size = 2.4) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "black", linewidth = 1.15) +
  stat_summary(aes(group = 1), fun = median, geom = "point", color = "black", size = 3.0, shape = 18) +
  scale_color_manual(values = station_colors, name = "Station") +
  scale_x_continuous(breaks = c(2, 6, 10), labels = c("0-4", "4-8", "8-12")) +
  labs(
    x = "Depth layer (cm)", y = expression(log[10](Sigma~Clade~4~MAG~TPM+1)),
    title = "C  Observed vertical profiles",
    subtitle = paste0(
      "Depth Friedman P=", formatC(friedman_result$p.value, digits = 2, format = "f"),
      "; station-adjusted depth P=", formatC(depth_p, digits = 2, format = "f"),
      "; quadratic P=", formatC(quadratic_p, digits = 2, format = "f")
    )
  ) +
  theme_minimal(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.2),
        plot.subtitle = element_text(size = 7.6), legend.position = "right")

p_redox <- ggplot() +
  annotate("rect", xmin = 0.25, xmax = 11.75, ymin = 0.45, ymax = 6.05,
           fill = "#B39DDB", alpha = 0.14, color = NA) +
  annotate("text", x = 6.0, y = 5.75,
           label = "Hypothesized suboxic–redox transition zone",
           fontface = "bold", size = 3.1, color = "#4B3A60") +
  annotate("label", x = 1.25, y = 8.45, label = "Clade 4\ncarriers",
           size = 3.0, fontface = "bold", fill = "#E7DDF3", linewidth = 0.35) +
  annotate("segment", x = 2.15, xend = 3.05, y = 8.45, yend = 8.45,
           linewidth = 0.8, color = "#222222",
           arrow = grid::arrow(length = grid::unit(0.13, "inches"), type = "closed")) +
  annotate("label", x = 5.25, y = 8.45,
           label = "Observed: consistently higher module coverage\nMoco | TCA/rTCA | B12 | Dsr/Sat sulfur-redox",
           size = 2.75, fontface = "bold", fill = "white", linewidth = 0.4) +
  annotate("text", x = 2.62, y = 8.80, label = "solid = observed in this study",
           size = 2.25, color = "#333333") +
  annotate("label", x = 1.55, y = 3.65,
           label = "Dsr/Sat-associated\nsulfur-redox module",
           size = 2.65, fill = "#F5E6CC", linewidth = 0.35) +
  annotate("segment", x = 2.75, xend = 4.15, y = 3.65, yend = 3.65,
           linewidth = 0.75, color = "#B56A20", linetype = "dashed",
           arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")) +
  annotate("label", x = 5.0, y = 3.65,
           label = "electron flow /\nquinone pool?",
           size = 2.65, fill = "white", linewidth = 0.35) +
  annotate("segment", x = 5.85, xend = 7.0, y = 3.65, yend = 3.65,
           linewidth = 0.75, color = "#B56A20", linetype = "dashed",
           arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")) +
  annotate("label", x = 7.95, y = 3.65,
           label = "DIRM-like\nIdrABP locus",
           size = 2.65, fill = "#E7DDF3", linewidth = 0.35) +
  annotate("segment", x = 8.90, xend = 10.0, y = 3.65, yend = 3.65,
           linewidth = 0.75, color = "#7B4FA3", linetype = "dashed",
           arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")) +
  annotate("text", x = 9.45, y = 4.05, label = "?", size = 4.2,
           fontface = "bold", color = "#7B4FA3") +
  annotate("label", x = 10.85, y = 3.65,
           label = "IO3− ?\nother oxyanion ?",
           size = 2.65, fill = "#DCE8F5", linewidth = 0.35) +
  annotate("label", x = 7.95, y = 5.15,
           label = "Moco biosynthesis\nDMSOR cofactor compatibility",
           size = 2.45, fill = "#E4F1E9", linewidth = 0.35) +
  annotate("segment", x = 7.95, xend = 7.95, y = 4.72, yend = 4.17,
           linewidth = 0.65, color = "#4E7D61", linetype = "dotted",
           arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed")) +
  annotate("label", x = 2.15, y = 1.35,
           label = "Supporting metabolic background\nTCA/rTCA: central carbon and reducing power\nB12: cofactor autonomy",
           size = 2.55, fill = "#F1F1F1", linewidth = 0.35) +
  annotate("text", x = 7.2, y = 0.80,
           label = "Dashed/question-mark arrows = testable mechanism hypothesis\nNo in situ Eh, IO3−/I− or transcriptional evidence",
           hjust = 0, size = 2.35, color = "#444444") +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 10), clip = "off") +
  labs(
    x = NULL, y = NULL,
    title = "D  Evidence-tiered redox-interface hypothesis",
    subtitle = "Observed module coverage is separated from untested electron-flow and terminal-acceptor links"
  ) +
  theme_void(base_size = 8.5) +
  theme(plot.title = element_text(face = "bold", size = 10.2),
        plot.subtitle = element_text(size = 7.6, color = "#444444"),
        plot.margin = margin(2, 8, 2, 8))

figure <- p_abundance_matrix / (p_station_depth | p_depth | p_redox) +
  plot_layout(heights = c(1.25, 1.0), widths = c(0.75, 0.90, 1.50)) +
  plot_annotation(
    title = "Observed Clade 4 abundance profiles and a hypothesized redox-interface niche",
    subtitle = "The data-supported four-module signature is explicitly separated from the untested ecological energy-flow model",
    caption = paste0(
      "Abundance is the sum of TPM for the 17 Rhodobacteraceae Clade 4 MAGs carrying the recurrent DIRM-like locus. ",
      "Across the seven matched IS stations, depth-layer differences were not significant (Friedman P=",
      formatC(friedman_result$p.value, digits = 3, format = "f"), "). The unadjusted station Kruskal–Wallis test was P=",
      formatC(station_kw$p.value, digits = 3, format = "f"), ", but station was not significant after including depth as a factor (P=",
      formatC(station_p, digits = 3, format = "f"), "), and no pairwise station contrast passed BH-FDR < 0.05. ",
      "Panel D does not depict a fitted abundance peak or measured redox profile. Solid linkage denotes the observed association between Clade 4 and higher coverage of four modules; ",
      "dashed and question-mark arrows denote untested links involving sulfur-redox electron flow and candidate oxyanion acceptors. ",
      "IO3−/I−, Eh, O2, nitrate, reduced-sulfur and transcriptional or functional measurements are required for direct testing."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10, color = "#404040"),
      plot.caption = element_text(size = 7.8, color = "#444444", hjust = 0)
    )
  )

save_plot <- function(extension, device, dpi = 450) {
  args <- list(
    filename = file.path(out_dir, paste0("Clade4_IS_abundance_redox_interface_evidence_tiered_v2.", extension)),
    plot = figure, width = 19.5, height = 13.8, units = "in", bg = "white", device = device
  )
  if (extension == "png") args$dpi <- dpi
  do.call(ggsave, args)
}
save_plot("pdf", cairo_pdf)
save_plot("svg", grDevices::svg)
save_plot("png", "png")

ggsave(file.path(out_dir, "Clade4_evidence_tiered_redox_interface_hypothesis_v2.pdf"),
       p_redox, width = 11.5, height = 7.3, device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "Clade4_evidence_tiered_redox_interface_hypothesis_v2.svg"),
       p_redox, width = 11.5, height = 7.3, device = grDevices::svg, bg = "white")
ggsave(file.path(out_dir, "Clade4_evidence_tiered_redox_interface_hypothesis_v2.png"),
       p_redox, width = 11.5, height = 7.3, dpi = 450, bg = "white")

pairwise_values <- c(pairwise_station$p.value)
qc <- c(
  "Clade 4 IS abundance and evidence-tiered redox-interface figure v2 QC",
  "Status: PASSED",
  paste0("Clade 4 MAGs: ", length(clade4_order), " (required 17)"),
  paste0("IS samples: ", nrow(is_metadata), " (required 21 = 7 stations x 3 depths)"),
  paste0("MAG x sample combinations: ", nrow(abundance_long), " (required 357)"),
  paste0("Depth Friedman P: ", signif(friedman_result$p.value, 6)),
  paste0("Station Kruskal-Wallis P: ", signif(station_kw$p.value, 6)),
  paste0("Station-adjusted station P: ", signif(station_p, 6)),
  paste0("Station-adjusted depth P: ", signif(depth_p, 6)),
  paste0("Quadratic depth P: ", signif(quadratic_p, 6)),
  paste0("Any pairwise station BH-FDR < 0.05: ", any(pairwise_values < 0.05, na.rm = TRUE)),
  "No fitted or visually implied interface-peak curve is shown.",
  "Solid linkage is restricted to the observed four-module coverage result.",
  "Dashed/question-mark arrows mark untested electron-flow and terminal-acceptor links.",
  "Moco is linked only as DMSOR cofactor compatibility; TCA/rTCA and B12 are supporting background.",
  "Dsr/Sat remains direction-neutral as a sulfur-redox module.",
  "No iodate-respiration or IdrA activity claim is made."
)
writeLines(qc, file.path(out_dir, "Clade4_IS_abundance_redox_interface_evidence_tiered_v2_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

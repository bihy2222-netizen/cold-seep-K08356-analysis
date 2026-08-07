suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "Clade4_redox_potential_schematic_20260806")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

potential_data <- data.frame(
  couple = c(
    "IO3− / I−", "NO3− / NO2−", "As(V) / As(III)",
    "menaquinone / menaquinol", "NAD+ / NADH"
  ),
  midpoint_mV = c(700, 420, 60, -80, -320),
  class = c("candidate", "context", "context", "carrier", "donor"),
  x_start = c(1.25, 1.25, 1.25, 1.25, 1.25),
  x_end = c(5.30, 4.80, 4.55, 4.55, 4.55),
  stringsAsFactors = FALSE
)

class_colors <- c(
  candidate = "#7B4FA3", context = "#547AA5", carrier = "#5D6D7E", donor = "#B56A20"
)

p_ladder <- ggplot() +
  annotate("rect", xmin = 0.65, xmax = 5.65, ymin = -500, ymax = -100,
           fill = "#E6A85C", alpha = 0.13, color = "#C47A20", linewidth = 0.45,
           linetype = "dashed") +
  annotate("text", x = 3.15, y = -130,
           label = "Dsr/Sat-associated sulfur-redox window",
           fontface = "bold", size = 3.5, color = "#8A531B") +
  annotate("text", x = 3.15, y = -180,
           label = "Direction unresolved; potential range is substrate- and condition-dependent",
           size = 2.8, color = "#6A4A2A") +
  annotate("segment", x = 1.35, xend = 4.85, y = -270, yend = -270,
           linewidth = 1.0, color = "#C47A20", linetype = "dashed") +
  annotate("point", x = 1.35, y = -270, shape = 21, size = 4.5,
           fill = "#F0B35F", color = "#8A531B") +
  annotate("point", x = 4.85, y = -270, shape = 21, size = 4.5,
           fill = "#F0B35F", color = "#8A531B") +
  annotate("text", x = 1.35, y = -315, label = "reduced sulfur\n(e.g., HS−)", size = 2.9) +
  annotate("text", x = 4.85, y = -315, label = "oxidized sulfur\nintermediates", size = 2.9) +
  geom_segment(
    data = potential_data,
    aes(x = x_start, xend = x_end, y = midpoint_mV, yend = midpoint_mV, color = class),
    linewidth = 1.15,
    arrow = grid::arrow(length = grid::unit(0.11, "inches"), type = "closed")
  ) +
  geom_point(
    data = potential_data,
    aes(x = x_start, y = midpoint_mV, fill = class),
    shape = 21, size = 4.2, color = "#303030", stroke = 0.55
  ) +
  geom_point(
    data = potential_data,
    aes(x = x_end, y = midpoint_mV, fill = class),
    shape = 21, size = 4.2, color = "#303030", stroke = 0.55
  ) +
  geom_text(
    data = potential_data,
    aes(x = (x_start + x_end) / 2, y = midpoint_mV + 38, label = couple),
    fontface = "bold", size = 3.6, color = "#202020"
  ) +
  annotate("text", x = 5.45, y = 735, label = "candidate acceptor ?",
           hjust = 1, size = 3.1, color = "#7B4FA3", fontface = "bold") +
  annotate("segment", x = 6.25, xend = 6.25, y = -300, yend = 675,
           linewidth = 1.25, color = "#2B2B2B", linetype = "dashed",
           arrow = grid::arrow(length = grid::unit(0.16, "inches"), type = "closed")) +
  annotate("text", x = 6.43, y = 190,
           label = "positive ΔE\n(hypothesized)", angle = 90,
           size = 3.2, fontface = "bold", color = "#2B2B2B") +
  annotate("text", x = 6.43, y = -350,
           label = "ΔG = −nFΔE", angle = 90, size = 3.0, color = "#444444") +
  scale_color_manual(values = class_colors, guide = "none") +
  scale_fill_manual(values = class_colors, guide = "none") +
  scale_y_continuous(
    limits = c(-550, 820), breaks = c(-500, -250, 0, 250, 500, 750),
    labels = function(x) paste0(x, " mV"), expand = c(0, 0)
  ) +
  scale_x_continuous(limits = c(0, 6.9), expand = c(0, 0)) +
  labs(
    x = NULL, y = "Illustrative standard midpoint potential, E°′ (mV)",
    title = "A  Conceptual redox-potential ladder",
    subtitle = "Vertical positions are approximate literature values/ranges—not measured Eh in the South China Sea samples"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#D8D8D8", linewidth = 0.45),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12.5),
    plot.subtitle = element_text(size = 8.7, color = "#444444"),
    plot.margin = margin(6, 12, 6, 6)
  )

p_mechanism <- ggplot() +
  annotate("rect", xmin = 0.2, xmax = 9.8, ymin = 0.25, ymax = 5.85,
           fill = "#B39DDB", alpha = 0.12, color = NA) +
  annotate("text", x = 5.0, y = 5.52,
           label = "Hypothesized suboxic–redox transition zone",
           fontface = "bold", size = 3.5, color = "#4B3A60") +
  annotate("label", x = 1.55, y = 8.45, label = "Clade 4 carriers",
           size = 3.2, fontface = "bold", fill = "#E7DDF3", linewidth = 0.4) +
  annotate("segment", x = 2.65, xend = 3.75, y = 8.45, yend = 8.45,
           linewidth = 0.95, color = "#222222",
           arrow = grid::arrow(length = grid::unit(0.14, "inches"), type = "closed")) +
  annotate("label", x = 6.35, y = 8.45,
           label = "Observed: higher genomic coverage\nMoco | TCA/rTCA | B12 | Dsr/Sat sulfur-redox",
           size = 3.0, fontface = "bold", fill = "white", linewidth = 0.4) +
  annotate("text", x = 3.20, y = 8.82, label = "solid = observed association",
           size = 2.45, color = "#333333") +
  annotate("label", x = 1.35, y = 3.70,
           label = "low-potential\nsulfur-redox pool",
           size = 2.9, fill = "#F5E6CC", linewidth = 0.35) +
  annotate("segment", x = 2.35, xend = 3.40, y = 3.70, yend = 3.70,
           linewidth = 0.85, color = "#B56A20", linetype = "dashed",
           arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")) +
  annotate("label", x = 4.35, y = 3.70, label = "electron flow /\nquinone pool ?",
           size = 2.9, fill = "white", linewidth = 0.35) +
  annotate("segment", x = 5.25, xend = 6.25, y = 3.70, yend = 3.70,
           linewidth = 0.85, color = "#B56A20", linetype = "dashed",
           arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")) +
  annotate("label", x = 7.10, y = 3.70, label = "DIRM-like\nIdrABP locus",
           size = 2.9, fill = "#E7DDF3", linewidth = 0.35) +
  annotate("segment", x = 7.95, xend = 8.65, y = 3.70, yend = 3.70,
           linewidth = 0.85, color = "#7B4FA3", linetype = "dashed",
           arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")) +
  annotate("text", x = 8.30, y = 4.05, label = "?", size = 4.4,
           fontface = "bold", color = "#7B4FA3") +
  annotate("label", x = 9.25, y = 3.70, label = "IO3− ?\nother oxyanion ?",
           size = 2.8, fill = "#DCE8F5", linewidth = 0.35) +
  annotate("label", x = 7.10, y = 5.00,
           label = "Moco biosynthesis\nDMSOR cofactor compatibility",
           size = 2.55, fill = "#E4F1E9", linewidth = 0.35) +
  annotate("segment", x = 7.10, xend = 7.10, y = 4.62, yend = 4.14,
           linewidth = 0.65, color = "#4E7D61", linetype = "dotted",
           arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed")) +
  annotate("label", x = 2.55, y = 1.25,
           label = "Supporting metabolic background\nTCA/rTCA: central carbon and reducing power\nB12: cofactor autonomy",
           size = 2.75, fill = "#F1F1F1", linewidth = 0.35) +
  annotate("text", x = 6.0, y = 0.78,
           label = "Dashed/question-mark arrows = testable hypothesis\nNo measured Eh, IO3−/I−, transcription or activity evidence",
           hjust = 0, size = 2.5, color = "#444444") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10), clip = "off") +
  labs(
    x = NULL, y = NULL,
    title = "B  Evidence-tiered ecological energy-flow hypothesis",
    subtitle = "The potential ladder motivates testable links but does not identify the DIRM-like substrate"
  ) +
  theme_void(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5),
    plot.subtitle = element_text(size = 8.7, color = "#444444"),
    plot.margin = margin(6, 6, 6, 12)
  )

figure <- (p_ladder | p_mechanism) +
  plot_layout(widths = c(1.0, 1.15)) +
  plot_annotation(
    title = "Conceptual redox-potential framework for the Clade 4 / DIRM-like lineage",
    subtitle = "Approximate redox couples provide energetic context for a hypothesis linking sulfur-redox metabolism to candidate high-potential oxyanion acceptors",
    caption = paste0(
      "This is an original conceptual schematic, not a reconstruction of an evolutionary timeline. Midpoint potentials are approximate standard/literature values and are sensitive to pH, concentrations and reaction definitions; ",
      "they are not in situ Eh measurements. The Dsr/Sat-associated module is direction-neutral. The IO3−/I− link is explicitly hypothetical, and the DIRM-like locus is not identified as an iodate reductase without substrate, product, expression or activity evidence. ",
      "Solid linkage denotes the observed association between Clade 4 and higher coverage of four modules; dashed and question-mark links denote untested mechanism steps."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 10, color = "#404040"),
      plot.caption = element_text(size = 8.0, color = "#444444", hjust = 0)
    )
  )

save_plot <- function(extension, device, dpi = 450) {
  args <- list(
    filename = file.path(out_dir, paste0("Clade4_conceptual_redox_potential_framework.", extension)),
    plot = figure, width = 18.5, height = 10.8, units = "in", bg = "white", device = device
  )
  if (extension == "png") args$dpi <- dpi
  do.call(ggsave, args)
}
save_plot("pdf", cairo_pdf)
save_plot("svg", grDevices::svg)
save_plot("png", "png")

write.csv(potential_data, file.path(out_dir, "redox_potential_values_used.csv"), row.names = FALSE)

qc <- c(
  "Clade 4 conceptual redox-potential framework QC",
  "Status: PASSED",
  "The figure is original and does not reproduce the reference geological timeline.",
  "Potential values are explicitly labeled illustrative/approximate and not measured sample Eh.",
  "Dsr/Sat sulfur-redox direction remains unresolved.",
  "IO3-/I- and other oxyanion links use question marks and dashed arrows.",
  "Solid linkage is restricted to observed higher module coverage.",
  "TCA/rTCA and B12 are supporting background, not direct electron-transfer steps.",
  "No iodate-respiration or IdrA activity claim is made."
)
writeLines(qc, file.path(out_dir, "Clade4_conceptual_redox_potential_framework_QC.txt"))
cat(paste(qc, collapse = "\n"), "\n")

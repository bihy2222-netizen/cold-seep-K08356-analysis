suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

out_dir <- "/Users/catherine/Downloads/MAG-bin-tpm/K08356_MAG_gene_island_arrows_20260708"
coords_file <- file.path(out_dir, "K08356_four_branch_gene_neighborhood_coordinates.csv")
fig_subdir <- Sys.getenv("K08356_FIG_SUBDIR", "four_branch_gene_island_figures")
show_domain_markers <- tolower(Sys.getenv("K08356_SHOW_DOMAIN_MARKERS", "true")) %in% c("1", "true", "yes", "y")
show_gene_text_labels <- tolower(Sys.getenv("K08356_SHOW_GENE_TEXT_LABELS", "true")) %in% c("1", "true", "yes", "y")
show_facet_y_labels <- tolower(Sys.getenv("K08356_SHOW_FACET_Y_LABELS", "false")) %in% c("1", "true", "yes", "y")
legend_position <- Sys.getenv("K08356_LEGEND_POSITION", "bottom")
selected_ids_file <- Sys.getenv("K08356_SELECTED_IDS_FILE", "")
fig_dir <- file.path(out_dir, fig_subdir)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

branch_order <- c(
  "IdrA-associated",
  "canonical aioA-associated",
  "aioA-like-associated",
  "unknown AioA-like / uncertain DMSOR"
)

class_colors <- c(
  "Target IdrA-like K08356" = "#5B8DB8",
  "Target canonical aioA" = "#2CA25F",
  "Target aioA-like K08356" = "#009E73",
  "Target uncertain AioA-like" = "#CC79A7",
  "Rieske/AioB" = "#B9A7D8",
  "Cytochrome/peroxidase" = "#9E9E9E",
  "Arsenic resistance" = "#A8A8A8",
  "Transporter" = "#BDBDBD",
  "Nitrogen metabolism" = "#969696",
  "Sox/sulfur oxidation" = "#8A8A8A",
  "Redox/oxidoreductase" = "#C7C7C7",
  "Annotated other" = "#D6D6D6",
  "Hypothetical/unknown" = "#EFEFEF"
)

domain_marker_colors <- c(
  "Molybdopterin oxidoreductase / K08356" = "#E41A1C",
  "Rieske [2Fe-2S] / AioB" = "#FFB000",
  "Cytochrome c peroxidase" = "#009E3D"
)

coords <- read_csv(coords_file, show_col_types = FALSE) %>%
  mutate(
    Feature = factor(Feature, levels = branch_order),
    Molecule_short = paste0(Original_MAG_ID, " | ", Corrected_habitat),
    Gene_class = factor(Gene_class, levels = names(class_colors)),
    Label_to_plot = ifelse(Gene_label == "" | is.na(Gene_label), NA_character_, Gene_label),
    Domain_marker = case_when(
      Is_target == "yes" ~ "Molybdopterin oxidoreductase / K08356",
      Gene_class == "Rieske/AioB" ~ "Rieske [2Fe-2S] / AioB",
      Gene_class == "Cytochrome/peroxidase" ~ "Cytochrome c peroxidase",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Molecule_short) %>%
  mutate(
    Region_start = min(Start),
    Start_plot = Start - Region_start + 1,
    End_plot = End - Region_start + 1,
    Mid_plot = (Start_plot + End_plot) / 2
  ) %>%
  ungroup()

if (nzchar(selected_ids_file)) {
  selected_ids <- read_lines(selected_ids_file)
  selected_ids <- selected_ids[nzchar(selected_ids)]
  coords <- coords %>% filter(Original_MAG_ID %in% selected_ids)
}

make_arrow_poly <- function(df, y_map, head_frac = 0.18) {
  pieces <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    xmin <- r$Start_plot
    xmax <- r$End_plot
    y <- y_map[[r$Molecule_key]]
    h <- 0.16
    len <- max(xmax - xmin, 1)
    head <- min(len * head_frac, len * 0.45, 450)
    if (r$Direction >= 0) {
      xs <- c(xmin, xmax - head, xmax - head, xmax, xmax - head, xmax - head, xmin)
    } else {
      xs <- c(xmax, xmin + head, xmin + head, xmin, xmin + head, xmin + head, xmax)
    }
    ys <- c(y - h, y - h, y - h * 1.55, y, y + h * 1.55, y + h, y + h)
    pieces[[i]] <- data.frame(
      row_id = r$row_id,
      Feature = r$Feature,
      Molecule_key = r$Molecule_key,
      Gene_class = r$Gene_class,
      x = xs,
      y = ys
    )
  }
  bind_rows(pieces)
}

make_plot <- function(dat, title, subtitle, height, filename_base, facet = TRUE) {
  mol_levels <- dat %>%
    distinct(Feature, Original_MAG_ID, Corrected_habitat, Molecule_short) %>%
    arrange(Feature, Corrected_habitat, Original_MAG_ID) %>%
    pull(Molecule_short)
  dat <- dat %>%
    mutate(
      Molecule_short = factor(Molecule_short, levels = rev(unique(mol_levels))),
      Molecule_key = as.character(Molecule_short),
      row_id = row_number()
    )
  y_map <- setNames(seq_along(levels(dat$Molecule_short)), levels(dat$Molecule_short))
  poly_df <- make_arrow_poly(dat, y_map)
  label_df <- dat %>%
    filter(!is.na(Label_to_plot)) %>%
    mutate(y = y_map[as.character(Molecule_short)] + 0.62)
  marker_df <- dat %>%
    filter(!is.na(Domain_marker)) %>%
    mutate(y = y_map[as.character(Molecule_short)] + 0.36)
  baseline_df <- dat %>%
    group_by(Feature, Molecule_short, Molecule_key) %>%
    summarise(xmin = min(Start_plot), xmax = max(End_plot), y = y_map[as.character(first(Molecule_short))], .groups = "drop")

  p <- ggplot() +
    geom_segment(
      data = baseline_df,
      aes(x = xmin, xend = xmax, y = y, yend = y),
      color = "grey72",
      linewidth = 0.28
    ) +
    geom_polygon(
      data = poly_df,
      aes(x = x, y = y, group = row_id, fill = Gene_class),
      color = "grey25",
      linewidth = 0.16
    )
  if (show_domain_markers) {
    p <- p +
      geom_point(
        data = marker_df,
        aes(x = Mid_plot, y = y, color = Domain_marker),
        size = 1.45,
        stroke = 0.35
      )
  }
  if (show_gene_text_labels) {
    p <- p +
      geom_text(
        data = label_df,
        aes(x = Mid_plot, y = y, label = Label_to_plot),
        family = "Times",
        size = 1.45
      )
  }
  p <- p +
    scale_fill_manual(values = class_colors, drop = TRUE) +
    scale_y_continuous(
      breaks = seq_along(levels(dat$Molecule_short)),
      labels = levels(dat$Molecule_short),
      expand = expansion(add = 1.15)
    ) +
    scale_x_continuous(labels = function(x) paste0(x / 1000, "k")) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Relative genomic position in target-centered neighborhood (bp)",
      y = NULL,
      fill = "Gene/function class",
      caption = "Colored arrows indicate K08356 branch targets; grey arrows indicate neighboring ORFs. Circle markers denote key functional domains."
    ) +
    theme_classic(base_family = "Times", base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 5.2),
      axis.text.x = element_text(size = 8),
      axis.title.x = element_text(face = "bold"),
      legend.position = legend_position,
      legend.title = element_text(face = "bold"),
      plot.caption = element_text(size = 8, hjust = 0.5),
      panel.grid.major.x = element_line(color = "grey92", linetype = "dashed")
    )
  if (show_domain_markers) {
    p <- p +
      scale_color_manual(values = domain_marker_colors, na.translate = FALSE, name = "Functional domain marker")
  }
  if (facet) {
    p <- p +
      facet_grid(Feature ~ ., scales = "free_y", space = "free_y") +
      theme(
        strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.3),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.y = if (show_facet_y_labels) element_text(size = 4.8) else element_blank(),
        axis.ticks.y = if (show_facet_y_labels) element_line(linewidth = 0.2) else element_blank()
      )
  }
  ggsave(file.path(fig_dir, paste0(filename_base, ".png")), p, width = 14, height = height, dpi = 450, bg = "white")
  ggsave(file.path(fig_dir, paste0(filename_base, ".pdf")), p, width = 14, height = height, device = cairo_pdf, bg = "white", family = "Times")
  if (filename_base == "K08356_four_branch_MAG_gene_island_arrows_all") {
    ggsave(file.path(fig_dir, paste0(filename_base, ".svg")), p, width = 14, height = height, device = svg, bg = "white", family = "Times")
    ggsave(file.path(fig_dir, paste0(filename_base, ".eps")), p, width = 14, height = height, device = cairo_ps, bg = "white", family = "Times")
    ggsave(file.path(fig_dir, paste0(filename_base, ".tif")), p, width = 14, height = height, dpi = 450, compression = "lzw", bg = "white")
  }
}

make_plot(
  coords,
  "K08356 MAG gene neighborhoods across four protein-branch groups",
  "Target-centered neighborhoods from downloaded sample-level CDS files; up to 10 upstream and 10 downstream ORFs are shown",
  height = 26,
  filename_base = "K08356_four_branch_MAG_gene_island_arrows_all",
  facet = TRUE
)

for (br in branch_order) {
  dat <- coords %>% filter(as.character(Feature) == br)
  h <- max(4.8, 1.25 + 0.58 * n_distinct(dat$Molecule_short))
  safe <- str_replace_all(br, "[^A-Za-z0-9]+", "_") %>% str_replace_all("_+$", "")
  make_plot(
    dat,
    paste0("K08356 MAG gene neighborhoods: ", br),
    "Target-centered neighborhoods; up to 10 upstream and 10 downstream ORFs are shown",
    height = h,
    filename_base = paste0("K08356_MAG_gene_island_arrows_", safe),
    facet = FALSE
  )
}

summary <- coords %>%
  group_by(Feature, Original_MAG_ID, Corrected_habitat, Target_contig) %>%
  summarise(
    n_orfs_plotted = n(),
    target_found = any(Is_target == "yes"),
    n_labeled_neighbors = sum(!is.na(Label_to_plot)),
    .groups = "drop"
  )
write_csv(summary, file.path(out_dir, "K08356_four_branch_gene_neighborhood_plot_summary.csv"))

message("Wrote gene-island figures to: ", fig_dir)

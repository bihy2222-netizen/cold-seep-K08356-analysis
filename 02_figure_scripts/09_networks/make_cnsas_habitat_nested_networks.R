suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(writexl)
})

# CNSAs functional-gene co-occurrence networks by habitat.
# Raw abundance source: all_koid_tpm-final.xlsx (contig-level TPM).
# Gene-to-module mapping and KO aggregation are performed by make_cnsas_figures.R.

base_dir <- getwd()
raw_file <- file.path(base_dir, "all_koid_tpm-final.xlsx")
prep_script <- file.path(base_dir, "make_cnsas_figures.R")
long_file <- file.path(base_dir, "outputs", "cnsas_contig_figures_story_final",
                       "CNSAs_contig_TPM_long_table.csv")
out_dir <- file.path(base_dir, "outputs", "cnsas_habitat_nested_networks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_file)) stop("Input file not found: ", raw_file)
if (!file.exists(long_file) || file.mtime(long_file) < file.mtime(raw_file)) {
  if (!file.exists(prep_script)) stop("Preprocessing script not found: ", prep_script)
  message("Preparing the gene-level table from all_koid_tpm-final.xlsx ...")
  source(prep_script, local = new.env(parent = globalenv()))
}

dat <- read.csv(long_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS")))

required <- c("module", "gene", "sample", "habitat", "TPM", "log10_TPM_plus1")
if (!all(required %in% names(dat))) {
  stop("Prepared table lacks columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
}

module_order <- c("Methane", "Nitrogen", "Sulfur", "Arsenic")
module_colors <- c(Methane = "#005B96", Nitrogen = "#007F5F",
                   Sulfur = "#D8A600", Arsenic = "#5A189A")
module_labels <- c(Methane = "Carbon", Nitrogen = "Nitrogen",
                   Sulfur = "Sulfur", Arsenic = "Arsenic")
habitat_labels <- c(IS = "Ice sediment (IS)", AS = "Active sediment (AS)",
                    ES = "Estuarine sediment (ES)", NS = "Nearshore sediment (NS)")

# A fixed nested-circle layout makes node positions directly comparable among panels.
gene_meta <- dat %>%
  group_by(gene, module) %>%
  summarise(overall_mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(module = factor(module, levels = module_order)) %>%
  arrange(module, gene)

centers <- tibble(
  module = factor(module_order, levels = module_order),
  center_x = c(-7.2, 0, 7.2, 0),
  center_y = c(0, 5.2, 0, -5.2)
)

node_pos <- gene_meta %>%
  left_join(centers, by = "module") %>%
  group_by(module) %>%
  arrange(gene, .by_group = TRUE) %>%
  mutate(node_rank = row_number(), nodes_in_module = n(),
         angle = pi / 2 + 2 * pi * (node_rank - 1) / nodes_in_module,
         local_radius = 1.55 + 0.055 * nodes_in_module,
         x = center_x + local_radius * cos(angle),
         y = center_y + local_radius * sin(angle)) %>%
  ungroup()

habitat_node_abundance <- dat %>%
  group_by(habitat, gene) %>%
  summarise(habitat_mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  left_join(node_pos, by = "gene")

# One shared scale allows valid visual comparison while each panel uses its own
# within-habitat mean abundance.
global_size_limits <- c(0, max(habitat_node_abundance$habitat_mean_TPM, na.rm = TRUE))

cor_edges <- function(habitat_code) {
  wide <- dat %>%
    filter(habitat == habitat_code) %>%
    select(sample, gene, log10_TPM_plus1) %>%
    pivot_wider(names_from = gene, values_from = log10_TPM_plus1) %>%
    arrange(sample)
  gene_names <- setdiff(names(wide), "sample")
  rows <- lapply(combn(gene_names, 2, simplify = FALSE), function(z) {
    x <- wide[[z[1]]]; y <- wide[[z[2]]]
    if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) return(NULL)
    test <- suppressWarnings(cor.test(x, y, method = "pearson"))
    data.frame(from = z[1], to = z[2], r = unname(test$estimate), p = test$p.value)
  })
  bind_rows(rows) %>%
    mutate(FDR = p.adjust(p, method = "BH"),
           direction = ifelse(r >= 0, "Positive", "Negative"), abs_r = abs(r),
           habitat = habitat_code)
}

all_edges <- bind_rows(lapply(levels(dat$habitat), cor_edges))

# Keep the original supplied-code criterion. FDR is retained in the output table
# so a stricter multiple-testing filter can be applied without recalculating.
plot_edges <- all_edges %>% filter(abs_r >= 0.60, p < 0.05)

make_panel <- function(habitat_code) {
  panel_nodes <- habitat_node_abundance %>% filter(habitat == habitat_code)
  e <- plot_edges %>% filter(habitat == habitat_code) %>%
    left_join(node_pos %>% select(gene, x, y), by = c("from" = "gene")) %>%
    rename(x_from = x, y_from = y) %>%
    left_join(node_pos %>% select(gene, x, y), by = c("to" = "gene")) %>%
    rename(x_to = x, y_to = y)
  n_samples <- n_distinct(dat$sample[dat$habitat == habitat_code])

  ggplot() +
    geom_segment(data = e,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
          color = direction, linewidth = abs_r), alpha = 0.20, lineend = "round") +
    geom_point(data = panel_nodes,
      aes(x = x, y = y, size = habitat_mean_TPM, fill = module),
      shape = 21, color = "grey15", stroke = 0.35) +
    geom_text(data = panel_nodes,
              aes(x = x, y = y,
                  label = sprintf("%s\n%.2f", gene, habitat_mean_TPM)),
              size = 2.05, vjust = 1.75, lineheight = 0.82,
              check_overlap = FALSE) +
    geom_text(data = centers,
      aes(x = center_x, y = center_y, label = module_labels[as.character(module)]),
      fontface = "bold", size = 4.2) +
    scale_fill_manual(values = module_colors, name = "Functional module") +
    scale_color_manual(values = c(Positive = "#D95F5F", Negative = "#4C78A8"),
                       name = "Pearson direction") +
    scale_linewidth_continuous(range = c(0.25, 1.25), limits = c(0.60, 1),
                               name = "|r|") +
    # Circle area is directly proportional to the untransformed habitat mean TPM.
    scale_size_area(max_size = 16, limits = global_size_limits,
                    name = "Raw mean TPM") +
    coord_fixed(xlim = c(-10.2, 10.2), ylim = c(-8.0, 8.0), clip = "off") +
    labs(title = unname(habitat_labels[habitat_code]),
         subtitle = sprintf("n = %d samples; %d edges (|r| >= 0.60, P < 0.05)",
                            n_samples, nrow(e))) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 8.5, color = "grey35"),
          legend.position = "none")
}

panels <- lapply(levels(dat$habitat), make_panel)
legend_plot <- panels[[1]] + theme(legend.position = "right")
combined <- wrap_plots(panels, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Contig-level CNSAs functional-gene co-occurrence networks",
    subtitle = "Node label and area: untransformed habitat mean TPM; edges: within-habitat Pearson correlations",
    theme = theme(plot.title = element_text(face = "bold", size = 17),
                  plot.subtitle = element_text(size = 10))
  ) & theme(legend.position = "right")

ggsave(file.path(out_dir, "CNSAs_habitat_nested_circle_networks.pdf"), combined,
       width = 16, height = 12, device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "CNSAs_habitat_nested_circle_networks.svg"), combined,
       width = 16, height = 12, device = grDevices::svg, bg = "white")
ggsave(file.path(out_dir, "CNSAs_habitat_nested_circle_networks.png"), combined,
       width = 16, height = 12, dpi = 300, bg = "white")

write.csv(habitat_node_abundance,
          file.path(out_dir, "CNSAs_nested_network_nodes.csv"), row.names = FALSE)
write.csv(all_edges, file.path(out_dir, "CNSAs_all_habitat_pairwise_correlations.csv"), row.names = FALSE)
write.csv(plot_edges, file.path(out_dir, "CNSAs_plotted_habitat_edges.csv"), row.names = FALSE)
write_xlsx(list(nodes_by_habitat = habitat_node_abundance,
                plotted_edges = plot_edges, all_correlations = all_edges),
           file.path(out_dir, "CNSAs_habitat_nested_network_statistics.xlsx"))

message("Finished. Output directory: ", out_dir)
print(plot_edges %>% count(habitat, name = "plotted_edges"))

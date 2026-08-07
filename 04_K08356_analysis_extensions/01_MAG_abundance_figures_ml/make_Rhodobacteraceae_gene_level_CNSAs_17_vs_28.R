suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
source_dir <- file.path(base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected")
out_dir <- file.path(base_dir, "Rhodobacteraceae_gene_level_CNSAs_17_vs_28")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

assignment_file <- file.path(source_dir, "Rhodobacteraceae_MAG_group_assignments_corrected.csv")
ko_file <- "/Users/catherine/Downloads/MAG KEGG/all_bin_kofamscan_filtered_split.20260121200414199.xlsx"

assignments <- read.csv(assignment_file, check.names = FALSE, stringsAsFactors = FALSE)
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
  assignments %>%
    filter(comparison_group == "non-Clade 4") %>%
    arrange(genus, MAG) %>%
    pull(MAG)
)

mandatory_qc <- c(
  nrow(assignments) == 45,
  length(mag_order) == 45,
  length(unique(mag_order)) == 45,
  identical(mag_order[1:17], fixed_clade4),
  all(assignments$family == "f__Rhodobacteraceae"),
  all(mag_order %in% assignments$MAG),
  sum(assignments$comparison_group == "Clade 4 carrier") == 17,
  sum(assignments$comparison_group == "non-Clade 4") == 28
)
if (!all(mandatory_qc)) stop("Mandatory 45 MAG / 17 vs 28 QC failed. No output was generated.")

gene_catalog <- tribble(
  ~category, ~process, ~KO, ~gene_label,
  "Arsenic", "Reduction/detoxification", "K00537", "arsC (glutaredoxin) | K00537",
  "Arsenic", "Efflux", "K01551", "arsA ATPase | K01551",
  "Arsenic", "Efflux", "K03325", "arsB/Acr3 transporter | K03325",
  "Arsenic", "Reduction/detoxification", "K03741", "arsC (thioredoxin) | K03741",
  "Arsenic", "Regulation", "K03892", "arsR regulator | K03892",
  "Arsenic", "Methylation", "K07755", "arsM methyltransferase | K07755",
  "Clade-definition markers", "Aio/Idr-related marker", "K08355", "aioB-related small subunit | K08355",
  "Clade-definition markers", "Aio/Idr-related marker", "K08356", "AioA/IdrA-related DMSOR | K08356",
  "Arsenic", "Resistance", "K11811", "arsH resistance protein | K11811",
  "Arsenic", "Regulation", "K23988", "arsR regulator II | K23988",
  "Arsenic", "Arsenolysis", "K25223", "APGE exporter | K25223",
  "Arsenic", "Arsenolysis", "K25224", "arsenate-transferring GAPDH | K25224",
  "Carbon", "TCA cycle", "K00031", "icd isocitrate dehydrogenase | K00031",
  "Carbon", "Calvin cycle", "K01601", "cbbL RuBisCO large subunit | K01601",
  "Carbon", "Calvin cycle", "K01602", "cbbS RuBisCO small subunit | K01602",
  "Carbon", "TCA cycle", "K01647", "gltA citrate synthase | K01647",
  "Carbon", "Anaplerosis", "K01958", "pyc pyruvate carboxylase | K01958",
  "Carbon", "CO oxidation", "K03518", "coxS CO dehydrogenase | K03518",
  "Carbon", "CO oxidation", "K03519", "coxM CO dehydrogenase | K03519",
  "Carbon", "CO oxidation", "K03520", "coxL CO dehydrogenase | K03520",
  "Carbon", "Methanol oxidation", "K23995", "xoxF methanol dehydrogenase | K23995",
  "Nitrogen", "Denitrification", "K00368", "nirK nitrite reductase | K00368",
  "Nitrogen", "Nitrate reduction", "K00370", "narG/nxrA catalytic subunit | K00370",
  "Nitrogen", "Nitrate reduction", "K00371", "narH/nxrB electron subunit | K00371",
  "Nitrogen", "Assimilatory nitrate reduction", "K00372", "nasA nitrate reductase | K00372",
  "Nitrogen", "Nitrate reduction", "K00374", "narI gamma subunit | K00374",
  "Nitrogen", "Denitrification", "K00376", "nosZ nitrous-oxide reductase | K00376",
  "Nitrogen", "Denitrification", "K02305", "norC nitric-oxide reductase | K02305",
  "Nitrogen", "Denitrification", "K04561", "norB nitric-oxide reductase | K04561",
  "Nitrogen", "Nitrite reduction", "K15864", "NO-forming nitrite reductase | K15864",
  "Nitrogen", "Assimilatory nitrite reduction", "K26138", "nirD-like small subunit | K26138",
  "Nitrogen", "Assimilatory nitrite reduction", "K26139", "nirB-like large subunit | K26139",
  "Sulfur", "Assimilatory sulfate reduction", "K00381", "sir sulfite reductase | K00381",
  "Sulfur", "Sulfate activation", "K00958", "sat sulfate adenylyltransferase | K00958",
  "Sulfur", "Dissimilatory sulfite reduction", "K11180", "dsrA sulfite reductase | K11180",
  "Sulfur", "Dissimilatory sulfite reduction", "K11181", "dsrB sulfite reductase | K11181",
  "Sulfur", "Anaerobic sulfite reduction", "K16951", "asrB sulfite reductase | K16951",
  "Sulfur", "Sulfide oxidation", "K17218", "sqr sulfide:quinone oxidoreductase | K17218",
  "Sulfur", "Thiosulfate oxidation", "K17226", "soxY sulfur carrier | K17226",
  "Sulfur", "Thiosulfate oxidation", "K17227", "soxZ sulfur carrier | K17227",
  "Sulfur", "Sulfide oxidation", "K17229", "fccB sulfide dehydrogenase | K17229",
  "Sulfur", "Sulfide oxidation", "K17230", "fccA cytochrome subunit | K17230",
  "Sulfur", "Dsr-related", "K23077", "Dsr-related protein | K23077",
  "Sulfur", "Dsr-related", "K27196", "Dsr flavoprotein | K27196"
) %>%
  mutate(
    original_order = row_number(),
    category = factor(category, levels = c("Clade-definition markers", "Arsenic", "Carbon", "Nitrogen", "Sulfur"))
  ) %>%
  arrange(category, original_order) %>%
  mutate(
    gene_order = row_number()
  ) %>%
  select(-original_order)

raw_ko <- read_excel(ko_file, col_names = FALSE, .name_repair = "minimal")
names(raw_ko) <- c("MAG", "gene", "KO", "description")
raw_ko <- raw_ko %>% mutate(across(everything(), as.character))

observed_catalog <- raw_ko %>%
  filter(MAG %in% mag_order, KO %in% gene_catalog$KO) %>%
  distinct(KO)
missing_selected_kos <- setdiff(gene_catalog$KO, observed_catalog$KO)
if (length(missing_selected_kos) > 0) {
  stop(paste("Selected KOs absent from the 45 MAG data:", paste(missing_selected_kos, collapse = ", ")))
}

copy_counts <- raw_ko %>%
  filter(MAG %in% mag_order, KO %in% gene_catalog$KO) %>%
  distinct(MAG, gene, KO) %>%
  count(MAG, KO, name = "copy_count")

gene_matrix <- expand_grid(MAG = mag_order, KO = gene_catalog$KO) %>%
  left_join(copy_counts, by = c("MAG", "KO")) %>%
  mutate(copy_count = coalesce(copy_count, 0L)) %>%
  left_join(gene_catalog, by = "KO") %>%
  left_join(assignments %>% select(MAG, comparison_group, K08356_status, source_habitat, genus), by = "MAG") %>%
  mutate(
    MAG = factor(MAG, levels = mag_order),
    gene_label = factor(gene_label, levels = rev(gene_catalog$gene_label)),
    present = copy_count > 0
  )

cliffs_delta <- function(left, right) {
  comparisons <- outer(left, right, FUN = "-")
  (sum(comparisons > 0) - sum(comparisons < 0)) / length(comparisons)
}

gene_statistics <- gene_matrix %>%
  mutate(comparison_group = factor(comparison_group, levels = c("Clade 4 carrier", "non-Clade 4"))) %>%
  group_by(category, process, KO, gene_label, gene_order) %>%
  group_modify(~{
    left <- .x %>% filter(comparison_group == "Clade 4 carrier") %>% pull(copy_count)
    right <- .x %>% filter(comparison_group == "non-Clade 4") %>% pull(copy_count)
    left_present <- sum(left > 0)
    right_present <- sum(right > 0)
    fisher_matrix <- matrix(c(left_present, length(left) - left_present,
                              right_present, length(right) - right_present), nrow = 2, byrow = TRUE)
    fisher_result <- fisher.test(fisher_matrix)
    wilcox_p <- suppressWarnings(wilcox.test(left, right, exact = FALSE)$p.value)
    if (!is.finite(wilcox_p)) wilcox_p <- 1
    tibble(
      clade4_n = length(left), non_clade4_n = length(right),
      clade4_positive = left_present, non_clade4_positive = right_present,
      clade4_prevalence_pct = 100 * left_present / length(left),
      non_clade4_prevalence_pct = 100 * right_present / length(right),
      prevalence_difference_pp = 100 * (left_present / length(left) - right_present / length(right)),
      fisher_odds_ratio = unname(fisher_result$estimate), fisher_p = fisher_result$p.value,
      clade4_mean_copy = mean(left), non_clade4_mean_copy = mean(right),
      mean_copy_difference = mean(left) - mean(right),
      clade4_median_copy = median(left), non_clade4_median_copy = median(right),
      cliffs_delta = cliffs_delta(left, right), wilcox_p = wilcox_p
    )
  }) %>%
  ungroup() %>%
  mutate(
    fisher_fdr = p.adjust(fisher_p, method = "BH"),
    wilcox_fdr = p.adjust(wilcox_p, method = "BH")
  ) %>%
  arrange(gene_order)

is_only_statistics <- gene_matrix %>%
  filter(source_habitat == "IS") %>%
  mutate(comparison_group = factor(comparison_group, levels = c("Clade 4 carrier", "non-Clade 4"))) %>%
  group_by(category, process, KO, gene_label, gene_order) %>%
  group_modify(~{
    left <- .x %>% filter(comparison_group == "Clade 4 carrier") %>% pull(copy_count)
    right <- .x %>% filter(comparison_group == "non-Clade 4") %>% pull(copy_count)
    left_present <- sum(left > 0)
    right_present <- sum(right > 0)
    fisher_matrix <- matrix(c(left_present, length(left) - left_present,
                              right_present, length(right) - right_present), nrow = 2, byrow = TRUE)
    fisher_result <- fisher.test(fisher_matrix)
    wilcox_p <- suppressWarnings(wilcox.test(left, right, exact = FALSE)$p.value)
    if (!is.finite(wilcox_p)) wilcox_p <- 1
    tibble(
      clade4_n = length(left), non_clade4_n = length(right),
      clade4_positive = left_present, non_clade4_positive = right_present,
      clade4_prevalence_pct = 100 * left_present / length(left),
      non_clade4_prevalence_pct = 100 * right_present / length(right),
      prevalence_difference_pp = 100 * (left_present / length(left) - right_present / length(right)),
      fisher_odds_ratio = unname(fisher_result$estimate), fisher_p = fisher_result$p.value,
      clade4_mean_copy = mean(left), non_clade4_mean_copy = mean(right),
      mean_copy_difference = mean(left) - mean(right),
      clade4_median_copy = median(left), non_clade4_median_copy = median(right),
      cliffs_delta = cliffs_delta(left, right), wilcox_p = wilcox_p
    )
  }) %>%
  ungroup() %>%
  mutate(
    fisher_fdr = p.adjust(fisher_p, method = "BH"),
    wilcox_fdr = p.adjust(wilcox_p, method = "BH"),
    sensitivity_set = "IS-only: 17 Clade 4 vs 23 non-Clade 4"
  ) %>%
  arrange(gene_order)

module_catalog <- tribble(
  ~module, ~category, ~required_KOs,
  "Calvin-cycle RuBisCO pair (cbbLS)", "Carbon", "K01601;K01602",
  "Aerobic CO dehydrogenase complex (coxLMS)", "Carbon", "K03518;K03519;K03520",
  "Dissimilatory sulfite reductase pair (dsrAB)", "Sulfur", "K11180;K11181",
  "Nitric-oxide reductase pair (norBC)", "Nitrogen", "K04561;K02305",
  "Sox sulfur-carrier pair (soxYZ)", "Sulfur", "K17226;K17227"
)
module_ko_long <- module_catalog %>%
  separate_rows(required_KOs, sep = ";") %>%
  rename(KO = required_KOs)

module_completeness <- module_ko_long %>%
  left_join(gene_matrix %>%
              transmute(MAG = as.character(MAG), KO, present = copy_count > 0),
            by = "KO") %>%
  group_by(module, category, MAG) %>%
  summarise(n_required = n(), n_present = sum(present),
            module_completeness_pct = 100 * n_present / n_required,
            complete = n_present == n_required, .groups = "drop") %>%
  left_join(assignments %>% select(MAG, comparison_group, source_habitat, genus), by = "MAG")

summarise_module_comparison <- function(data, analysis_set) {
  data %>%
    group_by(module, category) %>%
    group_modify(~{
      left <- .x %>% filter(comparison_group == "Clade 4 carrier") %>% pull(complete)
      right <- .x %>% filter(comparison_group == "non-Clade 4") %>% pull(complete)
      fisher_result <- fisher.test(matrix(c(sum(left), length(left) - sum(left),
                                            sum(right), length(right) - sum(right)),
                                          nrow = 2, byrow = TRUE))
      tibble(
        analysis_set = analysis_set,
        clade4_n = length(left), non_clade4_n = length(right),
        clade4_complete = sum(left), non_clade4_complete = sum(right),
        clade4_complete_pct = 100 * mean(left), non_clade4_complete_pct = 100 * mean(right),
        completeness_difference_pp = 100 * (mean(left) - mean(right)),
        fisher_odds_ratio = unname(fisher_result$estimate), fisher_p = fisher_result$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(fisher_fdr = p.adjust(fisher_p, method = "BH"))
}

module_statistics <- bind_rows(
  summarise_module_comparison(module_completeness, "All source habitats: 17 vs 28"),
  summarise_module_comparison(module_completeness %>% filter(source_habitat == "IS"),
                              "IS-only: 17 vs 23")
)

module_robustness <- module_statistics %>%
  select(module, category, analysis_set, completeness_difference_pp, fisher_fdr) %>%
  pivot_wider(names_from = analysis_set,
              values_from = c(completeness_difference_pp, fisher_fdr),
              names_sep = " | ") %>%
  mutate(
    interpretation = case_when(
      .data[["fisher_fdr | IS-only: 17 vs 23"]] < 0.05 ~ "Robust in IS-only sensitivity analysis",
      .data[["fisher_fdr | All source habitats: 17 vs 28"]] < 0.05 ~
        "Full comparison only; not robust to IS-only habitat sensitivity",
      TRUE ~ "Not significant after BH correction"
    )
  )

write.csv(gene_catalog %>% mutate(category = as.character(category)),
          file.path(out_dir, "selected_CNSAs_KO_catalog.csv"), row.names = FALSE)
write.csv(gene_matrix %>% mutate(MAG = as.character(MAG), gene_label = as.character(gene_label), category = as.character(category)),
          file.path(out_dir, "CNSAs_gene_copy_matrix_45_MAGs_long.csv"), row.names = FALSE)
write.csv(gene_statistics %>% mutate(category = as.character(category), gene_label = as.character(gene_label)),
          file.path(out_dir, "CNSAs_gene_statistics_17_vs_28.csv"), row.names = FALSE)
write.csv(is_only_statistics %>% mutate(category = as.character(category), gene_label = as.character(gene_label)),
          file.path(out_dir, "CNSAs_gene_statistics_IS_only_17_vs_23.csv"), row.names = FALSE)
write.csv(module_completeness,
          file.path(out_dir, "CNS_module_completeness_45_MAGs.csv"), row.names = FALSE)
write.csv(module_statistics,
          file.path(out_dir, "CNS_module_completeness_statistics_full_and_IS_only.csv"), row.names = FALSE)
write.csv(module_robustness,
          file.path(out_dir, "CNS_module_completeness_robustness_interpretation.csv"), row.names = FALSE)
write.csv(data.frame(column_order = seq_along(mag_order), MAG = mag_order,
                     comparison_group = assignments$comparison_group[match(mag_order, assignments$MAG)]),
          file.path(out_dir, "CNSAs_MAG_column_order_45.csv"), row.names = FALSE)

group_colors <- c("Clade 4 carrier" = "#7651A8", "non-Clade 4" = "#2E8B57")
category_colors <- c("Clade-definition markers" = "#6C757D", "Arsenic" = "#B85C38",
                     "Carbon" = "#4C78A8", "Nitrogen" = "#7A5195", "Sulfur" = "#D99B2B")

annotation_df <- assignments %>%
  filter(MAG %in% mag_order) %>%
  mutate(MAG = factor(MAG, levels = mag_order))

p_group <- ggplot(annotation_df, aes(MAG, y = 1, fill = comparison_group)) +
  geom_tile(width = 0.96, height = 0.75) +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 2.2) +
  annotate("text", x = 9, y = 1.72, label = "Clade 4 carriers (n=17)", fontface = "bold", size = 4.1) +
  annotate("text", x = 31.5, y = 1.72, label = "non-Clade 4 Rhodobacteraceae (n=28)", fontface = "bold", size = 4.1) +
  scale_fill_manual(values = group_colors) +
  coord_cartesian(ylim = c(0.55, 1.95), clip = "off") +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(14, 8, 0, 210))

max_copy <- max(gene_matrix$copy_count)
size_breaks <- sort(unique(c(1, 2, 3, max_copy)))
size_breaks <- size_breaks[size_breaks <= max_copy]

p_heat <- ggplot(gene_matrix, aes(MAG, gene_label)) +
  geom_tile(fill = "#FAFAF8", color = "#E7E7E3", linewidth = 0.25) +
  geom_point(data = gene_matrix %>% filter(copy_count > 0),
             aes(size = copy_count, fill = comparison_group),
             shape = 21, color = "white", stroke = 0.25, alpha = 0.95) +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 2.3) +
  scale_fill_manual(values = group_colors, name = "MAG group") +
  scale_size_area(max_size = 6.6, breaks = size_breaks, name = "KO copy count") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = NULL, y = NULL,
    title = "Selected C/N/S/As-related KO profiles across 45 Rhodobacteraceae MAGs",
    subtitle = "Bubble area represents KO copy count; K08355/K08356 are lineage-definition markers rather than independent metabolic enrichments"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.6, color = "#242424"),
    axis.text.y = element_text(size = 8.4, color = "#242424"),
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10.5, color = "#4A4A4A", margin = margin(b = 8)),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(4, 8, 8, 8)
  )

heatmap_figure <- p_group / p_heat + plot_layout(heights = c(0.05, 0.95))

save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, units = "in", bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".svg")), plot, width = width, height = height,
         device = grDevices::svg, units = "in", bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 400, units = "in", bg = "white")
}

save_figure(heatmap_figure, "Rhodobacteraceae_CNSAs_gene_bubble_heatmap_17_vs_28", 19, 18.5)

reference_palette <- c("#482878", "#31688E", "#1F9E89", "#6CCE59", "#FDE725")
reference_breaks <- seq(1, max_copy, by = 1)

reference_habitat_colors <- c("IS" = "#E69F00", "AS" = "#56B4E9")
reference_k08356_colors <- c("K08356-positive" = "#6A3D9A", "K08356-negative" = "#D9D9D9")
reference_annotation_cells <- bind_rows(
  annotation_df %>% transmute(MAG, annotation = "Clade 4 status", y = 3,
                              cell_color = unname(group_colors[comparison_group])),
  annotation_df %>% transmute(MAG, annotation = "K08356 status", y = 2,
                              cell_color = unname(reference_k08356_colors[K08356_status])),
  annotation_df %>% transmute(MAG, annotation = "MAG source habitat", y = 1,
                              cell_color = unname(reference_habitat_colors[source_habitat]))
)

p_reference_top <- ggplot(reference_annotation_cells, aes(MAG, y, fill = cell_color)) +
  geom_tile(width = 0.98, height = 0.88, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 1.8) +
  annotate("text", x = 9, y = 3.92, label = "Clade 4 carriers (n=17)",
           fontface = "bold", size = 4.0, color = group_colors[["Clade 4 carrier"]]) +
  annotate("text", x = 31.5, y = 3.92, label = "non-Clade 4 Rhodobacteraceae (n=28)",
           fontface = "bold", size = 4.0, color = group_colors[["non-Clade 4"]]) +
  scale_fill_identity() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = c(3, 2, 1),
                     labels = c("Clade 4 status", "K08356 status", "MAG source habitat")) +
  coord_cartesian(ylim = c(0.5, 4.18), clip = "off") +
  theme_minimal(base_size = 9.2) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(color = "black", size = 8.4),
        legend.position = "none", plot.margin = margin(15, 8, 0, 8))

category_ranges <- gene_catalog %>%
  mutate(category = as.character(category)) %>%
  group_by(category) %>%
  summarise(first_order = min(gene_order), last_order = max(gene_order), .groups = "drop") %>%
  mutate(
    category = factor(category, levels = names(category_colors)),
    ymin = nrow(gene_catalog) - last_order + 0.5,
    ymax = nrow(gene_catalog) - first_order + 1.5,
    ymid = (ymin + ymax) / 2,
    category_color = unname(category_colors[as.character(category)]),
    display_category = if_else(as.character(category) == "Clade-definition markers", "Clade markers",
                               as.character(category))
  ) %>%
  arrange(category)
category_counts <- gene_catalog %>% count(category, .drop = FALSE) %>% arrange(category)
category_separators <- head(cumsum(rev(category_counts$n)), -1) + 0.5

p_reference_heat <- ggplot(gene_matrix, aes(MAG, gene_label)) +
  geom_tile(fill = "white", color = "#ECECEC", linewidth = 0.22) +
  geom_rect(data = category_ranges,
            aes(xmin = -0.88, xmax = -0.12, ymin = ymin + 0.05, ymax = ymax - 0.05,
                fill = I(category_color)),
            inherit.aes = FALSE, alpha = 0.20, color = NA) +
  geom_point(data = gene_matrix %>% filter(copy_count > 0),
             aes(size = copy_count, color = copy_count), alpha = 0.98) +
  geom_hline(yintercept = category_separators, color = "white", linewidth = 2.3) +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 2.3) +
  geom_text(data = category_ranges,
            aes(x = -0.50, y = ymid, label = display_category),
            inherit.aes = FALSE, angle = 90, fontface = "bold", color = "black", size = 3.2) +
  scale_color_gradientn(
    colors = reference_palette,
    breaks = reference_breaks,
    limits = c(1, max_copy),
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(0.65, 5.4), breaks = reference_breaks,
    limits = c(1, max_copy), guide = "none"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(xlim = c(-0.95, 45.5), clip = "off") +
  labs(
    x = NULL, y = NULL,
    title = NULL, subtitle = NULL
  ) +
  theme_minimal(base_size = 10.2) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5, color = "#202020"),
    axis.text.y = element_text(size = 8.1, color = "#202020"),
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.2, color = "#4A4A4A", margin = margin(b = 8)),
    plot.margin = margin(4, 8, 8, 42)
  )

gradient_steps <- 100
gradient_legend <- tibble(
  ymin = seq(15.8, 18.8, length.out = gradient_steps + 1)[-(gradient_steps + 1)],
  ymax = seq(15.8, 18.8, length.out = gradient_steps + 1)[-1],
  fill_color = colorRampPalette(reference_palette)(gradient_steps)
)
bubble_legend <- tibble(value = reference_breaks, x = seq(0.35, 1.85, length.out = length(reference_breaks)), y = 14.3)
annotation_legend <- bind_rows(
  tibble(section = "Clade 4 status", item = names(group_colors), color = unname(group_colors), y = c(12.55, 11.9)),
  tibble(section = "K08356 status", item = names(reference_k08356_colors), color = unname(reference_k08356_colors), y = c(10.35, 9.7)),
  tibble(section = "MAG source habitat", item = names(reference_habitat_colors), color = unname(reference_habitat_colors), y = c(8.15, 7.5)),
  tibble(section = "Functional category", item = names(category_colors), color = unname(category_colors),
         y = c(6.15, 5.6, 5.05, 4.5, 3.95))
)
section_titles <- tibble(
  section = c("Clade 4 status", "K08356 status", "MAG source habitat", "Functional category"),
  y = c(13.25, 11.05, 8.85, 6.85)
)

p_reference_legend <- ggplot() +
  geom_rect(data = gradient_legend,
            aes(xmin = 0.18, xmax = 0.46, ymin = ymin, ymax = ymax, fill = I(fill_color))) +
  geom_text(data = tibble(value = reference_breaks,
                          y = 15.8 + 3 * (reference_breaks - min(reference_breaks)) /
                            max(1, max(reference_breaks) - min(reference_breaks))),
            aes(x = 0.58, y = y, label = value), hjust = 0, size = 2.8) +
  annotate("text", x = 0.15, y = 19.45, label = "KO copy count", hjust = 0,
           fontface = "bold", size = 3.25) +
  geom_point(data = bubble_legend, aes(x, y, size = value), color = "#1F9E89") +
  geom_text(data = bubble_legend, aes(x, y = 13.7, label = value), size = 2.7) +
  scale_size_continuous(range = c(0.65, 5.4), limits = c(1, max_copy), guide = "none") +
  annotate("text", x = 0.15, y = 15.0, label = "Bubble size", hjust = 0,
           fontface = "bold", size = 3.25) +
  geom_tile(data = annotation_legend, aes(x = 0.28, y = y, fill = I(color)),
            width = 0.25, height = 0.42, color = "#666666", linewidth = 0.2) +
  geom_text(data = annotation_legend, aes(x = 0.48, y = y, label = item),
            hjust = 0, size = 2.65, color = "black") +
  geom_text(data = section_titles, aes(x = 0.15, y = y, label = section),
            hjust = 0, fontface = "bold", size = 3.15) +
  coord_cartesian(xlim = c(0.05, 2.55), ylim = c(3.45, 19.8), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 5, 8, 4))

reference_main <- p_reference_top / p_reference_heat + plot_layout(heights = c(0.105, 0.895))
reference_figure <- (reference_main | p_reference_legend) +
  plot_layout(widths = c(0.90, 0.10)) +
  plot_annotation(
    title = "Selected C/N/S/As-related KO profiles across Rhodobacteraceae MAGs",
    subtitle = "Selected C/N/S/As-related KO profiles; blank cells indicate no detected copy",
    caption = paste0(
      "K08355 and K08356 are Clade-definition / Aio-Idr-related DMSOR markers and are shown only as internal validation. ",
      "All other selected KOs had BH-FDR > 0.05 in both the full 17 vs 28 and IS-only 17 vs 23 comparisons; directional patterns are descriptive trends, not significant enrichment."
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10.2, color = "#4A4A4A"),
                  plot.caption = element_text(size = 8.2, color = "#444444", hjust = 0))
  )
save_figure(reference_figure, "Rhodobacteraceae_CNSAs_reference_style_dot_bubble_17_vs_28", 19, 18.5)

summary_long <- bind_rows(
  gene_statistics %>% transmute(category, gene_label, gene_order,
                                metric = "Prevalence difference (percentage points)",
                                difference = prevalence_difference_pp, fdr = fisher_fdr),
  gene_statistics %>% transmute(category, gene_label, gene_order,
                                metric = "Mean KO copy-number difference (copies per MAG)",
                                difference = mean_copy_difference, fdr = wilcox_fdr)
) %>%
  mutate(
    direction = case_when(
      difference > 0 ~ "Clade 4 carrier",
      difference < 0 ~ "non-Clade 4",
      TRUE ~ "No difference"
    ),
    significance = if_else(fdr < 0.05, "FDR < 0.05", "Not significant"),
    gene_label = factor(gene_label, levels = rev(gene_catalog$gene_label)),
    metric = factor(metric, levels = c("Prevalence difference (percentage points)",
                                      "Mean KO copy-number difference (copies per MAG)"))
  )

p_difference <- ggplot(summary_long, aes(difference, gene_label)) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "#555555") +
  geom_segment(aes(x = 0, xend = difference, yend = gene_label, color = direction), linewidth = 1.2) +
  geom_point(aes(fill = direction, shape = significance), size = 3.2, color = "white", stroke = 0.45) +
  facet_grid(category ~ metric, scales = "free", space = "free_y") +
  scale_color_manual(values = c(group_colors, "No difference" = "#BDBDBD"), guide = "none") +
  scale_fill_manual(values = c(group_colors, "No difference" = "#BDBDBD"), name = "Higher in") +
  scale_shape_manual(values = c("FDR < 0.05" = 24, "Not significant" = 21), name = "Multiple-testing result") +
  labs(
    x = NULL, y = NULL,
    title = "Descriptive KO-level differences between Rhodobacteraceae groups",
    subtitle = "K08355/K08356 are definition markers; all other selected KOs remain non-significant after BH correction"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    axis.text.y = element_text(size = 7.7, color = "#242424"),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.2, color = "#4A4A4A", margin = margin(b = 8)),
    legend.position = "bottom",
    plot.margin = margin(8, 10, 8, 8)
  )

save_figure(p_difference, "Rhodobacteraceae_CNSAs_two_difference_summary_17_vs_28", 16, 18)

full_non_marker_sig <- gene_statistics %>%
  filter(category != "Clade-definition markers", fisher_fdr < 0.05 | wilcox_fdr < 0.05) %>%
  pull(KO)
is_non_marker_sig <- is_only_statistics %>%
  filter(category != "Clade-definition markers", fisher_fdr < 0.05 | wilcox_fdr < 0.05) %>%
  pull(KO)
full_only_modules <- module_robustness %>%
  filter(interpretation == "Full comparison only; not robust to IS-only habitat sensitivity") %>%
  pull(module)

qc_lines <- c(
  "Rhodobacteraceae gene-level C/N/S/As comparison QC",
  "Status: PASSED",
  "",
  sprintf("Figure columns: %d (required 45)", length(mag_order)),
  sprintf("First group: %d Clade 4 carriers (required 17)", sum(assignments$comparison_group == "Clade 4 carrier")),
  sprintf("Second group: %d non-Clade 4 MAGs (required 28)", sum(assignments$comparison_group == "non-Clade 4")),
  sprintf("Unique MAG names: %d", length(unique(mag_order))),
  sprintf("Selected functional genes: %d", nrow(gene_catalog)),
  sprintf("Clade-definition markers: %d", sum(gene_catalog$category == "Clade-definition markers")),
  sprintf("Arsenic genes: %d", sum(gene_catalog$category == "Arsenic")),
  sprintf("Carbon genes: %d", sum(gene_catalog$category == "Carbon")),
  sprintf("Nitrogen genes: %d", sum(gene_catalog$category == "Nitrogen")),
  sprintf("Sulfur genes: %d", sum(gene_catalog$category == "Sulfur")),
  sprintf("Selected KOs absent from 45 MAG data: %d", length(missing_selected_kos)),
  sprintf("IS-only sensitivity comparison: %d Clade 4 vs %d non-Clade 4 MAGs",
          sum(assignments$comparison_group == "Clade 4 carrier" & assignments$source_habitat == "IS"),
          sum(assignments$comparison_group == "non-Clade 4" & assignments$source_habitat == "IS")),
  "Copy number definition: number of distinct gene identifiers annotated to each KO within each MAG.",
  "Prevalence difference: percentage of MAGs with >=1 copy in Clade 4 carriers minus non-Clade 4.",
  "Copy-number difference: mean KO copies per MAG in Clade 4 carriers minus non-Clade 4.",
  "K08355/K08356 are classified as Clade-definition / Aio-Idr-related DMSOR markers, not arsenic functional genes.",
  "Zero-copy cells are blank in both bubble heatmaps.",
  sprintf("Non-marker KOs passing BH-FDR < 0.05 in full comparison: %d", length(full_non_marker_sig)),
  sprintf("Non-marker KOs passing BH-FDR < 0.05 in IS-only comparison: %d", length(is_non_marker_sig)),
  paste0("Module pairs significant only in the full comparison but not IS-only: ",
         if_else(length(full_only_modules) == 0, "none", paste(full_only_modules, collapse = "; "))),
  "Module-pair results are descriptive and must not be called robust enrichment when IS-only sensitivity is non-significant.",
  "MAG completeness and shared-genus/phylogenetic sensitivity are not adjusted in this script and remain required for confirmatory inference.",
  "Existing Figure 4 v3 files were not modified.",
  "",
  "First 17 MAGs:",
  paste(mag_order[1:17], collapse = "\n"),
  "",
  "Last 28 MAGs:",
  paste(mag_order[18:45], collapse = "\n")
)
writeLines(qc_lines, file.path(out_dir, "CNSAs_gene_analysis_QC.txt"))

interpretation_notes <- c(
  "Statistical interpretation for the Rhodobacteraceae KO-profile analysis",
  "",
  "1. K08355 and K08356 are Clade-definition / Aio-Idr-related DMSOR markers. Their significance is internal validation and not independent metabolic enrichment.",
  "2. K08356 is not treated as a canonical arsenite-oxidation marker solely from KO annotation.",
  "3. All non-marker C/N/S/As KOs have BH-FDR >= 0.05 in both the 17 vs 28 and IS-only 17 vs 23 comparisons.",
  "4. Directional differences among non-marker KOs are descriptive tendencies only.",
  "5. Blank heatmap cells indicate zero detected copies; colored bubbles indicate one or more copies.",
  "6. Complete dsrAB and norBC pairs pass BH correction only in the full 17 vs 28 module comparison, but not in IS-only 17 vs 23; these signals are not robust to habitat sensitivity.",
  "7. A single KO does not establish a complete pathway. The supplied module table evaluates only five predefined gene pairs/complexes.",
  "8. MAG completeness and shared-genus/phylogenetic dependence are not controlled here and should be addressed before confirmatory ecological claims."
)
writeLines(interpretation_notes, file.path(out_dir, "STATISTICAL_INTERPRETATION_README.txt"))

cat(paste(qc_lines[1:18], collapse = "\n"), "\n")
cat("\nOutputs written to: ", out_dir, "\n", sep = "")

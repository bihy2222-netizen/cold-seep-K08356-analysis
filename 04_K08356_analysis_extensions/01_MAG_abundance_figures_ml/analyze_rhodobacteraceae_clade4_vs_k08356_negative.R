suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "Rhodobacteraceae_Clade4_vs_K08356_negative_20260805")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

taxonomy_file <- "/Users/catherine/Downloads/MAG 宿主是样品里有哪些 MAG，这些 MAG 的相对丰度么/merged_tpm_matrix 处理哦吼.xlsx"
ko_file <- "/Users/catherine/Downloads/MAG KEGG/all_bin_kofamscan_filtered_split.20260121200414199.xlsx"
module_file <- "/Users/catherine/Downloads/MAG 代谢图/metabolism_summary-师兄给的原始数据.xlsx"
clade_file <- "/Users/catherine/Downloads/MAG 代谢图/K08356_final_clade_mapping_48.csv"

analysis_sheets <- c("MISC", "carbon utilization", "carbon utilization (Woodcroft)",
                     "Transporters", "N", "Organic Nitrogen")
category_rules <- c(
  Carbon = "Methane|Methanogenesis|Acetyl-CoA pathway, CO2|Wood-Ljungdahl|Citrate cycle|Glycolysis|Pyruvate oxidation|TCA /Reductive TCA",
  Nitrogen = "Nitrogen fixation|Nitrification|Denitrification|nitrate reduction|Nitrate assimilation|nitrite \\+ ammonia|Urea cycle|Urea transport",
  Sulfur = "sulfate reduction|Thiosulfate oxidation|Sulfate/thiosulfate transport|Sulfonate transport",
  Arsenic = "^Arsenate$|Arsenic|arsenite|arsenate",
  `Vitamin B12` = "Cobalamin|Vitamin B12|cobinamide",
  Motility = "^Flagellar Assembly$|Chemotaxis"
)

habitat_from_mag <- function(x) case_when(
  str_detect(x, "^SY") ~ "IS",
  str_detect(x, "^(S1_|S2_|S4_|SQ_)") ~ "AS",
  str_detect(x, "^(C1_|C2_|C3_|S13_|S14_|S15_|ES_)") ~ "ES",
  str_detect(x, "^(S3_|NS_|R2111_)") ~ "NS",
  TRUE ~ "Unknown"
)

taxonomy <- read_excel(taxonomy_file, sheet = "Sheet3") %>%
  filter(str_detect(as.character(family), regex("Rhodobacteraceae", ignore_case = TRUE))) %>%
  transmute(MAG = user_genome, family, genus = str_remove(as.character(genus), "^g__"),
            species = str_remove(as.character(species), "^s__")) %>%
  mutate(genus = if_else(is.na(genus) | genus == "", "Unclassified", genus),
         species = if_else(is.na(species) | species == "", "Unclassified", species),
         habitat = habitat_from_mag(MAG))

clades <- read.csv(clade_file, check.names = FALSE)
clade4_mags <- unique(clades$MAG[clades$final_clade == "Clade 4"])
mapped_k08356_mags <- unique(clades$MAG)

raw_ko <- read_excel(ko_file, col_names = FALSE) %>%
  setNames(c("MAG", "gene", "KO", "description")) %>%
  mutate(across(everything(), as.character)) %>%
  filter(!is.na(MAG), !is.na(KO))
raw_rhodo <- raw_ko %>% filter(MAG %in% taxonomy$MAG)
raw_k08356_mags <- unique(raw_rhodo$MAG[raw_rhodo$KO == "K08356"])

taxonomy <- taxonomy %>%
  mutate(group = case_when(
    MAG %in% clade4_mags ~ "Clade 4 carrier",
    !MAG %in% union(mapped_k08356_mags, raw_k08356_mags) ~ "K08356-negative",
    TRUE ~ "Other K08356 clade"
  ))
comparison <- taxonomy %>% filter(group %in% c("Clade 4 carrier", "K08356-negative"))

module_catalog <- bind_rows(lapply(analysis_sheets, function(sheet_name) {
  sheet_data <- read_excel(module_file, sheet = sheet_name, .name_repair = "unique")
  tibble(source_sheet = sheet_name,
         module = replace_na(as.character(sheet_data[[3]]), "Unclassified"),
         gene_id = str_trim(as.character(sheet_data[[1]]))) %>%
    filter(!is.na(gene_id), gene_id != "") %>% distinct()
})) %>% distinct()

ko_counts <- raw_rhodo %>% count(MAG, KO, name = "gene_count")
module_sizes <- module_catalog %>% count(source_sheet, module, name = "gene_set_size")
detected <- module_catalog %>%
  inner_join(ko_counts, by = c("gene_id" = "KO")) %>%
  filter(MAG %in% comparison$MAG, gene_count > 0) %>%
  group_by(source_sheet, module, MAG) %>%
  summarise(genes_detected = n_distinct(gene_id), total_gene_count = sum(gene_count), .groups = "drop")

scores <- crossing(MAG = comparison$MAG, module_sizes) %>%
  left_join(detected, by = c("MAG", "source_sheet", "module")) %>%
  mutate(across(c(genes_detected, total_gene_count), ~replace_na(.x, 0)),
         coverage_pct = 100 * genes_detected / gene_set_size,
         present = genes_detected > 0) %>%
  left_join(comparison, by = "MAG")

stats <- scores %>%
  group_by(source_sheet, module) %>%
  group_modify(~{
    carriers <- filter(.x, group == "Clade 4 carrier")
    negatives <- filter(.x, group == "K08356-negative")
    contingency <- matrix(c(sum(carriers$present), sum(!carriers$present),
                            sum(negatives$present), sum(!negatives$present)), nrow = 2, byrow = TRUE)
    fisher_result <- fisher.test(contingency)
    wilcox_result <- suppressWarnings(wilcox.test(carriers$coverage_pct, negatives$coverage_pct, exact = FALSE))
    tibble(
      n_clade4 = nrow(carriers), n_negative = nrow(negatives),
      clade4_prevalence_pct = 100 * mean(carriers$present),
      negative_prevalence_pct = 100 * mean(negatives$present),
      prevalence_difference_pp = 100 * (mean(carriers$present) - mean(negatives$present)),
      fisher_odds_ratio = unname(fisher_result$estimate), fisher_p = fisher_result$p.value,
      clade4_mean_coverage_pct = mean(carriers$coverage_pct),
      negative_mean_coverage_pct = mean(negatives$coverage_pct),
      mean_coverage_difference_pp = mean(carriers$coverage_pct) - mean(negatives$coverage_pct),
      mann_whitney_p = wilcox_result$p.value
    )
  }) %>% ungroup() %>%
  mutate(fisher_fdr = p.adjust(fisher_p, method = "BH"),
         mann_whitney_fdr = p.adjust(mann_whitney_p, method = "BH"))

focused <- bind_rows(lapply(names(category_rules), function(category_name) {
  stats %>% filter(str_detect(module, regex(category_rules[[category_name]], ignore_case = TRUE))) %>%
    mutate(category = category_name)
})) %>% distinct(source_sheet, module, .keep_all = TRUE) %>%
  filter(clade4_prevalence_pct > 0 | negative_prevalence_pct > 0) %>%
  mutate(rank_score = abs(mean_coverage_difference_pp) + abs(prevalence_difference_pp) / 2)

plot_stats <- focused %>% group_by(category) %>% slice_max(rank_score, n = 4, with_ties = FALSE) %>%
  ungroup() %>% arrange(factor(category, levels = names(category_rules)), desc(mean_coverage_difference_pp)) %>%
  mutate(module_label = str_trunc(module, 58), module_label = factor(module_label, levels = rev(unique(module_label))))

plot_data <- scores %>% inner_join(select(plot_stats, source_sheet, module, module_label, category),
                                   by = c("source_sheet", "module"))
mag_order <- comparison %>%
  mutate(group_order = match(group, c("Clade 4 carrier", "K08356-negative"))) %>%
  arrange(group_order, genus, habitat, MAG) %>% pull(MAG)
plot_data <- plot_data %>%
  mutate(MAG = factor(MAG, levels = mag_order),
         mag_label = factor(paste0(as.character(MAG), "\n", genus),
                            levels = paste0(mag_order, "\n", comparison$genus[match(mag_order, comparison$MAG)])))

group_annotation <- comparison %>% filter(MAG %in% mag_order) %>%
  mutate(MAG = factor(MAG, levels = mag_order), y = 1) %>%
  ggplot(aes(MAG, y, fill = group)) + geom_tile() +
  scale_fill_manual(values = c("Clade 4 carrier" = "#D95F02", "K08356-negative" = "#4C78A8")) +
  theme_void() + theme(legend.position = "top", legend.title = element_blank())

habitat_annotation <- comparison %>% filter(MAG %in% mag_order) %>%
  mutate(MAG = factor(MAG, levels = mag_order), y = 1) %>%
  ggplot(aes(MAG, y, fill = habitat)) + geom_tile() +
  scale_fill_manual(values = c(IS = "#F2A900", AS = "#B77BFF", ES = "#18B7B5", NS = "#18C77A", Unknown = "#999999")) +
  theme_void() + theme(legend.position = "top", legend.title = element_blank())

heatmap_plot <- ggplot(plot_data, aes(mag_label, module_label, fill = coverage_pct)) +
  geom_tile(color = "white", linewidth = 0.18) +
  scale_fill_gradientn(colors = c("#FFFFE5", "#78C679", "#238443", "#004529"), limits = c(0, 100),
                       name = "Gene-set\ncoverage (%)") +
  labs(x = "Rhodobacteraceae MAG (genus shown below)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.2),
        axis.text.y = element_text(size = 7.4), panel.grid = element_blank())

effect_plot <- ggplot(plot_stats, aes(mean_coverage_difference_pp, module_label,
                                      fill = mean_coverage_difference_pp >= 0)) +
  geom_vline(xintercept = 0, linewidth = 0.4) + geom_col(width = 0.72) +
  geom_text(aes(label = ifelse(mann_whitney_fdr < 0.05, "*", "")),
            hjust = ifelse(plot_stats$mean_coverage_difference_pp >= 0, -0.2, 1.2), size = 4) +
  scale_fill_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "#4C78A8"), guide = "none") +
  labs(x = "Mean coverage difference\n(Clade 4 − negative, percentage points)", y = NULL) +
  theme_bw(base_size = 9) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.minor = element_blank())

counts <- table(comparison$group)
final_plot <- (group_annotation / habitat_annotation / (heatmap_plot | effect_plot + plot_layout(widths = c(1, 0.34))) +
                 plot_layout(heights = c(0.55, 0.55, 10))) +
  plot_annotation(
    title = "Metabolic comparison within Rhodobacteraceae: Clade 4 carriers vs K08356-negative MAGs",
    subtitle = sprintf("MAG-equivalent module coverage; n=%d carriers and n=%d negatives. * BH-FDR < 0.05 for Mann–Whitney test.",
                       counts[["Clade 4 carrier"]], counts[["K08356-negative"]]),
    theme = theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 10))
  )

ggsave(file.path(out_dir, "Rhodobacteraceae_Clade4_vs_K08356_negative_heatmap.png"), final_plot,
       width = 18, height = 12, dpi = 450, bg = "white")
ggsave(file.path(out_dir, "Rhodobacteraceae_Clade4_vs_K08356_negative_heatmap.pdf"), final_plot,
       width = 18, height = 12, bg = "white")
ggsave(file.path(out_dir, "Rhodobacteraceae_Clade4_vs_K08356_negative_heatmap.svg"), final_plot,
       width = 18, height = 12, device = grDevices::svg, bg = "white")

write.csv(taxonomy %>% arrange(group, habitat, genus, MAG), file.path(out_dir, "Rhodobacteraceae_MAG_group_assignments.csv"), row.names = FALSE)
write.csv(scores, file.path(out_dir, "Rhodobacteraceae_all_module_scores_by_MAG.csv"), row.names = FALSE)
write.csv(stats %>% arrange(mann_whitney_fdr, fisher_fdr, module), file.path(out_dir, "Rhodobacteraceae_module_difference_statistics.csv"), row.names = FALSE)
write.csv(plot_stats, file.path(out_dir, "Rhodobacteraceae_heatmap_selected_modules.csv"), row.names = FALSE)
write.csv(plot_data, file.path(out_dir, "Rhodobacteraceae_heatmap_long_data.csv"), row.names = FALSE)

summary_table <- tibble(
  metric = c("Rhodobacteraceae MAGs", "Clade 4 carriers", "K08356-negative MAGs", "Other K08356-clade MAGs", "Modules tested", "Focused modules plotted"),
  value = c(nrow(taxonomy), sum(taxonomy$group == "Clade 4 carrier"), sum(taxonomy$group == "K08356-negative"),
            sum(taxonomy$group == "Other K08356 clade"), nrow(stats), nrow(plot_stats))
)
write.csv(summary_table, file.path(out_dir, "analysis_summary.csv"), row.names = FALSE)
print(summary_table)
print(plot_stats %>% select(category, module, mean_coverage_difference_pp, mann_whitney_fdr,
                            prevalence_difference_pp, fisher_fdr) %>% head(15))

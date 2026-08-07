suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

taxonomy_tpm_file <- "/Users/catherine/Downloads/MAG 宿主是样品里有哪些 MAG，这些 MAG 的相对丰度么/merged_tpm_matrix 处理哦吼.xlsx"
abundance_file <- file.path(base_dir, "mag在各个样品中的丰度.xls")
ko_file <- "/Users/catherine/Downloads/MAG KEGG/all_bin_kofamscan_filtered_split.20260121200414199.xlsx"
clade_file <- "/Users/catherine/Downloads/MAG 代谢图/K08356_final_clade_mapping_48.csv"
habitat_file <- "/Users/catherine/Downloads/Rstudio/sample_info.csv"
original_sample_info_file <- "/Users/catherine/Downloads/宏基因组学分析过程文件/冷泉样品信息表-1012重命名版本.xlsx"

fixed_clade4 <- c(
  "SY365BB-8-12_bin8", "SY457BB-8-12_bin4", "SY459WG-0-4_bin37",
  "SY368YW-8-12_bin33", "SY456YB-8-12_bin34", "SY459WG-8-12_bin2",
  "SY365BB-8-12_bin10", "SY366YW-8-12_bin2", "SY456YB-0-4_bin8",
  "SY457BB-8-12_bin15", "SY457BB-4-8_bin29", "SY457BB-8-12_bin25",
  "SY456YB-8-12_bin27", "SY457BB-0-4_bin35", "SY457BB-0-4_bin6",
  "SY365BB-8-12_bin26", "SY457BB-4-8_bin33"
)

taxonomy_tpm <- read_excel(taxonomy_tpm_file, sheet = "Sheet3")
abundance <- read_excel(abundance_file, sheet = "merged_tpm_matrix", .name_repair = "unique")
names(abundance)[1] <- "MAG"
sample_columns <- names(abundance)[-1]
if (length(sample_columns) != 56) stop("The specified MAG TPM workbook must contain 56 sample columns.")

sample_metadata_raw <- read.csv(habitat_file, check.names = FALSE, stringsAsFactors = FALSE)
names(sample_metadata_raw)[names(sample_metadata_raw) == ""] <- "row_id"
sample_metadata <- sample_metadata_raw %>%
  transmute(
    Sample = Sample,
    Site = Site,
    Sample_type = Type,
    Habitat = Group,
    Sample = str_replace(Sample, "^SY457YW", "SY457BB"),
    Sample = str_replace(Sample, "^SY459YW", "SY459WG"),
    Sample = str_replace(Sample, "^SQ_58_-8-12$", "SQ_58_8-12"),
    Sample = str_replace(Sample, "^SQ_81_-8-12$", "SQ_81_8-12")
  )
if (!setequal(sample_metadata$Sample, sample_columns)) {
  stop("Corrected formal sample metadata does not exactly match the 56 TPM sample columns.")
}

rhodo <- taxonomy_tpm %>%
  filter(family == "f__Rhodobacteraceae") %>%
  rename(MAG = user_genome) %>%
  mutate(source_sample = str_remove(MAG, "_bin[^_]+$")) %>%
  left_join(sample_metadata %>% select(source_sample = Sample, source_habitat = Habitat), by = "source_sample")

rhodo_abundance <- abundance %>% filter(MAG %in% rhodo$MAG)

raw_ko <- read_excel(ko_file, col_names = FALSE) %>%
  setNames(c("MAG", "gene", "KO", "description")) %>%
  mutate(across(everything(), as.character))
ko_qc <- raw_ko %>% filter(MAG %in% rhodo$MAG) %>% count(MAG, name = "raw_KO_annotation_rows")
k08356_mags <- raw_ko %>% filter(MAG %in% rhodo$MAG, KO == "K08356") %>% distinct(MAG) %>% pull(MAG)

old_clades <- read.csv(clade_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(MAG %in% rhodo$MAG) %>%
  distinct(MAG, .keep_all = TRUE) %>%
  select(MAG, old_clade = final_clade)

group_assignments <- rhodo %>%
  select(MAG, family, genus, species, source_sample, source_habitat) %>%
  left_join(old_clades, by = "MAG") %>%
  left_join(ko_qc, by = "MAG") %>%
  left_join(rhodo_abundance %>%
              mutate(TPM_data_available = if_all(all_of(sample_columns), ~!is.na(.x))) %>%
              select(MAG, TPM_data_available), by = "MAG") %>%
  mutate(
    K08356_status = if_else(MAG %in% k08356_mags, "K08356-positive", "K08356-negative"),
    clade = case_when(
      MAG %in% fixed_clade4 ~ "Clade 4",
      !is.na(old_clade) ~ old_clade,
      K08356_status == "K08356-positive" ~ "Other/uncertain K08356 clade",
      TRUE ~ "None"
    ),
    comparison_group = if_else(MAG %in% fixed_clade4, "Clade 4 carrier", "non-Clade 4"),
    module_data_available = !is.na(raw_KO_annotation_rows) & raw_KO_annotation_rows > 0,
    TPM_data_available = coalesce(TPM_data_available, FALSE)
  ) %>%
  select(MAG, family, genus, species, source_habitat, K08356_status, clade,
         comparison_group, module_data_available, TPM_data_available,
         source_sample, raw_KO_annotation_rows)

qc <- c(
  Rhodobacteraceae_total = nrow(group_assignments),
  Clade4_carrier = sum(group_assignments$comparison_group == "Clade 4 carrier"),
  non_Clade4 = sum(group_assignments$comparison_group == "non-Clade 4"),
  duplicated_MAG = sum(duplicated(group_assignments$MAG)),
  ungrouped_MAG = sum(is.na(group_assignments$comparison_group) | group_assignments$comparison_group == ""),
  non_Rhodobacteraceae_MAG = sum(group_assignments$family != "f__Rhodobacteraceae"),
  missing_source_habitat = sum(is.na(group_assignments$source_habitat)),
  missing_module_data = sum(!group_assignments$module_data_available),
  missing_TPM_data = sum(!group_assignments$TPM_data_available),
  formal_TPM_samples_matched = length(sample_columns)
)

required_qc <- c(
  qc[["Rhodobacteraceae_total"]] == 45,
  qc[["Clade4_carrier"]] == 17,
  qc[["non_Clade4"]] == 28,
  qc[["duplicated_MAG"]] == 0,
  qc[["ungrouped_MAG"]] == 0,
  qc[["non_Rhodobacteraceae_MAG"]] == 0,
  qc[["missing_source_habitat"]] == 0,
  setequal(group_assignments$MAG[group_assignments$comparison_group == "Clade 4 carrier"], fixed_clade4)
)
if (!all(required_qc)) stop("Mandatory 45 MAG / 17 vs 28 grouping QC failed; no figure may be produced.")

assignments_out <- group_assignments %>%
  mutate(group_order = match(comparison_group, c("Clade 4 carrier", "non-Clade 4")),
         fixed_order = if_else(comparison_group == "Clade 4 carrier", match(MAG, fixed_clade4), Inf)) %>%
  arrange(group_order, fixed_order, genus, MAG) %>%
  select(-group_order, -fixed_order)
write.csv(assignments_out, file.path(out_dir, "Rhodobacteraceae_MAG_group_assignments_corrected.csv"), row.names = FALSE)
writeLines(assignments_out$MAG[assignments_out$comparison_group == "Clade 4 carrier"],
           file.path(out_dir, "Clade4_carrier_17_MAGs.txt"))
writeLines(assignments_out$MAG[assignments_out$comparison_group == "non-Clade 4"],
           file.path(out_dir, "non_Clade4_28_MAGs.txt"))

previously_missing <- assignments_out %>%
  filter(MAG %in% c("SY457BB-0-4_bin35", "SY368YW-4-8_bin24")) %>%
  select(MAG, comparison_group, module_data_available, raw_KO_annotation_rows, TPM_data_available)
write.csv(previously_missing, file.path(out_dir, "previously_missing_MAG_data_QC.csv"), row.names = FALSE)

habitat_counts <- assignments_out %>% count(comparison_group, source_habitat, name = "n_MAG") %>%
  group_by(comparison_group) %>% mutate(within_group_pct = 100 * n_MAG / sum(n_MAG)) %>% ungroup()
write.csv(habitat_counts, file.path(out_dir, "group_habitat_counts_QC.csv"), row.names = FALSE)

report <- c(
  "Figure 4 Rhodobacteraceae corrected grouping QC",
  "Status: PASSED; final Figure 4 generated from corrected 17 vs 28 grouping.",
  "",
  sprintf("Rhodobacteraceae total: %d (required 45)", qc[["Rhodobacteraceae_total"]]),
  sprintf("Clade 4 carrier: %d (required 17)", qc[["Clade4_carrier"]]),
  sprintf("non-Clade 4: %d (required 28)", qc[["non_Clade4"]]),
  sprintf("Duplicated MAG: %d", qc[["duplicated_MAG"]]),
  sprintf("Ungrouped MAG: %d", qc[["ungrouped_MAG"]]),
  sprintf("Non-Rhodobacteraceae MAG: %d", qc[["non_Rhodobacteraceae_MAG"]]),
  sprintf("Missing source habitat: %d", qc[["missing_source_habitat"]]),
  sprintf("Missing raw module annotation: %d", qc[["missing_module_data"]]),
  sprintf("Missing TPM data: %d", qc[["missing_TPM_data"]]),
  sprintf("Formal sample metadata matched: %d samples", qc[["formal_TPM_samples_matched"]]),
  paste0("MAG TPM source: ", abundance_file),
  paste0("Normalized 56-sample metadata source: ", habitat_file),
  paste0("Original sample information source: ", original_sample_info_file),
  "Metadata name reconciliation: SY457YW->SY457BB; SY459YW->SY459WG; SQ_58_-8-12->SQ_58_8-12; SQ_81_-8-12->SQ_81_8-12.",
  "",
  "Previously missing module MAGs:",
  paste(capture.output(print(previously_missing, n = Inf)), collapse = "\n"),
  "",
  "Group x source habitat:",
  paste(capture.output(print(habitat_counts, n = Inf)), collapse = "\n"),
  "",
  "The obsolete 20 vs 23 grouping is not used in this corrected QC."
)
writeLines(report, file.path(out_dir, "Rhodobacteraceae_analysis_QC.txt"))

cat(paste(report, collapse = "\n"), "\n")
cat("\nClade 4 carrier (17):\n", paste(assignments_out$MAG[assignments_out$comparison_group == "Clade 4 carrier"], collapse = "\n"), "\n")
cat("\nnon-Clade 4 (28):\n", paste(assignments_out$MAG[assignments_out$comparison_group == "non-Clade 4"], collapse = "\n"), "\n")

suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

module_file <- "/Users/catherine/Downloads/MAG 代谢图/metabolism_summary-师兄给的原始数据.xlsx"
analysis_sheets <- c("MISC", "carbon utilization", "carbon utilization (Woodcroft)",
                     "Transporters", "N", "Organic Nitrogen")

module_catalog <- bind_rows(lapply(analysis_sheets, function(sheet_name) {
  sheet_data <- read_excel(module_file, sheet = sheet_name, .name_repair = "unique")
  tibble(source_sheet = sheet_name,
         module = coalesce(as.character(sheet_data[[3]]), "Unclassified"),
         gene_id = str_trim(as.character(sheet_data[[1]]))) %>%
    filter(!is.na(gene_id), gene_id != "") %>% distinct()
})) %>% distinct()

ko_counts <- raw_ko %>% filter(MAG %in% assignments_out$MAG) %>% count(MAG, KO, name = "gene_count")
module_sizes <- module_catalog %>% count(source_sheet, module, name = "gene_set_size")
detected <- module_catalog %>%
  inner_join(ko_counts, by = c("gene_id" = "KO"), relationship = "many-to-many") %>%
  group_by(source_sheet, module, MAG) %>%
  summarise(genes_detected = n_distinct(gene_id), total_gene_count = sum(gene_count), .groups = "drop")

module_scores <- crossing(MAG = assignments_out$MAG, module_sizes) %>%
  left_join(detected, by = c("MAG", "source_sheet", "module")) %>%
  mutate(across(c(genes_detected, total_gene_count), ~coalesce(.x, 0)),
         module_coverage_pct = 100 * genes_detected / gene_set_size,
         module_positive = genes_detected > 0) %>%
  left_join(assignments_out %>% select(MAG, comparison_group, source_habitat, genus, species), by = "MAG")

rank_biserial <- function(x, y) {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  test <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  rank_sum <- unname(test$statistic)
  u_stat <- rank_sum - length(x) * (length(x) + 1) / 2
  2 * u_stat / (length(x) * length(y)) - 1
}

module_statistics <- function(data, analysis_label) {
  data %>% group_by(source_sheet, module) %>%
    group_modify(~{
      carriers <- filter(.x, comparison_group == "Clade 4 carrier")
      controls <- filter(.x, comparison_group == "non-Clade 4")
      contingency <- matrix(c(sum(carriers$module_positive), sum(!carriers$module_positive),
                              sum(controls$module_positive), sum(!controls$module_positive)), nrow = 2, byrow = TRUE)
      fisher_result <- fisher.test(contingency)
      mw <- suppressWarnings(wilcox.test(carriers$module_coverage_pct, controls$module_coverage_pct, exact = FALSE))
      tibble(
        analysis = analysis_label,
        n_clade4 = nrow(carriers), n_non_clade4 = nrow(controls),
        clade4_positive_MAG = sum(carriers$module_positive),
        non_clade4_positive_MAG = sum(controls$module_positive),
        clade4_prevalence_pct = 100 * mean(carriers$module_positive),
        non_clade4_prevalence_pct = 100 * mean(controls$module_positive),
        prevalence_difference_pp = 100 * (mean(carriers$module_positive) - mean(controls$module_positive)),
        fisher_odds_ratio = unname(fisher_result$estimate), fisher_p = fisher_result$p.value,
        clade4_mean_coverage_pct = mean(carriers$module_coverage_pct),
        non_clade4_mean_coverage_pct = mean(controls$module_coverage_pct),
        clade4_median_coverage_pct = median(carriers$module_coverage_pct),
        non_clade4_median_coverage_pct = median(controls$module_coverage_pct),
        mean_coverage_difference_pp = mean(carriers$module_coverage_pct) - mean(controls$module_coverage_pct),
        median_coverage_difference_pp = median(carriers$module_coverage_pct) - median(controls$module_coverage_pct),
        mann_whitney_u_p = mw$p.value,
        rank_biserial_effect = rank_biserial(carriers$module_coverage_pct, controls$module_coverage_pct)
      )
    }) %>% ungroup() %>%
    mutate(fisher_BH_FDR = p.adjust(fisher_p, method = "BH"),
           mann_whitney_BH_FDR = p.adjust(mann_whitney_u_p, method = "BH"))
}

main_module_stats <- module_statistics(module_scores, "All Rhodobacteraceae: 17 vs 28")
is_module_stats <- module_statistics(filter(module_scores, source_habitat == "IS"), "IS-only: 17 vs 23")
all_module_stats <- bind_rows(main_module_stats, is_module_stats)
write.csv(all_module_stats, file.path(out_dir, "Rhodobacteraceae_all_module_statistics_17_vs_28.csv"), row.names = FALSE)

mechanism_rules <- tribble(
  ~category, ~pattern,
  "Carbon metabolism", "Citrate cycle|reductive citrate|Reductive acetyl-CoA|Wood-Ljungdahl|Glycolysis|Pyruvate oxidation",
  "Nitrogen metabolism", "Denitrification|Dissimilatory nitrate reduction|nitrite \\+ ammonia|Nitrogen fixation|Nitrate assimilation|Urea transport",
  "Sulfur metabolism", "Dissimilatory sulfate reduction|Assimilatory sulfate reduction|Thiosulfate oxidation by SOX|Sulfate/thiosulfate transport|Sulfonate transport",
  "Arsenic-related", "Arsenate|Arsenic|arsenite|arsenate",
  "Cobalamin/B12", "Cobalamin biosynthesis|Cobalamin salvage|Vitamin B12",
  "Motility/colonization", "Flagellar Assembly|Chemotaxis",
  "Transport systems", "Sulfate/thiosulfate transport|Sulfonate transport|Urea transport|Nitrate/nitrite transport",
  "Redox-related", "Cytochrome c oxidase|aerobic respiration|Anaerobic respiration|DMSO reduction"
)

preselected <- bind_rows(lapply(seq_len(nrow(mechanism_rules)), function(index) {
  main_module_stats %>%
    filter(str_detect(module, regex(mechanism_rules$pattern[index], ignore_case = TRUE))) %>%
    mutate(category = mechanism_rules$category[index], selection_basis = "Predefined mechanism")
})) %>% distinct(source_sheet, module, .keep_all = TRUE) %>%
  filter(clade4_prevalence_pct > 0 | non_clade4_prevalence_pct > 0)

predefined_priority <- c(
  "Incomplete reductive citrate cycle", "Citrate cycle, second carbon oxidation",
  "Reductive acetyl-CoA pathway", "Denitrification, nitrate => nitrogen",
  "Dissimilatory nitrate reduction", "nitrite + ammonia => nitrogen", "Nitrogen fixation",
  "Dissimilatory sulfate reduction", "Assimilatory sulfate reduction",
  "Thiosulfate oxidation by SOX complex", "Sulfate/thiosulfate transport system",
  "Sulfonate transport system", "Arsenate", "Cobalamin biosynthesis",
  "Cobalamin salvage", "Flagellar Assembly", "Chemotaxis"
)

selected_predefined <- bind_rows(lapply(predefined_priority, function(term) {
  preselected %>% filter(str_detect(module, fixed(term, ignore_case = TRUE))) %>%
    arrange(desc(abs(mean_coverage_difference_pp))) %>% slice_head(n = 1)
})) %>% distinct(source_sheet, module, .keep_all = TRUE)

fdr_extras <- main_module_stats %>%
  filter(mann_whitney_BH_FDR < 0.05,
         !paste(source_sheet, module) %in% paste(selected_predefined$source_sheet, selected_predefined$module)) %>%
  arrange(desc(abs(rank_biserial_effect))) %>% slice_head(n = 3) %>%
  mutate(category = "Additional FDR-supported", selection_basis = "Large effect and BH-FDR < 0.05")

selected_stats <- bind_rows(selected_predefined, fdr_extras) %>%
  distinct(source_sheet, module, .keep_all = TRUE) %>%
  left_join(is_module_stats %>%
              select(source_sheet, module, is_mean_coverage_difference_pp = mean_coverage_difference_pp,
                     is_rank_biserial_effect = rank_biserial_effect,
                     is_mann_whitney_BH_FDR = mann_whitney_BH_FDR),
            by = c("source_sheet", "module")) %>%
  mutate(module_label = str_trunc(module, 58))
write.csv(selected_stats, file.path(out_dir, "Rhodobacteraceae_selected_module_statistics_17_vs_28.csv"), row.names = FALSE)

habitat_counts_full <- assignments_out %>%
  count(comparison_group, source_habitat, name = "n_MAG") %>%
  complete(comparison_group, source_habitat = c("IS", "AS", "ES", "NS"), fill = list(n_MAG = 0)) %>%
  group_by(comparison_group) %>% mutate(within_group_pct = 100 * n_MAG / sum(n_MAG)) %>% ungroup()
habitat_table <- xtabs(n_MAG ~ comparison_group + source_habitat, habitat_counts_full)
habitat_fisher <- fisher.test(habitat_table, simulate.p.value = TRUE, B = 99999)
expected <- outer(rowSums(habitat_table), colSums(habitat_table)) / sum(habitat_table)
standardized_residuals <- (habitat_table - expected) / sqrt(expected *
  (1 - rowSums(habitat_table) / sum(habitat_table)) %o%
    (1 - colSums(habitat_table) / sum(habitat_table)))
habitat_stats <- habitat_counts_full %>%
  mutate(expected_count = as.vector(t(expected)), standardized_residual = as.vector(t(standardized_residuals)),
         overall_fisher_p = habitat_fisher$p.value)
write.csv(habitat_stats, file.path(out_dir, "Rhodobacteraceae_source_habitat_statistics.csv"), row.names = FALSE)

tpm_long_mag <- abundance %>%
  filter(MAG %in% assignments_out$MAG) %>%
  pivot_longer(all_of(sample_columns), names_to = "Sample", values_to = "TPM") %>%
  left_join(assignments_out %>% select(MAG, comparison_group), by = "MAG") %>%
  left_join(sample_metadata %>% select(Sample, habitat = Habitat, station = Site, sample_type = Sample_type), by = "Sample")

tpm_by_sample_group <- tpm_long_mag %>%
  group_by(Sample, habitat, station, sample_type, comparison_group) %>%
  summarise(total_TPM = sum(TPM), mean_TPM_per_MAG = mean(TPM),
            prevalent_MAG_count = sum(TPM > 0), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(total_Rhodobacteraceae_TPM = sum(total_TPM),
         fraction_of_total_Rhodobacteraceae_TPM = if_else(total_Rhodobacteraceae_TPM > 0,
                                                           total_TPM / total_Rhodobacteraceae_TPM, 0)) %>%
  ungroup()

tpm_wide <- tpm_by_sample_group %>%
  select(Sample, habitat, station, sample_type, comparison_group, total_TPM, mean_TPM_per_MAG,
         prevalent_MAG_count, fraction_of_total_Rhodobacteraceae_TPM) %>%
  pivot_wider(names_from = comparison_group,
              values_from = c(total_TPM, mean_TPM_per_MAG, prevalent_MAG_count,
                              fraction_of_total_Rhodobacteraceae_TPM), names_sep = "__") %>%
  transmute(Sample, habitat, station, sample_type,
            Clade4_total_TPM = `total_TPM__Clade 4 carrier`,
            nonClade4_total_TPM = `total_TPM__non-Clade 4`,
            Clade4_mean_TPM_per_MAG = `mean_TPM_per_MAG__Clade 4 carrier`,
            nonClade4_mean_TPM_per_MAG = `mean_TPM_per_MAG__non-Clade 4`,
            Clade4_prevalent_MAG_count = `prevalent_MAG_count__Clade 4 carrier`,
            nonClade4_prevalent_MAG_count = `prevalent_MAG_count__non-Clade 4`,
            Clade4_fraction_of_total_Rhodobacteraceae_TPM = `fraction_of_total_Rhodobacteraceae_TPM__Clade 4 carrier`)
write.csv(tpm_wide, file.path(out_dir, "Rhodobacteraceae_TPM_56samples_17_vs_28_long.csv"), row.names = FALSE)

tpm_habitat_summary <- tpm_by_sample_group %>%
  group_by(comparison_group, habitat) %>%
  summarise(n_samples = n(),
            total_TPM_mean = mean(total_TPM), total_TPM_median = median(total_TPM), total_TPM_sd = sd(total_TPM),
            mean_per_MAG_mean = mean(mean_TPM_per_MAG), mean_per_MAG_median = median(mean_TPM_per_MAG),
            fraction_mean = mean(fraction_of_total_Rhodobacteraceae_TPM),
            fraction_median = median(fraction_of_total_Rhodobacteraceae_TPM), .groups = "drop")
write.csv(tpm_habitat_summary, file.path(out_dir, "Rhodobacteraceae_TPM_habitat_summary.csv"), row.names = FALSE)

pairwise_effect <- function(x, y) rank_biserial(x, y)
tpm_metrics <- c("total_TPM", "mean_TPM_per_MAG", "fraction_of_total_Rhodobacteraceae_TPM")
tpm_kw <- bind_rows(lapply(tpm_metrics, function(metric_name) {
  tpm_by_sample_group %>% group_by(comparison_group) %>%
    group_modify(~{
      result <- kruskal.test(.x[[metric_name]] ~ factor(.x$habitat, levels = c("IS", "AS", "ES", "NS")))
      tibble(test = "Kruskal-Wallis among habitats", metric = metric_name,
             contrast = "IS vs AS vs ES vs NS", statistic = unname(result$statistic), p_value = result$p.value,
             effect_size = NA_real_)
    })
}))

habitat_pairs <- combn(c("IS", "AS", "ES", "NS"), 2, simplify = FALSE)
tpm_pairwise <- bind_rows(lapply(tpm_metrics, function(metric_name) {
  bind_rows(lapply(unique(tpm_by_sample_group$comparison_group), function(group_name) {
    bind_rows(lapply(habitat_pairs, function(pair) {
      x <- tpm_by_sample_group %>% filter(comparison_group == group_name, habitat == pair[1]) %>% pull(all_of(metric_name))
      y <- tpm_by_sample_group %>% filter(comparison_group == group_name, habitat == pair[2]) %>% pull(all_of(metric_name))
      result <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
      tibble(comparison_group = group_name, test = "Pairwise Mann-Whitney", metric = metric_name,
             contrast = paste(pair, collapse = " vs "), statistic = unname(result$statistic),
             p_value = result$p.value, effect_size = pairwise_effect(x, y))
    }))
  }))
}))
tpm_statistics <- bind_rows(tpm_kw, tpm_pairwise) %>%
  group_by(test, metric, comparison_group) %>% mutate(BH_FDR = p.adjust(p_value, method = "BH")) %>% ungroup()
write.csv(tpm_statistics, file.path(out_dir, "Rhodobacteraceae_TPM_habitat_statistics.csv"), row.names = FALSE)

group_colors <- c("Clade 4 carrier" = "#D55E00", "non-Clade 4" = "#3B6FB6")
habitat_colors <- c(IS = "#2E8B57", AS = "#8E63CE", ES = "#E3B341", NS = "#3E7CB1")
genus_levels <- sort(unique(assignments_out$genus))
genus_colors <- setNames(hcl.colors(length(genus_levels), "Dynamic"), genus_levels)

mag_order <- assignments_out %>%
  mutate(group_order = match(comparison_group, c("Clade 4 carrier", "non-Clade 4")),
         fixed_order = if_else(comparison_group == "Clade 4 carrier", match(MAG, fixed_clade4), Inf)) %>%
  arrange(group_order, fixed_order, genus, source_habitat, MAG) %>% pull(MAG)
module_order <- selected_stats %>% arrange(match(category, unique(category)), desc(mean_coverage_difference_pp)) %>% pull(module)

axis_qc <- c(
  Figure4A_actual_columns = length(mag_order),
  first17_all_Clade4 = all(mag_order[1:17] %in% fixed_clade4),
  last28_all_non_Clade4 = all(!mag_order[18:45] %in% fixed_clade4),
  all45_Rhodobacteraceae = all(assignments_out$family[match(mag_order, assignments_out$MAG)] == "f__Rhodobacteraceae"),
  duplicated_MAG_names = anyDuplicated(mag_order),
  missing_MAG_names = sum(is.na(mag_order) | mag_order == ""),
  sample_names_mixed_into_axis = sum(mag_order %in% sample_columns)
)
if (!(axis_qc[["Figure4A_actual_columns"]] == 45 &&
      axis_qc[["first17_all_Clade4"]] && axis_qc[["last28_all_non_Clade4"]] &&
      axis_qc[["all45_Rhodobacteraceae"]] && axis_qc[["duplicated_MAG_names"]] == 0 &&
      axis_qc[["missing_MAG_names"]] == 0 && axis_qc[["sample_names_mixed_into_axis"]] == 0)) {
  stop("Figure 4A mandatory 45-MAG axis QC failed.")
}
write.csv(assignments_out %>%
            mutate(column_order = match(MAG, mag_order)) %>% arrange(column_order) %>%
            select(column_order, MAG, comparison_group, K08356_status, source_habitat, genus),
          file.path(out_dir, "Figure4_MAG_column_order_45.csv"), row.names = FALSE)
writeLines(c(
  "Figure 4A axis QC: PASSED",
  sprintf("Figure 4A actual columns: %d", axis_qc[["Figure4A_actual_columns"]]),
  sprintf("First 17 all Clade 4: %s", axis_qc[["first17_all_Clade4"]]),
  sprintf("Last 28 all non-Clade 4: %s", axis_qc[["last28_all_non_Clade4"]]),
  sprintf("All 45 Rhodobacteraceae: %s", axis_qc[["all45_Rhodobacteraceae"]]),
  sprintf("Duplicated MAG names: %d", axis_qc[["duplicated_MAG_names"]]),
  sprintf("Missing MAG names: %d", axis_qc[["missing_MAG_names"]]),
  sprintf("56-sample names mixed into Figure 4A axis: %d", axis_qc[["sample_names_mixed_into_axis"]]),
  "", "First 17 columns:", mag_order[1:17], "", "Last 28 columns:", mag_order[18:45]
), file.path(out_dir, "Figure4_axis_QC.txt"))

heat_data <- module_scores %>%
  inner_join(selected_stats %>% select(source_sheet, module, module_label, category), by = c("source_sheet", "module")) %>%
  mutate(MAG = factor(MAG, levels = mag_order),
         module_label = factor(module_label, levels = rev(selected_stats$module_label)))

annotation_data <- assignments_out %>% filter(MAG %in% mag_order) %>% mutate(MAG = factor(MAG, levels = mag_order), y = 1)
annotation_theme <- theme_void() + theme(legend.position = "top", legend.title = element_blank(),
                                         legend.text = element_text(size = 6.5), legend.key.size = unit(3, "mm"),
                                         axis.title.y = element_text(size = 7, angle = 0, hjust = 1, vjust = 0.5,
                                                                     margin = margin(r = 5)))
p_group_titles <- ggplot() +
  annotate("text", x = 9, y = 1, label = "Clade 4 carriers (n=17)", fontface = "bold", size = 3.2) +
  annotate("text", x = 31.5, y = 1, label = "non-Clade 4 Rhodobacteraceae (n=28)", fontface = "bold", size = 3.2) +
  xlim(0.5, 45.5) + ylim(0.5, 1.5) + theme_void()
p_group_annot <- ggplot(annotation_data, aes(MAG, y, fill = comparison_group)) + geom_tile() +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 1.15) +
  scale_fill_manual(values = group_colors) + labs(y = "Clade status") + annotation_theme + theme(legend.position = "none")
p_k08356_annot <- ggplot(annotation_data, aes(MAG, y, fill = K08356_status)) + geom_tile() +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 1.15) +
  scale_fill_manual(values = c("K08356-positive" = "#6A3D9A", "K08356-negative" = "#D9D9D9")) +
  labs(y = "K08356 status") + annotation_theme + theme(legend.position = "none")
p_habitat_annot <- ggplot(annotation_data, aes(MAG, y, fill = source_habitat)) + geom_tile() +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 1.15) +
  scale_fill_manual(values = habitat_colors, drop = FALSE) + labs(y = "Source habitat") + annotation_theme + theme(legend.position = "none")
p_genus_annot <- ggplot(annotation_data, aes(MAG, y, fill = genus)) + geom_tile() +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 1.15) +
  scale_fill_manual(values = genus_colors) + labs(y = "Genus") + annotation_theme + theme(legend.position = "none")

p_heat <- ggplot(heat_data, aes(MAG, module_label, fill = module_coverage_pct)) +
  geom_tile(color = "white", linewidth = 0.16) +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 1.2) +
  scale_fill_gradientn(colors = c("#F7FCF0", "#C7E9C0", "#41AB5D", "#005A32"), limits = c(0, 100),
                       name = "Module coverage (%)",
                       guide = guide_colorbar(direction = "horizontal", barwidth = unit(40, "mm"),
                                              barheight = unit(3, "mm"), title.position = "top")) +
  scale_x_discrete(drop = FALSE, labels = setNames(mag_order, mag_order)) +
  labs(x = NULL, y = NULL, title = "A  Metabolic-module coverage across 45 Rhodobacteraceae MAGs") + theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.2),
        axis.text.y = element_text(size = 6.7), panel.grid = element_blank(),
        legend.position = "bottom", legend.justification = "center",
        plot.title = element_text(face = "bold", size = 9), plot.margin = margin(0, 3, 2, 2))

p_effect <- selected_stats %>% mutate(module_label = factor(module_label, levels = levels(heat_data$module_label))) %>%
  ggplot(aes(mean_coverage_difference_pp, module_label, fill = mean_coverage_difference_pp >= 0)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey25") + geom_col(width = 0.7) +
  geom_text(aes(label = if_else(mann_whitney_BH_FDR < 0.05, "*", "")),
            hjust = ifelse(selected_stats$mean_coverage_difference_pp >= 0, -0.15, 1.15), size = 3.4) +
  scale_fill_manual(values = c(`TRUE` = group_colors[["Clade 4 carrier"]], `FALSE` = group_colors[["non-Clade 4"]]), guide = "none") +
  labs(x = "Mean coverage difference\n(Clade 4 - non-Clade 4, pp)", y = NULL,
       title = "B  Module effects") +
  theme_bw(base_size = 8) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                  panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 9),
                                  plot.margin = margin(0, 2, 2, 1))

heat_effect_row <- (p_heat | p_effect) + plot_layout(widths = c(4.2, 1.25))
p_ab <- wrap_plots(p_group_titles, p_group_annot, p_k08356_annot, p_habitat_annot, p_genus_annot, heat_effect_row,
                   ncol = 1, heights = c(0.38, 0.27, 0.27, 0.27, 0.27, 8.7))

p_habitat <- ggplot(habitat_counts_full,
                     aes(comparison_group, within_group_pct, fill = factor(source_habitat, levels = c("IS", "AS", "ES", "NS")))) +
  geom_col(width = 0.68, color = "white", linewidth = 0.4) +
  geom_text(aes(label = if_else(n_MAG > 0, sprintf("%d\n(%.1f%%)", n_MAG, within_group_pct), "")),
            position = position_stack(vjust = 0.5), size = 3.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = habitat_colors, drop = FALSE, name = "Source habitat") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "MAG composition (%)", title = "C  MAG source-habitat distribution") +
  theme_classic(base_size = 9) + theme(axis.text.x = element_text(face = "bold"), legend.position = "top")

p_tpm <- tpm_by_sample_group %>%
  mutate(habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS"))) %>%
  ggplot(aes(habitat, log10(total_TPM + 1), fill = comparison_group, color = comparison_group)) +
  geom_boxplot(position = position_dodge(width = 0.72), width = 0.62, outlier.shape = NA, alpha = 0.26, linewidth = 0.55) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72), size = 1.35, alpha = 0.75) +
  scale_fill_manual(values = group_colors) + scale_color_manual(values = group_colors) +
  labs(x = "Sample habitat", y = expression(log[10]("Total MAG TPM" + 1)),
       title = "D  Rhodobacteraceae MAG abundance across 56 samples", fill = NULL, color = NULL) +
  theme_classic(base_size = 9) + theme(legend.position = "top")

bottom_row <- wrap_plots(p_habitat, p_tpm, nrow = 1, widths = c(0.86, 1.5))
figure4 <- wrap_plots(p_ab, bottom_row, ncol = 1, heights = c(3.8, 1.45)) +
  plot_annotation(
    title = "Metabolic and ecological comparison within Rhodobacteraceae: Clade 4 carriers vs non-Clade 4 MAGs",
    subtitle = "MAG-level module coverage (17 vs 28), source habitat, and abundance across 56 metagenomes",
    caption = "Stars indicate Mann-Whitney BH-FDR < 0.05. Module statistics use MAGs as biological units; TPM statistics use samples as biological units. IS-only sensitivity analysis: 17 vs 23 MAGs.",
    theme = theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 10),
                  plot.caption = element_text(size = 7.5, color = "grey30"))
  )

ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28.pdf"), figure4,
       width = 18, height = 15, device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28.png"), figure4,
       width = 18, height = 15, dpi = 450, bg = "white")
ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28.svg"), figure4,
       width = 18, height = 15, device = grDevices::svg, bg = "white")

ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28_v2.pdf"), figure4,
       width = 18, height = 15.5, device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28_v2.png"), figure4,
       width = 18, height = 15.5, dpi = 450, bg = "white")
ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28_v2.svg"), figure4,
       width = 18, height = 15.5, device = grDevices::svg, bg = "white")

p_mean <- tpm_by_sample_group %>% mutate(habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS"))) %>%
  ggplot(aes(habitat, log10(mean_TPM_per_MAG + 1), fill = comparison_group, color = comparison_group)) +
  geom_boxplot(position = position_dodge(width = 0.72), outlier.shape = NA, alpha = 0.25) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72), size = 1.4, alpha = 0.75) +
  scale_fill_manual(values = group_colors) + scale_color_manual(values = group_colors) +
  labs(x = "Sample habitat", y = expression(log[10]("Mean TPM per MAG" + 1)), fill = NULL, color = NULL) + theme_classic(base_size = 10)
p_fraction <- tpm_wide %>% mutate(habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS"))) %>%
  ggplot(aes(habitat, Clade4_fraction_of_total_Rhodobacteraceae_TPM, fill = habitat)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35) + geom_jitter(width = 0.12, size = 1.4, alpha = 0.75) +
  scale_fill_manual(values = habitat_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Sample habitat", y = "Clade 4 fraction of total\nRhodobacteraceae TPM") + theme_classic(base_size = 10)
supp_tpm <- (p_mean | p_fraction) +
  plot_annotation(title = "Supplementary TPM views for Rhodobacteraceae Clade 4 comparison",
                  subtitle = "Per-MAG normalized abundance and Clade 4 contribution to total family-level TPM")
ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28_supp_TPM.svg"), supp_tpm,
       width = 13, height = 5.5, device = grDevices::svg, bg = "white")
ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28_supp_TPM.pdf"), supp_tpm,
       width = 13, height = 5.5, device = cairo_pdf, bg = "white")

final_report <- c(
  readLines(file.path(out_dir, "Rhodobacteraceae_analysis_QC.txt")),
  "", "Final analysis outputs:",
  sprintf("Module scores: %d MAGs x %d modules", n_distinct(module_scores$MAG), n_distinct(paste(module_scores$source_sheet, module_scores$module))),
  sprintf("Main module statistics n: %d vs %d", unique(main_module_stats$n_clade4), unique(main_module_stats$n_non_clade4)),
  sprintf("IS-only module statistics n: %d vs %d", unique(is_module_stats$n_clade4), unique(is_module_stats$n_non_clade4)),
  sprintf("TPM coverage: %d MAGs x %d samples", n_distinct(tpm_long_mag$MAG), n_distinct(tpm_long_mag$Sample)),
  "Within-family comparison reduces family-level taxonomic confounding, but genus-level and habitat-associated effects cannot be completely excluded.",
  "Interpretation: Within Rhodobacteraceae, Clade 4 carriage was associated with differences in selected metabolic-module coverage and with a distinct habitat-abundance pattern in this dataset.",
  "The obsolete 20 vs 23 analysis has been fully replaced by the corrected 17 vs 28 analysis."
)
writeLines(final_report, file.path(out_dir, "Rhodobacteraceae_analysis_QC.txt"))
cat("\nFinal Figure 4 and statistics written to:", out_dir, "\n")

## Figure 4 v3: heatmap-only layout requested after visual review.
## All tiles and annotation rows share one numeric x coordinate system, which
## guarantees exact alignment of the 45 MAG columns.
bubble_group_colors <- c("Clade 4 carrier" = "#7B2CBF", "non-Clade 4" = "#2E8B57")
display_modules <- selected_stats$module_label
n_modules <- length(display_modules)

heatmap_cells <- heat_data %>%
  mutate(x = as.integer(MAG),
         module_order = match(as.character(module_label), display_modules),
         y = n_modules - module_order + 1,
         bubble_color = unname(bubble_group_colors[comparison_group])) %>%
  select(x, y, module_coverage_pct, bubble_color)

annotation_positions <- c("Genus" = n_modules + 1,
                          "MAG source habitat" = n_modules + 2,
                          "K08356 status" = n_modules + 3,
                          "Clade 4 status" = n_modules + 4)
annotation_cells <- bind_rows(
  annotation_data %>% transmute(x = as.integer(MAG), y = annotation_positions[["Clade 4 status"]],
                                cell_color = unname(bubble_group_colors[comparison_group])),
  annotation_data %>% transmute(x = as.integer(MAG), y = annotation_positions[["K08356 status"]],
                                cell_color = if_else(K08356_status == "K08356-positive", "#6A3D9A", "#D9D9D9")),
  annotation_data %>% transmute(x = as.integer(MAG), y = annotation_positions[["MAG source habitat"]],
                                cell_color = unname(habitat_colors[source_habitat])),
  annotation_data %>% transmute(x = as.integer(MAG), y = annotation_positions[["Genus"]],
                                cell_color = unname(genus_colors[genus]))
)

module_y <- n_modules - seq_along(display_modules) + 1
y_breaks <- c(unname(annotation_positions[c("Clade 4 status", "K08356 status", "MAG source habitat", "Genus")]), module_y)
y_labels <- c("Clade 4 status", "K08356 status", "MAG source habitat", "Genus", display_modules)

p_heat_only <- ggplot() +
  geom_tile(data = heatmap_cells, aes(x, y), width = 1, height = 1,
            fill = "white", color = "#E5E5E5", linewidth = 0.16) +
  geom_point(data = heatmap_cells,
             aes(x, y, size = module_coverage_pct, color = bubble_color), alpha = 0.9) +
  scale_size_area(max_size = 8.2, limits = c(0, 100), guide = "none") +
  scale_color_identity() +
  geom_tile(data = annotation_cells,
            aes(x, y, fill = cell_color), width = 1, height = 1,
            color = "white", linewidth = 0.18) +
  scale_fill_identity() +
  geom_vline(xintercept = 17.5, color = "white", linewidth = 1.5) +
  geom_hline(yintercept = n_modules + 0.5, color = "white", linewidth = 1.5) +
  annotate("text", x = 9, y = n_modules + 5.05, label = "Clade 4 carriers (n=17)",
           fontface = "bold", size = 4.1) +
  annotate("text", x = 31.5, y = n_modules + 5.05,
           label = "non-Clade 4 Rhodobacteraceae (n=28)", fontface = "bold", size = 4.1) +
  scale_x_continuous(breaks = seq_along(mag_order), labels = mag_order,
                     limits = c(0.5, 45.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels,
                     limits = c(0.5, n_modules + 5.55), expand = c(0, 0)) +
  labs(x = "Rhodobacteraceae MAG", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.1, color = "black"),
        axis.text.y = element_text(size = 8.2, color = "black"),
        axis.ticks = element_blank(),
        panel.border = element_rect(linewidth = 0.6, color = "black"),
        plot.margin = margin(5, 7, 2, 5))

p_coverage_legend <- ggplot(tibble(value = c(0, 25, 50, 75, 100), x = 1:5, y = 1),
                            aes(x, y, size = value)) +
  geom_point(shape = 21, fill = "#8C8C8C", color = "#4D4D4D", stroke = 0.3) +
  geom_text(aes(y = 0.55, label = value), size = 2.8) +
  scale_size_area(max_size = 8.2, limits = c(0, 100), guide = "none") +
  coord_cartesian(xlim = c(0.5, 5.5), ylim = c(0.35, 1.45), clip = "off") +
  labs(title = "Bubble area: module / gene-set coverage (%)", x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        plot.margin = margin(2, 10, 2, 2))

legend_items <- bind_rows(
  tibble(section = "Clade 4 status", item = c("Clade 4 carrier", "non-Clade 4"),
         color = unname(bubble_group_colors[c("Clade 4 carrier", "non-Clade 4")])),
  tibble(section = "K08356 status", item = c("K08356-positive", "K08356-negative"),
         color = c("#6A3D9A", "#D9D9D9")),
  tibble(section = "MAG source habitat", item = names(habitat_colors), color = unname(habitat_colors)),
  tibble(section = "Genus", item = genus_levels, color = unname(genus_colors[genus_levels]))
) %>%
  group_by(section) %>% mutate(item_index = row_number()) %>% ungroup() %>%
  mutate(display_item = if_else(section == "Genus", str_remove(item, "^g__"), item),
         display_item = if_else(section == "Genus" & (is.na(display_item) | display_item == ""), "Unclassified", display_item),
         section_order = match(section, c("Clade 4 status", "K08356 status", "MAG source habitat", "Genus")),
         x = case_when(
           section_order == 1 ~ item_index,
           section_order == 2 ~ item_index + 3,
           section_order == 3 ~ item_index + 6,
           TRUE ~ ((item_index - 1) %% 8) + 1
         ),
         y = case_when(
           section_order <= 3 ~ 3,
           TRUE ~ 1.65 - floor((item_index - 1) / 8) * 0.65
         ),
         section_x = case_when(section_order == 1 ~ 1, section_order == 2 ~ 4,
                               section_order == 3 ~ 7, TRUE ~ 1),
         label_x = x + 0.16)

section_labels <- legend_items %>% distinct(section, section_order, section_x) %>%
  mutate(y = if_else(section == "Genus", 2.15, 3.52))
p_annotation_legend <- ggplot(legend_items) +
  geom_tile(aes(x, y, fill = color), width = 0.12, height = 0.28, color = "grey45", linewidth = 0.2) +
  geom_text(aes(label_x, y, label = display_item), hjust = 0, size = 2.55) +
  geom_text(data = section_labels, aes(section_x, y, label = section), inherit.aes = FALSE,
            hjust = 0, fontface = "bold", size = 3.1) +
  scale_fill_identity() + coord_cartesian(xlim = c(0.8, 11.4), ylim = c(0.25, 3.8), clip = "off") +
  theme_void() + theme(plot.margin = margin(2, 5, 2, 5))

unified_legend <- wrap_plots(p_coverage_legend, p_annotation_legend, nrow = 1, widths = c(0.9, 3.8))
figure4_v3 <- wrap_plots(p_heat_only, unified_legend, ncol = 1, heights = c(8.8, 1.55)) +
  plot_annotation(
    title = "Metabolic-module profiles of Clade 4 and non-Clade 4 Rhodobacteraceae MAGs",
    subtitle = "Columns are 45 individual MAGs; bubble area represents module coverage",
    caption = "Purple bubbles: 17 Clade 4 carriers; green bubbles: 28 non-Clade 4 MAGs. Clade 4 status, K08356 status, MAG source habitat, and genus are shown as aligned annotation rows.",
    theme = theme(plot.title = element_text(face = "bold", size = 17),
                  plot.subtitle = element_text(size = 10.5),
                  plot.caption = element_text(size = 8.2, color = "grey30"))
  )

ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28_v3_heatmap_only.pdf"), figure4_v3,
       width = 17.5, height = 13.5, device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28_v3_heatmap_only.png"), figure4_v3,
       width = 17.5, height = 13.5, dpi = 450, bg = "white")
ggsave(file.path(out_dir, "Figure4_Rhodobacteraceae_17_vs_28_v3_heatmap_only.svg"), figure4_v3,
       width = 17.5, height = 13.5, device = grDevices::svg, bg = "white")

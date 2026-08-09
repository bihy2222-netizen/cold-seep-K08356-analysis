suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
out_dir <- file.path(base_dir, "Figure4_predefined_module_reanalysis_17_vs_28")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

assignment_file <- file.path(base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected",
                             "Rhodobacteraceae_MAG_group_assignments_corrected.csv")
ko_file <- "/Users/catherine/Downloads/MAG KEGG/all_bin_kofamscan_filtered_split.20260121200414199.xlsx"
module_file <- "/Users/catherine/Downloads/MAG 代谢图/metabolism_summary-师兄给的原始数据.xlsx"

assignments <- read.csv(assignment_file, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(assignments) != 45 ||
    sum(assignments$comparison_group == "Clade 4 carrier") != 17 ||
    sum(assignments$comparison_group == "non-Clade 4") != 28 ||
    sum(assignments$source_habitat == "IS" & assignments$comparison_group == "Clade 4 carrier") != 17 ||
    sum(assignments$source_habitat == "IS" & assignments$comparison_group == "non-Clade 4") != 23) {
  stop("Mandatory 45 MAG / 17 vs 28 / IS-only 17 vs 23 QC failed.")
}

raw_ko <- read_excel(ko_file, col_names = FALSE, .name_repair = "minimal")
names(raw_ko) <- c("MAG", "gene", "KO", "description")
raw_ko <- raw_ko %>%
  mutate(across(everything(), as.character)) %>%
  filter(MAG %in% assignments$MAG) %>%
  distinct(MAG, gene, KO, .keep_all = TRUE)
ko_presence <- raw_ko %>% distinct(MAG, KO) %>% mutate(present = TRUE)

module_sheets <- c("MISC", "carbon utilization", "carbon utilization (Woodcroft)",
                   "Transporters", "N", "Organic Nitrogen")
workbook_catalog <- bind_rows(lapply(module_sheets, function(sheet_name) {
  sheet_data <- read_excel(module_file, sheet = sheet_name, .name_repair = "unique")
  tibble(source_sheet = sheet_name,
         KO = str_trim(as.character(sheet_data[[1]])),
         source_module = as.character(sheet_data[[3]])) %>%
    filter(!is.na(KO), str_detect(KO, "^K[0-9]{5}$"),
           !is.na(source_module), source_module != "") %>%
    distinct()
}))

workbook_kos <- function(pattern) {
  workbook_catalog %>%
    filter(str_detect(source_module, regex(pattern, ignore_case = TRUE))) %>%
    distinct(KO) %>%
    arrange(KO) %>%
    pull(KO)
}

make_steps <- function(module, category, route, step_prefix, alternatives, definition_source) {
  tibble(
    module = module,
    category = category,
    route = route,
    step_id = paste0(step_prefix, seq_along(alternatives)),
    alternative_KOs = alternatives,
    definition_source = definition_source
  )
}

make_ko_set_steps <- function(module, category, kos, definition_source) {
  make_steps(module, category, "main", "step_", kos, definition_source)
}

definitions <- bind_rows(
  make_steps("Denitrification", "Nitrogen", "main", "reaction_", c(
    "K00370|K02567", "K00368|K15864", "K04561", "K00376"
  ), "Curated reaction steps from workbook denitrification KO set"),
  make_steps("DNRA", "Nitrogen", "main", "reaction_", c(
    "K00370|K02567", "K00362|K03385|K15876|K26139"
  ), "Curated nitrate-to-ammonia reaction steps"),
  make_steps("NO reduction", "Nitrogen", "main", "subunit_", c(
    "K04561", "K02305"
  ), "NorBC subunit pair"),
  make_steps("N2O reduction", "Nitrogen", "main", "enzyme_", c(
    "K00376"
  ), "NosZ marker step"),
  make_steps("Assimilatory nitrate reduction", "Nitrogen", "main", "reaction_", c(
    "K00360|K00372|K10534", "K00366|K00367|K26139"
  ), "Curated assimilatory nitrate-to-ammonia reaction steps"),
  make_ko_set_steps("Dissimilatory sulfate/sulfite reduction", "Sulfur",
                    workbook_kos("^Dissimilatory sulfate reduction"),
                    "Original metabolism workbook KO set"),
  make_ko_set_steps("SOX thiosulfate oxidation", "Sulfur",
                    workbook_kos("^Thiosulfate oxidation by SOX complex"),
                    "Original metabolism workbook KO set"),
  make_steps("Sulfide oxidation", "Sulfur", "SQR", "step_", c(
    "K17218"
  ), "Alternative route: SQR"),
  make_steps("Sulfide oxidation", "Sulfur", "FccAB", "subunit_", c(
    "K17229", "K17230"
  ), "Alternative route: flavocytochrome c sulfide dehydrogenase"),
  make_ko_set_steps("B12 biosynthesis", "Cofactor",
                    workbook_kos("^Cobalamin biosynthesis"),
                    "Combined original workbook cobalamin biosynthesis KO sets"),
  make_ko_set_steps("B12 salvage", "Cofactor",
                    workbook_kos("^Cobalamin salvage$"),
                    "Original metabolism workbook KO set"),
  make_steps("Molybdopterin cofactor biosynthesis", "Cofactor", "main", "step_", c(
    "K03639", "K03637", "K03635", "K03636", "K03634", "K03750", "K03638"
  ), "Curated MoaA/MoaC/MoaD/MoaE/MoeB/MoeA/MobA steps"),
  make_ko_set_steps("Flagellar assembly", "Motility",
                    workbook_kos("^Flagellar Assembly$"),
                    "Original metabolism workbook KO set"),
  make_steps("Chemotaxis", "Motility", "main", "core_", c(
    "K03406", "K03407", "K03408", "K03412", "K00575", "K03413"
  ), "Curated MCP/CheA/CheW/CheB/CheR/CheY core"),
  make_steps("Arsenic resistance/detoxification", "Arsenic", "main", "component_", c(
    "K00537|K03741", "K01551", "K03325", "K03892|K23988", "K07755", "K11811", "K25223|K25224"
  ), "Curated resistance/detoxification steps; K08355/K08356 excluded"),
  make_steps("TCA/reductive TCA central-cycle coverage", "Carbon", "main", "reaction_", c(
    "K01647", "K01681|K01682", "K00030|K00031", "K00164|K00174|K00175|K00176|K00177",
    "K01902", "K01903", "K00239|K00240|K00241|K00242|K00244|K00245|K00246|K00247",
    "K01676|K01677|K01678|K01679", "K00024|K00025|K00026"
  ), "Curated central-cycle reaction steps with enzyme alternatives"),
  make_steps("CO oxidation", "Carbon", "main", "subunit_", c(
    "K03518", "K03519", "K03520"
  ), "CoxS/CoxM/CoxL complex"),
  make_ko_set_steps("Carbon fixation (Calvin cycle)", "Carbon",
                    workbook_kos("^Reductive pentose phosphate cycle \\(Calvin cycle\\)$"),
                    "Original metabolism workbook KO set")
)

if (n_distinct(definitions$module) != 17) stop("Exactly 17 predefined modules are required.")
if (any(str_detect(definitions$alternative_KOs, "K08355|K08356"))) {
  stop("K08355/K08356 must not participate in module scoring.")
}
if (any(is.na(definitions$alternative_KOs) | definitions$alternative_KOs == "")) {
  stop("At least one predefined module has no KO definition.")
}

step_alternatives <- definitions %>%
  separate_rows(alternative_KOs, sep = "\\|") %>%
  rename(KO = alternative_KOs)

step_scores <- crossing(MAG = assignments$MAG,
                        definitions %>% select(module, category, route, step_id) %>% distinct()) %>%
  left_join(step_alternatives %>%
              left_join(ko_presence, by = "KO", relationship = "many-to-many") %>%
              group_by(module, category, route, step_id, MAG) %>%
              summarise(step_present = any(coalesce(present, FALSE)), .groups = "drop"),
            by = c("MAG", "module", "category", "route", "step_id")) %>%
  mutate(step_present = coalesce(step_present, FALSE))

route_scores <- step_scores %>%
  group_by(MAG, module, category, route) %>%
  summarise(required_steps = n_distinct(step_id), detected_steps = sum(step_present),
            coverage = detected_steps / required_steps, .groups = "drop")

module_scores <- route_scores %>%
  group_by(MAG, module, category) %>%
  arrange(desc(coverage), desc(required_steps), route) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  rename(best_route = route) %>%
  mutate(
    module_coverage = 100 * coverage,
    module_prevalent = coverage >= 0.75,
    module_complete = coverage == 1
  ) %>%
  left_join(assignments %>% select(MAG, comparison_group, source_habitat, genus), by = "MAG")

cliffs_delta <- function(left, right) {
  comparisons <- outer(left, right, FUN = "-")
  (sum(comparisons > 0) - sum(comparisons < 0)) / length(comparisons)
}

summarise_modules <- function(data, analysis_set) {
  data %>%
    group_by(module, category) %>%
    group_modify(~{
      left <- .x %>% filter(comparison_group == "Clade 4 carrier")
      right <- .x %>% filter(comparison_group == "non-Clade 4")
      if (nrow(left) == 0 || nrow(right) == 0) stop("A comparison group is empty.")
      mw <- suppressWarnings(wilcox.test(left$module_coverage, right$module_coverage, exact = FALSE))
      mw_p <- if_else(is.finite(mw$p.value), mw$p.value, 1)
      contingency <- matrix(c(sum(left$module_prevalent), sum(!left$module_prevalent),
                              sum(right$module_prevalent), sum(!right$module_prevalent)),
                            nrow = 2, byrow = TRUE)
      fisher_result <- fisher.test(contingency)
      tibble(
        analysis_set = analysis_set,
        clade4_n = nrow(left), nonclade4_n = nrow(right),
        clade4_mean_coverage = mean(left$module_coverage),
        nonclade4_mean_coverage = mean(right$module_coverage),
        coverage_difference = mean(left$module_coverage) - mean(right$module_coverage),
        cliffs_delta = cliffs_delta(left$module_coverage, right$module_coverage),
        mann_whitney_p = mw_p,
        clade4_prevalence = mean(left$module_prevalent),
        nonclade4_prevalence = mean(right$module_prevalent),
        prevalence_difference = mean(left$module_prevalent) - mean(right$module_prevalent),
        fisher_odds_ratio = unname(fisher_result$estimate),
        fisher_p = fisher_result$p.value,
        clade4_complete_100pct = mean(left$module_complete),
        nonclade4_complete_100pct = mean(right$module_complete)
      )
    }) %>%
    ungroup() %>%
    mutate(
      coverage_fdr = p.adjust(mann_whitney_p, method = "BH"),
      fisher_fdr = p.adjust(fisher_p, method = "BH"),
      min_fdr = pmin(coverage_fdr, fisher_fdr)
    )
}

main_stats <- summarise_modules(module_scores, "Main: 17 vs 28")
is_stats <- summarise_modules(filter(module_scores, source_habitat == "IS"), "IS-only: 17 vs 23")

definition_summary <- definitions %>%
  separate_rows(alternative_KOs, sep = "\\|") %>%
  group_by(module, category) %>%
  summarise(
    module_size = n_distinct(alternative_KOs),
    required_steps = paste(route, ave(step_id, route, FUN = function(x) length(unique(x))), sep = ":") %>%
      unique() %>% paste(collapse = "; "),
    definition_sources = paste(unique(definition_source), collapse = "; "),
    .groups = "drop"
  )

ranked_table <- main_stats %>%
  select(module, category, clade4_mean_coverage, nonclade4_mean_coverage,
         coverage_difference, cliffs_delta, coverage_fdr,
         clade4_prevalence, nonclade4_prevalence, prevalence_difference, fisher_fdr,
         main_min_fdr = min_fdr) %>%
  left_join(is_stats %>%
              select(module, IS_only_difference = coverage_difference,
                     IS_only_coverage_fdr = coverage_fdr,
                     IS_only_fisher_fdr = fisher_fdr,
                     IS_only_fdr = min_fdr), by = "module") %>%
  left_join(definition_summary %>% select(module, module_size, required_steps), by = "module") %>%
  mutate(
    IS_only_direction_consistent = sign(coverage_difference) == sign(IS_only_difference) |
      coverage_difference == 0 | IS_only_difference == 0,
    main_significant = main_min_fdr < 0.05,
    IS_only_significant = IS_only_fdr < 0.05,
    robust_significant = main_significant & IS_only_significant & IS_only_direction_consistent,
    directional_trend = !main_significant & !IS_only_significant & IS_only_direction_consistent &
      (abs(cliffs_delta) >= 0.33 | abs(coverage_difference) >= 10)
  ) %>%
  arrange(main_min_fdr, desc(IS_only_direction_consistent), desc(abs(cliffs_delta))) %>%
  select(module, category, module_size, required_steps,
         clade4_mean_coverage, nonclade4_mean_coverage, coverage_difference,
         cliffs_delta, coverage_fdr, clade4_prevalence, nonclade4_prevalence,
         prevalence_difference, fisher_fdr, IS_only_difference, IS_only_fdr,
         everything())

full_statistics <- bind_rows(main_stats, is_stats) %>%
  left_join(definition_summary, by = c("module", "category"))

summary_counts <- c(
  main_significant_modules = sum(ranked_table$main_significant),
  is_only_significant_modules = sum(ranked_table$IS_only_significant),
  robust_in_both_modules = sum(ranked_table$robust_significant),
  directional_trend_non_significant_modules = sum(ranked_table$directional_trend)
)

write.csv(definitions, file.path(out_dir, "predefined_module_step_definitions.csv"), row.names = FALSE)
write.csv(step_alternatives, file.path(out_dir, "predefined_module_KO_alternatives_long.csv"), row.names = FALSE)
write.csv(module_scores, file.path(out_dir, "predefined_module_MAG_coverage_45x17.csv"), row.names = FALSE)
write.csv(full_statistics, file.path(out_dir, "predefined_module_statistics_full_and_IS_only.csv"), row.names = FALSE)
write.csv(ranked_table, file.path(out_dir, "predefined_module_statistics_ranked.csv"), row.names = FALSE)

significant_main <- ranked_table %>% filter(main_significant) %>% pull(module)
significant_is <- ranked_table %>% filter(IS_only_significant) %>% pull(module)
robust_modules <- ranked_table %>% filter(robust_significant) %>% pull(module)
trend_modules <- ranked_table %>% filter(directional_trend) %>% pull(module)

summary_report <- c(
  "Predefined Rhodobacteraceae module reanalysis",
  "Status: statistics complete; NO FIGURE GENERATED.",
  "",
  "Scoring rules:",
  "- Coverage is the percentage of predefined biological steps detected in each MAG.",
  "- Alternative KOs within one step are treated as equivalent.",
  "- Alternative sulfide-oxidation routes are scored separately and the higher route coverage is retained.",
  "- Module prevalence is predefined as coverage >= 75%; 100% completion is also retained in the detailed table.",
  "- K08355, K08356 and DIRM-like neighborhood genes do not contribute to any module score.",
  "- Main significance means coverage BH-FDR < 0.05 OR prevalence BH-FDR < 0.05.",
  "- Robust significance additionally requires IS-only significance and consistent coverage direction.",
  "- Directional trend means no FDR significance, consistent IS-only direction, and |Cliff delta| >= 0.33 or |coverage difference| >= 10 percentage points.",
  "",
  sprintf("Predefined modules: %d", n_distinct(definitions$module)),
  sprintf("MAG x module scores: %d (required 45 x 17 = 765)", nrow(module_scores)),
  sprintf("Main significant modules: %d", summary_counts[["main_significant_modules"]]),
  paste0("  ", ifelse(length(significant_main), paste(significant_main, collapse = "; "), "none")),
  sprintf("  Main coverage BH-FDR < 0.05: %d", sum(ranked_table$coverage_fdr < 0.05)),
  sprintf("  Main prevalence Fisher BH-FDR < 0.05: %d", sum(ranked_table$fisher_fdr < 0.05)),
  sprintf("IS-only significant modules: %d", summary_counts[["is_only_significant_modules"]]),
  paste0("  ", ifelse(length(significant_is), paste(significant_is, collapse = "; "), "none")),
  sprintf("  IS-only coverage BH-FDR < 0.05: %d", sum(ranked_table$IS_only_coverage_fdr < 0.05)),
  sprintf("  IS-only prevalence Fisher BH-FDR < 0.05: %d", sum(ranked_table$IS_only_fisher_fdr < 0.05)),
  sprintf("Robust in both analyses: %d", summary_counts[["robust_in_both_modules"]]),
  paste0("  ", ifelse(length(robust_modules), paste(robust_modules, collapse = "; "), "none")),
  sprintf("Directional trends but non-significant: %d", summary_counts[["directional_trend_non_significant_modules"]]),
  paste0("  ", ifelse(length(trend_modules), paste(trend_modules, collapse = "; "), "none")),
  "",
  "Important interpretation:",
  "- These results test predefined module profiles, not single-KO enrichment.",
  "- Workbook-derived KO-set coverage is descriptive and does not by itself prove pathway activity.",
  "- No main figure should be generated until the ranked table and module definitions are reviewed."
)
writeLines(summary_report, file.path(out_dir, "predefined_module_summary_counts.txt"))

result_preview <- ranked_table %>%
  transmute(
    module, category,
    main_coverage_difference_pp = round(coverage_difference, 3),
    main_cliffs_delta = round(cliffs_delta, 3),
    main_coverage_fdr = signif(coverage_fdr, 4),
    main_prevalence_difference = round(prevalence_difference, 3),
    main_fisher_fdr = signif(fisher_fdr, 4),
    IS_only_coverage_difference_pp = round(IS_only_difference, 3),
    IS_only_min_fdr = signif(IS_only_fdr, 4),
    robust_significant, directional_trend
  )
write.table(result_preview, file.path(out_dir, "predefined_module_results_preview.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

qc_report <- c(
  "Predefined module analysis QC",
  "Status: PASSED",
  sprintf("Rhodobacteraceae MAGs: %d", nrow(assignments)),
  sprintf("Clade 4 vs non-Clade 4: %d vs %d",
          sum(assignments$comparison_group == "Clade 4 carrier"),
          sum(assignments$comparison_group == "non-Clade 4")),
  sprintf("IS-only Clade 4 vs non-Clade 4: %d vs %d",
          sum(assignments$comparison_group == "Clade 4 carrier" & assignments$source_habitat == "IS"),
          sum(assignments$comparison_group == "non-Clade 4" & assignments$source_habitat == "IS")),
  sprintf("Predefined modules: %d", n_distinct(definitions$module)),
  sprintf("Module-score rows: %d", nrow(module_scores)),
  sprintf("Duplicate MAG-module rows: %d", sum(duplicated(module_scores[c("MAG", "module")]))),
  sprintf("K08355/K08356 in definitions: %d",
          sum(str_detect(definitions$alternative_KOs, "K08355|K08356"))),
  sprintf("Missing coverage values: %d", sum(is.na(module_scores$module_coverage))),
  sprintf("Coverage outside 0-100: %d",
          sum(module_scores$module_coverage < 0 | module_scores$module_coverage > 100)),
  "No PDF/SVG/PNG files are generated by this script."
)
writeLines(qc_report, file.path(out_dir, "predefined_module_analysis_QC.txt"))

cat(paste(summary_report, collapse = "\n"), "\n")

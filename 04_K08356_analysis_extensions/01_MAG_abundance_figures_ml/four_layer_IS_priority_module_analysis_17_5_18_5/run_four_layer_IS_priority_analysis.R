suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
source_dir <- file.path(base_dir, "Figure4_predefined_module_reanalysis_17_vs_28")
previous_dir <- file.path(base_dir, "K08356_negative_control_module_analysis_17_5_23")
out_dir <- file.path(base_dir, "four_layer_IS_priority_module_analysis_17_5_18_5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

assignments <- read.csv(
  file.path(base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected",
            "Rhodobacteraceae_MAG_group_assignments_corrected.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  mutate(
    four_layer = case_when(
      comparison_group == "Clade 4 carrier" & source_habitat == "IS" ~
        "IS_Clade4_K08356_positive",
      comparison_group == "non-Clade 4" & K08356_status == "K08356-positive" &
        source_habitat == "IS" ~ "IS_other_K08356_positive",
      K08356_status == "K08356-negative" & source_habitat == "IS" ~
        "IS_K08356_negative",
      K08356_status == "K08356-negative" & source_habitat == "AS" ~
        "AS_K08356_negative_reference",
      TRUE ~ NA_character_
    ),
    four_layer = factor(four_layer, levels = c(
      "IS_Clade4_K08356_positive",
      "IS_other_K08356_positive",
      "IS_K08356_negative",
      "AS_K08356_negative_reference"
    ))
  )

coverage <- read.csv(file.path(source_dir, "predefined_module_MAG_coverage_45x17.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE) %>%
  select(-any_of(c("comparison_group", "source_habitat", "genus"))) %>%
  left_join(assignments %>% select(MAG, four_layer, source_habitat, genus, K08356_status,
                                    comparison_group), by = "MAG")
definitions <- read.csv(file.path(source_dir, "predefined_module_step_definitions.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
supplementary_17v23 <- read.csv(
  file.path(previous_dir, "module_statistics_Clade4_17_vs_K08356negative_23.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)

layer_counts <- table(assignments$four_layer)
expected_counts <- c(
  "IS_Clade4_K08356_positive" = 17,
  "IS_other_K08356_positive" = 5,
  "IS_K08356_negative" = 18,
  "AS_K08356_negative_reference" = 5
)
mandatory_qc <- c(
  nrow(assignments) == 45,
  all(as.integer(layer_counts[names(expected_counts)]) == unname(expected_counts)),
  sum(is.na(assignments$four_layer)) == 0,
  sum(duplicated(assignments$MAG)) == 0,
  nrow(coverage) == 45 * 17,
  n_distinct(coverage$module) == 17,
  sum(duplicated(coverage[c("MAG", "module")])) == 0,
  all(assignments$source_habitat[assignments$four_layer == "AS_K08356_negative_reference"] == "AS"),
  all(assignments$K08356_status[assignments$four_layer == "AS_K08356_negative_reference"] ==
        "K08356-negative"),
  all(assignments$source_habitat[assignments$four_layer != "AS_K08356_negative_reference"] == "IS"),
  !any(str_detect(definitions$alternative_KOs, "K08355|K08356"))
)
if (!all(mandatory_qc)) stop("Mandatory four-layer QC failed.")

neutral_module_name <- function(module_name) {
  ifelse(module_name == "Dissimilatory sulfate/sulfite reduction",
         "Dsr/Sat-associated sulfur-redox module", module_name)
}
prepare_output <- function(data) {
  data %>%
    rename(module_definition_name = module) %>%
    mutate(module = neutral_module_name(module_definition_name), .before = module_definition_name)
}

cliffs_delta <- function(left, right) {
  comparisons <- outer(left, right, FUN = "-")
  (sum(comparisons > 0) - sum(comparisons < 0)) / length(comparisons)
}
bootstrap_difference <- function(left, right, iterations = 5000) {
  estimates <- replicate(iterations,
                         mean(sample(left, length(left), replace = TRUE)) -
                           mean(sample(right, length(right), replace = TRUE)))
  c(low = unname(quantile(estimates, 0.025)),
    high = unname(quantile(estimates, 0.975)))
}

compare_layers <- function(data, left_layer, right_layer, analysis_set, seed) {
  set.seed(seed)
  data %>%
    filter(four_layer %in% c(left_layer, right_layer)) %>%
    group_by(module, category) %>%
    group_modify(~{
      left <- .x %>% filter(four_layer == left_layer)
      right <- .x %>% filter(four_layer == right_layer)
      mw <- suppressWarnings(wilcox.test(left$module_coverage, right$module_coverage, exact = FALSE))
      mw_p <- if (is.finite(mw$p.value)) mw$p.value else 1
      fisher_result <- fisher.test(matrix(
        c(sum(left$module_prevalent), sum(!left$module_prevalent),
          sum(right$module_prevalent), sum(!right$module_prevalent)),
        nrow = 2, byrow = TRUE
      ))
      coverage_ci <- bootstrap_difference(left$module_coverage, right$module_coverage)
      prevalence_ci <- 100 * bootstrap_difference(as.numeric(left$module_prevalent),
                                                   as.numeric(right$module_prevalent))
      tibble(
        analysis_set = analysis_set,
        left_group = left_layer, right_group = right_layer,
        left_n = nrow(left), right_n = nrow(right),
        left_mean_coverage = mean(left$module_coverage),
        right_mean_coverage = mean(right$module_coverage),
        coverage_difference = mean(left$module_coverage) - mean(right$module_coverage),
        coverage_bootstrap_95CI_low = coverage_ci[["low"]],
        coverage_bootstrap_95CI_high = coverage_ci[["high"]],
        cliffs_delta = cliffs_delta(left$module_coverage, right$module_coverage),
        mann_whitney_p = mw_p,
        left_prevalence = mean(left$module_prevalent),
        right_prevalence = mean(right$module_prevalent),
        prevalence_difference = mean(left$module_prevalent) - mean(right$module_prevalent),
        prevalence_difference_pp = 100 * prevalence_difference,
        prevalence_bootstrap_95CI_low_pp = prevalence_ci[["low"]],
        prevalence_bootstrap_95CI_high_pp = prevalence_ci[["high"]],
        fisher_odds_ratio = unname(fisher_result$estimate),
        fisher_p = fisher_result$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(
      coverage_BH_FDR = p.adjust(mann_whitney_p, method = "BH"),
      fisher_BH_FDR = p.adjust(fisher_p, method = "BH"),
      min_BH_FDR = pmin(coverage_BH_FDR, fisher_BH_FDR),
      significant = min_BH_FDR < 0.05
    ) %>%
    arrange(min_BH_FDR, desc(abs(cliffs_delta)))
}

primary_17v18 <- compare_layers(
  coverage,
  "IS_Clade4_K08356_positive", "IS_K08356_negative",
  "Primary: IS Clade 4 positive 17 vs IS K08356-negative 18",
  20260808
)
habitat_18v5 <- compare_layers(
  coverage,
  "IS_K08356_negative", "AS_K08356_negative_reference",
  "Habitat sensitivity: IS K08356-negative 18 vs AS K08356-negative 5",
  20260809
)

is_three_group_data <- coverage %>%
  filter(four_layer %in% c("IS_Clade4_K08356_positive",
                           "IS_other_K08356_positive",
                           "IS_K08356_negative")) %>%
  droplevels()

three_group_kw <- is_three_group_data %>%
  group_by(module, category) %>%
  group_modify(~{
    summary <- .x %>%
      group_by(four_layer) %>%
      summarise(n = n(), mean_coverage = mean(module_coverage),
                median_coverage = median(module_coverage),
                prevalence = mean(module_prevalent), .groups = "drop")
    kw <- kruskal.test(module_coverage ~ four_layer, data = .x)
    tibble(
      clade4_n = summary$n[summary$four_layer == "IS_Clade4_K08356_positive"],
      other_positive_n = summary$n[summary$four_layer == "IS_other_K08356_positive"],
      negative_n = summary$n[summary$four_layer == "IS_K08356_negative"],
      clade4_mean_coverage = summary$mean_coverage[summary$four_layer == "IS_Clade4_K08356_positive"],
      other_positive_mean_coverage = summary$mean_coverage[summary$four_layer == "IS_other_K08356_positive"],
      negative_mean_coverage = summary$mean_coverage[summary$four_layer == "IS_K08356_negative"],
      clade4_prevalence = summary$prevalence[summary$four_layer == "IS_Clade4_K08356_positive"],
      other_positive_prevalence = summary$prevalence[summary$four_layer == "IS_other_K08356_positive"],
      negative_prevalence = summary$prevalence[summary$four_layer == "IS_K08356_negative"],
      kruskal_wallis_statistic = unname(kw$statistic),
      kruskal_wallis_p = kw$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(kruskal_wallis_BH_FDR = p.adjust(kruskal_wallis_p, method = "BH")) %>%
  arrange(kruskal_wallis_BH_FDR)

pair_specs <- tribble(
  ~left_layer, ~right_layer, ~contrast, ~seed,
  "IS_Clade4_K08356_positive", "IS_other_K08356_positive",
  "IS Clade 4 vs IS other K08356-positive", 20260810,
  "IS_Clade4_K08356_positive", "IS_K08356_negative",
  "IS Clade 4 vs IS K08356-negative", 20260811,
  "IS_other_K08356_positive", "IS_K08356_negative",
  "IS other K08356-positive vs IS K08356-negative", 20260812
)
three_group_pairwise <- bind_rows(lapply(seq_len(nrow(pair_specs)), function(index) {
  compare_layers(is_three_group_data,
                 pair_specs$left_layer[index], pair_specs$right_layer[index],
                 pair_specs$contrast[index], pair_specs$seed[index])
}))

write.csv(assignments %>% arrange(four_layer, MAG),
          file.path(out_dir, "four_layer_MAG_assignments_17_5_18_5.csv"), row.names = FALSE)
write.csv(prepare_output(primary_17v18),
          file.path(out_dir, "primary_IS_Clade4_17_vs_IS_K08356negative_18.csv"), row.names = FALSE)
write.csv(prepare_output(three_group_kw),
          file.path(out_dir, "IS_three_group_KruskalWallis_17_5_18.csv"), row.names = FALSE)
write.csv(prepare_output(three_group_pairwise),
          file.path(out_dir, "IS_three_group_pairwise_17_5_18.csv"), row.names = FALSE)
write.csv(prepare_output(habitat_18v5),
          file.path(out_dir, "habitat_sensitivity_ISnegative_18_vs_ASnegative_5.csv"), row.names = FALSE)
write.csv(supplementary_17v23,
          file.path(out_dir, "supplementary_uncontrolled_Clade4_17_vs_all_negative_23.csv"), row.names = FALSE)

primary_sig <- primary_17v18 %>% filter(significant) %>% pull(module) %>% neutral_module_name()
habitat_sig <- habitat_18v5 %>% filter(significant) %>% pull(module) %>% neutral_module_name()
kw_sig <- three_group_kw %>% filter(kruskal_wallis_BH_FDR < 0.05) %>%
  pull(module) %>% neutral_module_name()

primary_habitat_audit <- primary_17v18 %>%
  select(module, category,
         primary_difference = coverage_difference,
         primary_coverage_fdr = coverage_BH_FDR,
         primary_fisher_fdr = fisher_BH_FDR) %>%
  left_join(habitat_18v5 %>%
              select(module, habitat_difference = coverage_difference,
                     habitat_coverage_fdr = coverage_BH_FDR,
                     habitat_fisher_fdr = fisher_BH_FDR), by = "module") %>%
  mutate(
    primary_significant = pmin(primary_coverage_fdr, primary_fisher_fdr) < 0.05,
    habitat_significant = pmin(habitat_coverage_fdr, habitat_fisher_fdr) < 0.05,
    habitat_could_conflict = habitat_significant,
    interpretation = case_when(
      primary_significant & !habitat_significant ~ "Primary IS contrast significant; no negative-group habitat FDR signal",
      primary_significant & habitat_significant ~ "Primary significant, but negative-group habitat sensitivity also significant",
      !primary_significant & habitat_significant ~ "Habitat sensitivity only",
      TRUE ~ "Neither comparison significant"
    )
  ) %>%
  arrange(desc(primary_significant), primary_coverage_fdr) %>%
  prepare_output()
write.csv(primary_habitat_audit,
          file.path(out_dir, "primary_vs_habitat_sensitivity_audit.csv"), row.names = FALSE)

qc_lines <- c(
  "Four-layer IS-priority Rhodobacteraceae module analysis",
  "Status: PASSED; statistics generated; NO FIGURE CREATED.",
  "",
  sprintf("IS Clade 4 / K08356-positive: %d", layer_counts[["IS_Clade4_K08356_positive"]]),
  sprintf("IS other K08356-positive: %d", layer_counts[["IS_other_K08356_positive"]]),
  sprintf("IS K08356-negative: %d", layer_counts[["IS_K08356_negative"]]),
  sprintf("AS K08356-negative reference: %d", layer_counts[["AS_K08356_negative_reference"]]),
  sprintf("Total: %d", sum(layer_counts)),
  sprintf("Duplicate or unassigned MAGs: %d",
          sum(duplicated(assignments$MAG)) + sum(is.na(assignments$four_layer))),
  sprintf("Modules reused unchanged: %d", n_distinct(coverage$module)),
  sprintf("K08355/K08356 in module definitions: %d",
          sum(str_detect(definitions$alternative_KOs, "K08355|K08356"))),
  "Bootstrap iterations: 5000 per module and two-group contrast",
  "Prevalence threshold: module coverage >= 75%",
  "",
  sprintf("Primary IS 17 vs 18 significant modules: %d", length(primary_sig)),
  paste0("  ", ifelse(length(primary_sig), paste(primary_sig, collapse = "; "), "none")),
  sprintf("IS three-group Kruskal-Wallis significant modules: %d", length(kw_sig)),
  paste0("  ", ifelse(length(kw_sig), paste(kw_sig, collapse = "; "), "none")),
  sprintf("Negative-group habitat sensitivity 18 vs 5 significant modules: %d", length(habitat_sig)),
  paste0("  ", ifelse(length(habitat_sig), paste(habitat_sig, collapse = "; "), "none")),
  "",
  "Interpretation constraints:",
  "- The AS K08356-negative group is an environmental reference, not a fourth functional type.",
  "- IS Clade 4 vs AS K08356-negative is not tested or interpreted as a Clade 4 effect.",
  "- No Clade4_status x habitat interaction is fitted because AS contains no Clade 4 MAGs.",
  "- The five-member groups are descriptive and directional; no strong inference is assigned.",
  "- The 17 vs 23 analysis is retained only as an uncontrolled-habitat supplementary comparison.",
  "- Existing Figure 4 and earlier statistical tables were not modified."
)
writeLines(qc_lines, file.path(out_dir, "four_layer_IS_priority_QC.txt"))
cat(paste(qc_lines, collapse = "\n"), "\n")

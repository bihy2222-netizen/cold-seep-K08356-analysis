suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

base_dir <- "/Users/catherine/Downloads/MAG-bin-tpm"
source_dir <- file.path(base_dir, "Figure4_predefined_module_reanalysis_17_vs_28")
out_dir <- file.path(base_dir, "K08356_negative_control_module_analysis_17_5_23")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

assignment_file <- file.path(base_dir, "Figure4_Rhodobacteraceae_17_vs_28_corrected",
                             "Rhodobacteraceae_MAG_group_assignments_corrected.csv")
coverage_file <- file.path(source_dir, "predefined_module_MAG_coverage_45x17.csv")
original_stats_file <- file.path(source_dir, "predefined_module_statistics_ranked.csv")
definition_file <- file.path(source_dir, "predefined_module_step_definitions.csv")

assignments <- read.csv(assignment_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  mutate(
    three_group = case_when(
      comparison_group == "Clade 4 carrier" ~ "Clade 4 / K08356-positive",
      comparison_group == "non-Clade 4" & K08356_status == "K08356-positive" ~
        "other non-Clade 4 / K08356-positive",
      K08356_status == "K08356-negative" ~ "K08356-negative Rhodobacteraceae",
      TRUE ~ NA_character_
    ),
    three_group = factor(three_group, levels = c(
      "Clade 4 / K08356-positive",
      "other non-Clade 4 / K08356-positive",
      "K08356-negative Rhodobacteraceae"
    ))
  )
coverage <- read.csv(coverage_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  select(-any_of(c("comparison_group", "source_habitat", "genus"))) %>%
  left_join(assignments %>% select(MAG, three_group, source_habitat, genus, K08356_status,
                                    original_comparison_group = comparison_group), by = "MAG")
original_stats <- read.csv(original_stats_file, check.names = FALSE, stringsAsFactors = FALSE)
definitions <- read.csv(definition_file, check.names = FALSE, stringsAsFactors = FALSE)

group_counts <- table(assignments$three_group)
required_counts <- c(
  "Clade 4 / K08356-positive" = 17,
  "other non-Clade 4 / K08356-positive" = 5,
  "K08356-negative Rhodobacteraceae" = 23
)
is_counts <- table(droplevels(assignments$three_group[assignments$source_habitat == "IS"]))
required_is_counts <- c(
  "Clade 4 / K08356-positive" = 17,
  "other non-Clade 4 / K08356-positive" = 5,
  "K08356-negative Rhodobacteraceae" = 18
)

mandatory_qc <- c(
  nrow(assignments) == 45,
  all(as.integer(group_counts[names(required_counts)]) == unname(required_counts)),
  all(as.integer(is_counts[names(required_is_counts)]) == unname(required_is_counts)),
  sum(duplicated(assignments$MAG)) == 0,
  sum(is.na(assignments$three_group)) == 0,
  all(assignments$K08356_status[assignments$three_group == "K08356-negative Rhodobacteraceae"] ==
        "K08356-negative"),
  all(assignments$K08356_status[assignments$three_group ==
                                  "other non-Clade 4 / K08356-positive"] == "K08356-positive"),
  all(assignments$comparison_group[assignments$three_group ==
                                     "other non-Clade 4 / K08356-positive"] == "non-Clade 4"),
  n_distinct(coverage$module) == 17,
  nrow(coverage) == 45 * 17,
  sum(duplicated(coverage[c("MAG", "module")])) == 0,
  !any(str_detect(definitions$alternative_KOs, "K08355|K08356"))
)
if (!all(mandatory_qc)) stop("Mandatory three-group / module-definition QC failed.")

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

compare_two_groups <- function(data, control_group, analysis_set, seed) {
  set.seed(seed)
  output <- data %>%
    filter(three_group %in% c("Clade 4 / K08356-positive", control_group)) %>%
    group_by(module, category) %>%
    group_modify(~{
      left <- .x %>% filter(three_group == "Clade 4 / K08356-positive")
      right <- .x %>% filter(three_group == control_group)
      mw <- suppressWarnings(wilcox.test(left$module_coverage, right$module_coverage, exact = FALSE))
      mw_p <- if (is.finite(mw$p.value)) mw$p.value else 1
      contingency <- matrix(c(sum(left$module_prevalent), sum(!left$module_prevalent),
                              sum(right$module_prevalent), sum(!right$module_prevalent)),
                            nrow = 2, byrow = TRUE)
      fisher_result <- fisher.test(contingency)
      coverage_ci <- bootstrap_difference(left$module_coverage, right$module_coverage)
      prevalence_ci <- 100 * bootstrap_difference(as.numeric(left$module_prevalent),
                                                   as.numeric(right$module_prevalent))
      tibble(
        analysis_set = analysis_set,
        clade4_n = nrow(left), control_n = nrow(right), control_group = control_group,
        clade4_mean_coverage = mean(left$module_coverage),
        control_mean_coverage = mean(right$module_coverage),
        coverage_difference = mean(left$module_coverage) - mean(right$module_coverage),
        coverage_bootstrap_95CI_low = coverage_ci[["low"]],
        coverage_bootstrap_95CI_high = coverage_ci[["high"]],
        cliffs_delta = cliffs_delta(left$module_coverage, right$module_coverage),
        mann_whitney_p = mw_p,
        clade4_prevalence = mean(left$module_prevalent),
        control_prevalence = mean(right$module_prevalent),
        prevalence_difference = mean(left$module_prevalent) - mean(right$module_prevalent),
        prevalence_difference_pp = 100 * (mean(left$module_prevalent) - mean(right$module_prevalent)),
        prevalence_bootstrap_95CI_low_pp = prevalence_ci[["low"]],
        prevalence_bootstrap_95CI_high_pp = prevalence_ci[["high"]],
        fisher_odds_ratio = unname(fisher_result$estimate),
        fisher_p = fisher_result$p.value,
        clade4_complete_100pct = mean(left$module_complete),
        control_complete_100pct = mean(right$module_complete)
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
  output
}

main_negative <- compare_two_groups(
  coverage,
  "K08356-negative Rhodobacteraceae",
  "Clade 4 vs K08356-negative: 17 vs 23",
  20260806
)
is_negative <- compare_two_groups(
  coverage %>% filter(source_habitat == "IS"),
  "K08356-negative Rhodobacteraceae",
  "IS-only Clade 4 vs IS K08356-negative: 17 vs 18",
  20260807
)

safe_wilcox_p <- function(values, groups, left_group, right_group) {
  left <- values[groups == left_group]
  right <- values[groups == right_group]
  result <- suppressWarnings(wilcox.test(left, right, exact = FALSE))
  if (is.finite(result$p.value)) result$p.value else 1
}

three_group_stats <- coverage %>%
  group_by(module, category) %>%
  group_modify(~{
    group_summary <- .x %>%
      group_by(three_group) %>%
      summarise(n = n(), mean_coverage = mean(module_coverage),
                median_coverage = median(module_coverage),
                prevalence = mean(module_prevalent), .groups = "drop")
    kw <- kruskal.test(module_coverage ~ three_group, data = .x)
    clade_name <- "Clade 4 / K08356-positive"
    other_name <- "other non-Clade 4 / K08356-positive"
    negative_name <- "K08356-negative Rhodobacteraceae"
    tibble(
      clade4_n = group_summary$n[group_summary$three_group == clade_name],
      other_positive_n = group_summary$n[group_summary$three_group == other_name],
      negative_n = group_summary$n[group_summary$three_group == negative_name],
      clade4_mean_coverage = group_summary$mean_coverage[group_summary$three_group == clade_name],
      other_positive_mean_coverage = group_summary$mean_coverage[group_summary$three_group == other_name],
      negative_mean_coverage = group_summary$mean_coverage[group_summary$three_group == negative_name],
      clade4_median_coverage = group_summary$median_coverage[group_summary$three_group == clade_name],
      other_positive_median_coverage = group_summary$median_coverage[group_summary$three_group == other_name],
      negative_median_coverage = group_summary$median_coverage[group_summary$three_group == negative_name],
      clade4_prevalence = group_summary$prevalence[group_summary$three_group == clade_name],
      other_positive_prevalence = group_summary$prevalence[group_summary$three_group == other_name],
      negative_prevalence = group_summary$prevalence[group_summary$three_group == negative_name],
      kruskal_wallis_statistic = unname(kw$statistic),
      kruskal_wallis_p = kw$p.value,
      wilcox_p_clade4_vs_other_positive = safe_wilcox_p(.x$module_coverage, .x$three_group,
                                                        clade_name, other_name),
      wilcox_p_clade4_vs_negative = safe_wilcox_p(.x$module_coverage, .x$three_group,
                                                  clade_name, negative_name),
      wilcox_p_other_positive_vs_negative = safe_wilcox_p(.x$module_coverage, .x$three_group,
                                                           other_name, negative_name)
    )
  }) %>%
  ungroup() %>%
  mutate(
    kruskal_wallis_BH_FDR = p.adjust(kruskal_wallis_p, method = "BH"),
    wilcox_BH_FDR_clade4_vs_other_positive =
      p.adjust(wilcox_p_clade4_vs_other_positive, method = "BH"),
    wilcox_BH_FDR_clade4_vs_negative =
      p.adjust(wilcox_p_clade4_vs_negative, method = "BH"),
    wilcox_BH_FDR_other_positive_vs_negative =
      p.adjust(wilcox_p_other_positive_vs_negative, method = "BH"),
    other_positive_closest_to = case_when(
      abs(other_positive_mean_coverage - clade4_mean_coverage) <
        abs(other_positive_mean_coverage - negative_mean_coverage) ~ "Clade 4 level",
      abs(other_positive_mean_coverage - clade4_mean_coverage) >
        abs(other_positive_mean_coverage - negative_mean_coverage) ~ "K08356-negative level",
      TRUE ~ "Equidistant"
    ),
    other_positive_position = case_when(
      other_positive_mean_coverage > pmax(clade4_mean_coverage, negative_mean_coverage) ~ "Above both",
      other_positive_mean_coverage < pmin(clade4_mean_coverage, negative_mean_coverage) ~ "Below both",
      TRUE ~ "Between Clade 4 and K08356-negative"
    )
  ) %>%
  arrange(kruskal_wallis_BH_FDR)

stable_comparison <- original_stats %>%
  select(module, category,
         original_17v28_difference = coverage_difference,
         original_17v28_coverage_fdr = coverage_fdr,
         original_17v28_fisher_fdr = fisher_fdr) %>%
  left_join(main_negative %>%
              select(module, negative_17v23_difference = coverage_difference,
                     negative_17v23_coverage_fdr = coverage_BH_FDR,
                     negative_17v23_fisher_fdr = fisher_BH_FDR), by = "module") %>%
  left_join(is_negative %>%
              select(module, IS_negative_17v18_difference = coverage_difference,
                     IS_negative_17v18_coverage_fdr = coverage_BH_FDR,
                     IS_negative_17v18_fisher_fdr = fisher_BH_FDR), by = "module") %>%
  mutate(
    direction_consistent_all_three =
      sign(original_17v28_difference) == sign(negative_17v23_difference) &
      sign(original_17v28_difference) == sign(IS_negative_17v18_difference),
    coverage_FDR_significant_all_three =
      original_17v28_coverage_fdr < 0.05 & negative_17v23_coverage_fdr < 0.05 &
      IS_negative_17v18_coverage_fdr < 0.05,
    stable_direction_and_coverage_FDR =
      direction_consistent_all_three & coverage_FDR_significant_all_three
  ) %>%
  arrange(desc(stable_direction_and_coverage_FDR), original_17v28_coverage_fdr)

profile_matrix <- coverage %>%
  select(MAG, module, module_coverage) %>%
  pivot_wider(names_from = module, values_from = module_coverage)
profile_values <- as.matrix(profile_matrix[, -1, drop = FALSE])
scaled_profile <- scale(profile_values)
scaled_profile[is.na(scaled_profile)] <- 0
profile_scaled_long <- as.data.frame(scaled_profile) %>%
  mutate(MAG = profile_matrix$MAG) %>%
  pivot_longer(-MAG, names_to = "module", values_to = "scaled_coverage") %>%
  left_join(assignments %>% select(MAG, three_group), by = "MAG")
centroids <- profile_scaled_long %>%
  group_by(three_group, module) %>%
  summarise(centroid = mean(scaled_coverage), .groups = "drop") %>%
  pivot_wider(names_from = module, values_from = centroid)
centroid_matrix <- as.matrix(centroids[, -1, drop = FALSE])
rownames(centroid_matrix) <- as.character(centroids$three_group)
other_name <- "other non-Clade 4 / K08356-positive"
clade_name <- "Clade 4 / K08356-positive"
negative_name <- "K08356-negative Rhodobacteraceae"
distance_to_clade4 <- sqrt(sum((centroid_matrix[other_name, ] - centroid_matrix[clade_name, ])^2))
distance_to_negative <- sqrt(sum((centroid_matrix[other_name, ] - centroid_matrix[negative_name, ])^2))

neutral_module_name <- function(module_name) {
  ifelse(module_name == "Dissimilatory sulfate/sulfite reduction",
         "Dsr/Sat-associated sulfur-redox module", module_name)
}
prepare_output <- function(data) {
  data %>%
    rename(module_definition_name = module) %>%
    mutate(module = neutral_module_name(module_definition_name), .before = module_definition_name)
}

write.csv(assignments %>% arrange(three_group, MAG),
          file.path(out_dir, "three_group_MAG_assignments.csv"), row.names = FALSE)
write.csv(prepare_output(main_negative),
          file.path(out_dir, "module_statistics_Clade4_17_vs_K08356negative_23.csv"), row.names = FALSE)
write.csv(prepare_output(is_negative),
          file.path(out_dir, "module_statistics_ISonly_17_vs_18.csv"), row.names = FALSE)
write.csv(prepare_output(three_group_stats),
          file.path(out_dir, "module_statistics_three_groups_17_5_23.csv"), row.names = FALSE)
write.csv(prepare_output(stable_comparison),
          file.path(out_dir, "module_stability_across_17v28_17v23_IS17v18.csv"), row.names = FALSE)

other_mags <- assignments %>%
  filter(three_group == other_name) %>%
  select(MAG, source_habitat, genus, species, K08356_status, clade)
write.csv(other_mags, file.path(out_dir, "other_K08356_positive_5_MAGs.csv"), row.names = FALSE)

key_module_levels <- three_group_stats %>%
  filter(module %in% stable_comparison$module[stable_comparison$stable_direction_and_coverage_FDR]) %>%
  select(module, category, clade4_mean_coverage, other_positive_mean_coverage,
         negative_mean_coverage, kruskal_wallis_BH_FDR,
         wilcox_BH_FDR_clade4_vs_other_positive,
         wilcox_BH_FDR_clade4_vs_negative,
         wilcox_BH_FDR_other_positive_vs_negative,
         other_positive_closest_to, other_positive_position) %>%
  prepare_output()
write.csv(key_module_levels,
          file.path(out_dir, "other_K08356_positive_key_module_levels.csv"), row.names = FALSE)

main_sig <- main_negative %>% filter(significant) %>% pull(module) %>% neutral_module_name()
is_sig <- is_negative %>% filter(significant) %>% pull(module) %>% neutral_module_name()
stable_modules <- stable_comparison %>% filter(stable_direction_and_coverage_FDR) %>%
  pull(module) %>% neutral_module_name()
closer_counts <- table(three_group_stats$other_positive_closest_to)
kw_sig <- three_group_stats %>% filter(kruskal_wallis_BH_FDR < 0.05) %>%
  pull(module) %>% neutral_module_name()

qc_lines <- c(
  "K08356-negative control module analysis QC",
  "Status: PASSED; statistics generated; NO FIGURE CREATED.",
  "",
  sprintf("Three groups: %d + %d + %d = %d",
          group_counts[["Clade 4 / K08356-positive"]],
          group_counts[["other non-Clade 4 / K08356-positive"]],
          group_counts[["K08356-negative Rhodobacteraceae"]], sum(group_counts)),
  sprintf("IS-only groups: %d Clade 4 + %d other positive + %d K08356-negative",
          is_counts[["Clade 4 / K08356-positive"]],
          is_counts[["other non-Clade 4 / K08356-positive"]],
          is_counts[["K08356-negative Rhodobacteraceae"]]),
  sprintf("Duplicate MAGs: %d", sum(duplicated(assignments$MAG))),
  sprintf("K08356-positive MAGs in negative group: %d",
          sum(assignments$three_group == negative_name & assignments$K08356_status == "K08356-positive")),
  sprintf("K08356-negative MAGs in other-positive group: %d",
          sum(assignments$three_group == other_name & assignments$K08356_status == "K08356-negative")),
  sprintf("Predefined modules reused unchanged: %d", n_distinct(coverage$module)),
  sprintf("MAG x module rows reused: %d", nrow(coverage)),
  sprintf("K08355/K08356 in module definitions: %d",
          sum(str_detect(definitions$alternative_KOs, "K08355|K08356"))),
  "Bootstrap iterations per module and comparison: 5000",
  "Module prevalence threshold: coverage >= 75%",
  "",
  sprintf("17 vs 23 significant modules: %d", length(main_sig)),
  paste0("  ", ifelse(length(main_sig), paste(main_sig, collapse = "; "), "none")),
  sprintf("IS-only 17 vs 18 significant modules: %d", length(is_sig)),
  paste0("  ", ifelse(length(is_sig), paste(is_sig, collapse = "; "), "none")),
  sprintf("Direction + coverage FDR stable across 17v28, 17v23 and IS17v18: %d",
          length(stable_modules)),
  paste0("  ", ifelse(length(stable_modules), paste(stable_modules, collapse = "; "), "none")),
  "",
  "Other K08356-positive group profile position across 17 modules:",
  sprintf("  Euclidean distance to Clade 4 centroid: %.3f", distance_to_clade4),
  sprintf("  Euclidean distance to K08356-negative centroid: %.3f", distance_to_negative),
  paste0("  Overall closer to: ", ifelse(distance_to_clade4 < distance_to_negative,
                                         "Clade 4", "K08356-negative")),
  paste0("  Module-level closest counts: ",
         paste(names(closer_counts), as.integer(closer_counts), sep = "=", collapse = "; ")),
  sprintf("  Three-group Kruskal-Wallis BH-FDR < 0.05: %d", length(kw_sig)),
  paste0("    ", ifelse(length(kw_sig), paste(kw_sig, collapse = "; "), "none")),
  sprintf("  Pairwise Clade 4 vs other-positive BH-FDR < 0.05: %d",
          sum(three_group_stats$wilcox_BH_FDR_clade4_vs_other_positive < 0.05)),
  sprintf("  Pairwise Clade 4 vs K08356-negative BH-FDR < 0.05: %d",
          sum(three_group_stats$wilcox_BH_FDR_clade4_vs_negative < 0.05)),
  sprintf("  Pairwise other-positive vs K08356-negative BH-FDR < 0.05: %d",
          sum(three_group_stats$wilcox_BH_FDR_other_positive_vs_negative < 0.05)),
  "",
  "The five-member other-positive group is descriptive only; no strong inference is assigned to it.",
  "Existing 17 vs 28 Figure 4 and original statistical tables were not modified."
)
writeLines(qc_lines, file.path(out_dir, "K08356_negative_comparison_QC.txt"))
cat(paste(qc_lines, collapse = "\n"), "\n")

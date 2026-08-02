suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
out_dir <- file.path(script_dir, "results")

score_path <- file.path(
  out_dir, "K08356seq49_MAG48_metabolic_feature_scores.tsv"
)
quality_path <- file.path(out_dir, "K08356seq49_MAG48_quality_audit.tsv")
mapping_path <- file.path(
  out_dir, "K08356seq49_MAG48_clade_habitat_mapping.tsv"
)

required <- c(score_path, quality_path, mapping_path)
if (any(!file.exists(required))) {
  stop(
    "Missing completeness-audit input: ",
    paste(required[!file.exists(required)], collapse = ", ")
  )
}

scores <- read.delim(score_path, check.names = FALSE)
quality <- read.delim(quality_path, check.names = FALSE) %>%
  transmute(
    MAG = MAG_ID,
    quality_habitat = Habitat,
    completeness,
    contamination
  )
mapping <- read.delim(mapping_path, check.names = FALSE) %>%
  transmute(
    MAG,
    contig_gene,
    habitat,
    final_clade,
    final_clade_name
  )

if (nrow(scores) != 1920L ||
    n_distinct(scores$MAG) != 48L ||
    n_distinct(scores$feature_name) != 40L ||
    anyDuplicated(scores[c("MAG", "feature_name")])) {
  stop("Expected 48 MAG x 40 unique feature scores.")
}
if (nrow(quality) != 48L || anyDuplicated(quality$MAG) ||
    !setequal(scores$MAG, quality$MAG)) {
  stop("Quality table does not match the exact 48 scored MAGs.")
}

adjusted_scores <- scores %>%
  left_join(quality, by = "MAG") %>%
  mutate(
    score_semantics = if_else(
      source_type == "Previous algorithm: distinct-gene set coverage",
      "multi-gene partial reconstruction",
      "binary marker carriage"
    ),
    completeness_adjustment_factor = 100 / completeness,
    completeness_adjusted_score_pct = pmin(
      100,
      reconstruction_score_pct * completeness_adjustment_factor
    ),
    raw_ge_50 = reconstruction_score_pct >= 50,
    adjusted_ge_50 = completeness_adjusted_score_pct >= 50,
    threshold_flip_50 = raw_ge_50 != adjusted_ge_50,
    raw_ge_75 = reconstruction_score_pct >= 75,
    adjusted_ge_75 = completeness_adjusted_score_pct >= 75,
    threshold_flip_75 = raw_ge_75 != adjusted_ge_75,
    adjustment_note = if_else(
      score_semantics == "binary marker carriage",
      paste(
        "Presence remains 100 and absence remains 0; completeness cannot",
        "impute an undetected single marker"
      ),
      paste(
        "Sensitivity estimate under random gene loss:",
        "min(100, observed gene-set coverage / MAG completeness fraction)"
      )
    )
  )

safe_spearman <- function(x, y) {
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(c(rho = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(test$estimate), p_value = test$p.value)
}

feature_sensitivity <- adjusted_scores %>%
  group_by(research_category, feature_name, source_type, score_semantics) %>%
  group_modify(function(.x, .y) {
    correlation <- safe_spearman(
      .x$reconstruction_score_pct,
      .x$completeness
    )
    tibble(
      n_MAG = n_distinct(.x$MAG),
      mean_raw_score_pct = mean(.x$reconstruction_score_pct),
      mean_completeness_adjusted_score_pct =
        mean(.x$completeness_adjusted_score_pct),
      mean_score_delta_pct = mean(
        .x$completeness_adjusted_score_pct - .x$reconstruction_score_pct
      ),
      n_MAG_score_changed = sum(
        abs(.x$completeness_adjusted_score_pct -
              .x$reconstruction_score_pct) > 1e-10
      ),
      n_MAG_threshold_flip_50 = sum(.x$threshold_flip_50),
      n_MAG_threshold_flip_75 = sum(.x$threshold_flip_75),
      spearman_score_vs_completeness_rho = correlation[["rho"]],
      spearman_score_vs_completeness_p = correlation[["p_value"]]
    )
  }) %>%
  ungroup() %>%
  mutate(
    spearman_score_vs_completeness_BH = p.adjust(
      spearman_score_vs_completeness_p,
      method = "BH"
    )
  )

membership_adjusted <- adjusted_scores %>%
  inner_join(mapping, by = "MAG", relationship = "many-to-many") %>%
  group_by(
    research_category, feature_name, source_type, score_semantics,
    habitat, final_clade, final_clade_name, MAG
  ) %>%
  summarise(
    n_sequence_memberships = n_distinct(contig_gene),
    completeness = first(completeness),
    raw_score_pct = max(reconstruction_score_pct),
    completeness_adjusted_score_pct =
      max(completeness_adjusted_score_pct),
    raw_ge_50 = any(raw_ge_50),
    adjusted_ge_50 = any(adjusted_ge_50),
    raw_ge_75 = any(raw_ge_75),
    adjusted_ge_75 = any(adjusted_ge_75),
    .groups = "drop"
  )

clade_habitat_sensitivity <- membership_adjusted %>%
  group_by(
    research_category, feature_name, source_type, score_semantics,
    habitat, final_clade, final_clade_name
  ) %>%
  summarise(
    n_unique_MAG = n_distinct(MAG),
    completeness_mean = mean(completeness),
    completeness_median = median(completeness),
    mean_raw_score_pct = mean(raw_score_pct),
    mean_completeness_adjusted_score_pct =
      mean(completeness_adjusted_score_pct),
    mean_score_delta_pct =
      mean(completeness_adjusted_score_pct - raw_score_pct),
    raw_prevalence_ge_50_pct = 100 * mean(raw_ge_50),
    adjusted_prevalence_ge_50_pct = 100 * mean(adjusted_ge_50),
    prevalence_ge_50_delta_pct =
      adjusted_prevalence_ge_50_pct - raw_prevalence_ge_50_pct,
    raw_prevalence_ge_75_pct = 100 * mean(raw_ge_75),
    adjusted_prevalence_ge_75_pct = 100 * mean(adjusted_ge_75),
    prevalence_ge_75_delta_pct =
      adjusted_prevalence_ge_75_pct - raw_prevalence_ge_75_pct,
    .groups = "drop"
  )

model_data <- adjusted_scores %>%
  inner_join(
    mapping %>%
      filter(MAG != "S1_9-12_bin1") %>%
      distinct(MAG, habitat, final_clade),
    by = "MAG"
  ) %>%
  mutate(
    habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS")),
    final_clade = factor(final_clade, levels = paste("Clade", 1:4))
  )

fit_feature_model <- function(df) {
  if (n_distinct(df$reconstruction_score_pct) < 2L) {
    return(tibble(
      model_status = "invariant outcome",
      n_MAG = n_distinct(df$MAG),
      design_rank = NA_integer_,
      design_columns = NA_integer_,
      rank_deficient = NA,
      completeness_slope_per_10pct = NA_real_,
      completeness_p = NA_real_
    ))
  }
  fit <- tryCatch(
    lm(
      reconstruction_score_pct ~ completeness + habitat + final_clade,
      data = df
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(
      model_status = "model failed",
      n_MAG = n_distinct(df$MAG),
      design_rank = NA_integer_,
      design_columns = NA_integer_,
      rank_deficient = NA,
      completeness_slope_per_10pct = NA_real_,
      completeness_p = NA_real_
    ))
  }
  design <- model.matrix(fit)
  coefficients <- summary(fit)$coefficients
  tibble(
    model_status = "descriptive OLS; dual-copy MAG excluded",
    n_MAG = n_distinct(df$MAG),
    design_rank = fit$rank,
    design_columns = ncol(design),
    rank_deficient = fit$rank < ncol(design),
    completeness_slope_per_10pct =
      unname(coefficients["completeness", "Estimate"]) * 10,
    completeness_p = unname(coefficients["completeness", "Pr(>|t|)"])
  )
}

model_diagnostics <- model_data %>%
  group_by(research_category, feature_name, source_type, score_semantics) %>%
  group_modify(~ fit_feature_model(.x)) %>%
  ungroup() %>%
  mutate(
    completeness_BH = p.adjust(completeness_p, method = "BH"),
    interpretation_limit = paste(
      "Descriptive sensitivity only; sparse and structurally empty",
      "habitat-clade cells limit adjusted inference"
    )
  )

locked_features <- c(
  "Arsenate detoxification/resistance (arsC) marker carriage",
  "Respiratory arsenate reduction (arrA) marker carriage",
  "Arsenite oxidation (arxA/aioA) marker carriage",
  "Cobalamin biosynthesis, cobinamide => cobalamin"
)
locked_clade4_is_sensitivity <- clade_habitat_sensitivity %>%
  filter(
    habitat == "IS",
    final_clade == "Clade 4",
    feature_name %in% locked_features
  )
if (nrow(locked_clade4_is_sensitivity) != 4L ||
    any(locked_clade4_is_sensitivity$n_unique_MAG != 20L)) {
  stop("Locked Clade 4-IS completeness sensitivity is incomplete.")
}

write.table(
  adjusted_scores,
  file.path(out_dir, "completeness_adjusted_MAG_feature_scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  feature_sensitivity,
  file.path(out_dir, "completeness_adjustment_sensitivity_by_feature.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  clade_habitat_sensitivity,
  file.path(out_dir, "completeness_adjusted_clade_habitat_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  model_diagnostics,
  file.path(out_dir, "completeness_model_diagnostics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  locked_clade4_is_sensitivity,
  file.path(out_dir, "Clade4_IS_completeness_sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Completeness-adjusted MAG-feature rows:", nrow(adjusted_scores), "\n")
cat("Threshold flips at 50%:", sum(adjusted_scores$threshold_flip_50), "\n")
cat("Threshold flips at 75%:", sum(adjusted_scores$threshold_flip_75), "\n")
cat(
  "Feature models rank deficient:",
  sum(model_diagnostics$rank_deficient, na.rm = TRUE), "/",
  nrow(model_diagnostics), "\n"
)
cat(
  "Clade 4-IS locked rows retained:",
  nrow(locked_clade4_is_sensitivity), "\n"
)

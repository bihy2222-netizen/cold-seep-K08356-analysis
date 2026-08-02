suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
project_dir <- dirname(script_dir)
out_dir <- file.path(script_dir, "results")

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) paths[1] else hit[1]
}
classification_path <- Sys.getenv(
  "K08356_CLADE_TSV",
  first_existing(c(
    file.path(script_dir, "input", "four_clade_classification_49.tsv"),
    file.path(project_dir, "four_clade_classification_49.tsv")
  ))
)
neighborhood_path <- Sys.getenv(
  "K08356_NEIGHBORHOOD_TSV",
  first_existing(c(
    file.path(script_dir, "input", "gene_neighborhood_49_summary.tsv"),
    file.path(project_dir, "gene_neighborhood_49_summary.tsv")
  ))
)
tree_svg_path <- Sys.getenv(
  "K08356_TREE_SVG",
  first_existing(c(
    file.path(script_dir, "input", "small_tree_MAG49_four_clade_with_gene_islands.svg"),
    file.path(project_dir, "figures", "small_tree_MAG49_four_clade_with_gene_islands.svg")
  ))
)
mapping_path <- file.path(
  out_dir, "K08356seq49_MAG48_clade_habitat_mapping.tsv"
)
quality_path <- file.path(out_dir, "K08356seq49_MAG48_quality_audit.tsv")
score_path <- file.path(out_dir, "K08356seq49_MAG48_metabolic_feature_scores.tsv")
group_summary_path <- file.path(out_dir, "habitat_clade_metabolic_summary.tsv")
canonical_hmm_path <- file.path(
  out_dir, "canonical_AioA_METABOLIC_HMM_score_audit.tsv"
)
denominator_sensitivity_path <- file.path(
  out_dir, "gene_set_unavailable_KO_threshold_sensitivity_by_feature.tsv"
)
dual_copy_summary_path <- file.path(
  out_dir, "dual_copy_MAG_exclusion_sensitivity_summary.tsv"
)
completeness_sensitivity_path <- file.path(
  out_dir, "completeness_adjustment_sensitivity_by_feature.tsv"
)
locked_completeness_path <- file.path(
  out_dir, "Clade4_IS_completeness_sensitivity.tsv"
)
gtdb_taxonomy_path <- file.path(
  out_dir, "K08356_MAG48_GTDB_r226_taxonomy.tsv"
)
family_confounding_path <- file.path(
  out_dir, "GTDB_family_confounding_audit.tsv"
)
family_estimability_path <- file.path(
  out_dir, "GTDB_family_model_estimability.tsv"
)
rhodo_comparison_path <- file.path(
  out_dir, "Rhodobacteraceae_IS_Clade3_vs_Clade4_locked_features.tsv"
)

required <- c(
  classification_path, neighborhood_path, tree_svg_path, mapping_path,
  quality_path, score_path, group_summary_path, canonical_hmm_path,
  denominator_sensitivity_path, dual_copy_summary_path,
  completeness_sensitivity_path, locked_completeness_path,
  gtdb_taxonomy_path, family_confounding_path, family_estimability_path,
  rhodo_comparison_path
)
if (any(!file.exists(required))) {
  stop("Missing closeout-audit input: ", paste(required[!file.exists(required)], collapse = ", "))
}

classification <- read.delim(classification_path, check.names = FALSE)
canonical_ids <- sort(unique(classification$Sequence_ID))
if (length(canonical_ids) != 49L ||
    n_distinct(classification$MAG_ID) != 48L ||
    !identical(
      as.integer(table(factor(classification$clade_display, levels = c(
        "Canonical AioA", "Uncertain DMSOR", "IdrA phylogenetic",
        "DIRM-synteny IdrA"
      )))),
      c(3L, 7L, 19L, 20L)
    )) {
  stop("Canonical classification is not the expected 49/48, 3/7/19/20 set.")
}

set_audit <- function(artifact, path, observed_ids, note) {
  observed_ids <- sort(unique(observed_ids))
  tibble(
    artifact,
    path,
    sequence_IDs_expected = length(canonical_ids),
    sequence_IDs_found = length(intersect(canonical_ids, observed_ids)),
    missing_IDs = length(setdiff(canonical_ids, observed_ids)),
    extra_IDs = length(setdiff(observed_ids, canonical_ids)),
    status = if_else(setequal(canonical_ids, observed_ids), "verified", "mismatch"),
    note
  )
}

neighborhood <- read.delim(neighborhood_path, check.names = FALSE)
mapping <- read.delim(mapping_path, check.names = FALSE)
tree_svg <- paste(readLines(tree_svg_path, warn = FALSE), collapse = "\n")
tree_ids <- canonical_ids[vapply(
  canonical_ids, function(id) str_detect(tree_svg, fixed(id)), logical(1)
)]

artifact_audit <- bind_rows(
  set_audit(
    "Canonical classification", "../four_clade_classification_49.tsv",
    classification$Sequence_ID,
    "Authoritative 49-sequence table; clade counts 3/7/19/20"
  ),
  set_audit(
    "Four-clade tree and gene-island SVG",
    "../figures/small_tree_MAG49_four_clade_with_gene_islands.svg",
    tree_ids, "Every canonical Sequence_ID occurs in the SVG text"
  ),
  set_audit(
    "Gene-neighborhood summary", "../gene_neighborhood_49_summary.tsv",
    neighborhood$target_id,
    "Exact target_id set matches the canonical classification"
  ),
  set_audit(
    "Metabolism clade-habitat mapping",
    "K08356seq49_MAG48_clade_habitat_mapping.tsv", mapping$contig_gene,
    "Exact contig_gene set matches the canonical classification"
  ),
  tibble(
    artifact = "Abundance figure/table",
    path = "not found in joint_reference_search_49_update_20260801",
    sequence_IDs_expected = 49L,
    sequence_IDs_found = NA_integer_,
    missing_IDs = NA_integer_,
    extra_IDs = NA_integer_,
    status = "not verified",
    note = paste(
      "No abundance/TPM artifact is present in this current 49-sequence folder;",
      "do not claim synchronization until regenerated or supplied"
    )
  )
)
if (any(artifact_audit$status[1:4] != "verified")) {
  stop("One or more current 49-sequence artifacts failed ID synchronization.")
}

scores <- read.delim(score_path, check.names = FALSE)
if (nrow(scores) != 1920L ||
    n_distinct(scores$MAG) != 48L ||
    n_distinct(scores$feature_name) != 40L ||
    anyDuplicated(scores[c("MAG", "feature_name")])) {
  stop("Metabolic score matrix failed 48 x 40 uniqueness assertions.")
}

canonical_hmm <- read.delim(canonical_hmm_path, check.names = FALSE)
if (nrow(canonical_hmm) != 3L ||
    sum(canonical_hmm$METABOLIC_aioA_reported_in_original_tblout) != 2L ||
    sum(canonical_hmm$joint_screen_Combined_640_pass) != 3L) {
  stop("Canonical AioA evidence must remain 2/3 at METABOLIC -T800 and 3/3 at the independent joint-screen >=640.")
}

denominator_sensitivity <- read.delim(
  denominator_sensitivity_path, check.names = FALSE
)
if (sum(denominator_sensitivity$n_MAG_threshold_flip_50) != 41L ||
    sum(denominator_sensitivity$n_MAG_threshold_flip_75) != 10L) {
  stop("Unavailable-KO sensitivity no longer matches the audited 41/10 threshold flips.")
}

dual_copy_summary <- read.delim(dual_copy_summary_path, check.names = FALSE)
if (nrow(dual_copy_summary) != 1L ||
    dual_copy_summary$cells_becoming_empty != 80L ||
    dual_copy_summary$populated_cells_with_numeric_change != 0L ||
    dual_copy_summary$Clade4_IS_affected) {
  stop("Dual-copy MAG exclusion no longer matches the locked sensitivity result.")
}

group_summary <- read.delim(group_summary_path, check.names = FALSE)
clade4_is_locked_results <- group_summary %>%
  filter(
    habitat == "IS",
    final_clade == "Clade 4",
    feature_name %in% c(
      "Arsenate detoxification/resistance (arsC) marker carriage",
      "Respiratory arsenate reduction (arrA) marker carriage",
      "Arsenite oxidation (arxA/aioA) marker carriage",
      "Cobalamin biosynthesis, cobinamide => cobalamin"
    )
  ) %>%
  transmute(
    habitat,
    final_clade,
    feature_name,
    n_unique_MAG,
    n_MAG_ge_50 = round(n_unique_MAG * prevalence_score_ge_50_pct / 100),
    prevalence_ge_50_pct = prevalence_score_ge_50_pct,
    prevalence_ge_75_pct = prevalence_score_ge_75_pct,
    mean_score_or_binary_carriage_pct = mean_reconstruction_score_pct,
    interpretation = case_when(
      str_detect(feature_name, fixed("arsC")) ~
        "Widespread arsenic detoxification/resistance marker carriage",
      str_detect(feature_name, fixed("arrA")) ~
        "Respiratory arsenate reduction marker was uncommon",
      str_detect(feature_name, fixed("arxA/aioA")) ~
        "No independent arsenite-oxidation marker remained after focal exclusion",
      TRUE ~
        "All MAGs exceeded the legacy 50% partial-reconstruction threshold; not evidence of a complete pathway"
    )
  )
if (nrow(clade4_is_locked_results) != 4L ||
    !setequal(clade4_is_locked_results$n_MAG_ge_50, c(17L, 2L, 0L, 20L)) ||
    clade4_is_locked_results$mean_score_or_binary_carriage_pct[
      clade4_is_locked_results$feature_name ==
        "Cobalamin biosynthesis, cobinamide => cobalamin"
    ] != 79.375) {
  stop("Clade 4-IS locked arsenic/B12 result changed unexpectedly.")
}

quality <- read.delim(quality_path, check.names = FALSE)
if (nrow(quality) != 48L ||
    any(quality$completeness < 50) ||
    any(quality$contamination > 10)) {
  stop("MAG quality audit failed the expected 48-MAG 50/10 threshold.")
}

completeness_sensitivity <- read.delim(
  completeness_sensitivity_path, check.names = FALSE
)
locked_completeness <- read.delim(
  locked_completeness_path, check.names = FALSE
)
if (nrow(completeness_sensitivity) != 40L ||
    sum(completeness_sensitivity$n_MAG_threshold_flip_50) != 83L ||
    sum(completeness_sensitivity$n_MAG_threshold_flip_75) != 35L ||
    nrow(locked_completeness) != 4L ||
    any(
      locked_completeness$raw_prevalence_ge_50_pct !=
        locked_completeness$adjusted_prevalence_ge_50_pct
    )) {
  stop("Completeness sensitivity no longer matches the audited result.")
}

gtdb_taxonomy <- read.delim(gtdb_taxonomy_path, check.names = FALSE)
family_confounding <- read.delim(
  family_confounding_path, check.names = FALSE
)
family_estimability <- read.delim(
  family_estimability_path, check.names = FALSE
)
rhodo_comparison <- read.delim(rhodo_comparison_path, check.names = FALSE)
if (nrow(gtdb_taxonomy) != 48L ||
    anyDuplicated(gtdb_taxonomy$MAG) ||
    !setequal(gtdb_taxonomy$MAG, unique(mapping$MAG)) ||
    nrow(family_estimability) != 4L ||
    nrow(rhodo_comparison) != 4L) {
  stop("GTDB family audit does not match the exact 48-MAG analysis set.")
}

n_multihabitat_families <- sum(family_confounding$n_habitats > 1L)
family_habitat_df <- family_estimability$independently_estimable_df[
  family_estimability$term == "habitat"
]
family_clade_df <- family_estimability$independently_estimable_df[
  family_estimability$term == "final_clade"
]

confounding_status <- tibble(
  item = c(
    "MAG quality threshold", "Completeness balance",
    "Completeness sensitivity", "GTDB family-level host table",
    "Family/habitat/completeness-adjusted model"
  ),
  status = c(
    "complete", "imbalanced", "complete as sensitivity", "complete",
    "partially estimable; strongly confounded"
  ),
  evidence = c(
    "48/48 MAGs pass completeness >=50% and contamination <=10%",
    paste0(
      "Completeness range ", min(quality$completeness), "%-",
      max(quality$completeness),
      "%; habitat medians AS 60.79%, ES 71.22%, IS 88.25%, NS 79.02%"
    ),
    paste0(
      "Random-gene-loss scaling yields ",
      sum(completeness_sensitivity$n_MAG_threshold_flip_50),
      " flips at 50% and ",
      sum(completeness_sensitivity$n_MAG_threshold_flip_75),
      " at 75%; locked Clade 4-IS 50% prevalences are unchanged"
    ),
    paste0(
      "GTDB-Tk r226 classification matches 48/48 MAGs; ",
      n_distinct(gtdb_taxonomy$family, na.rm = TRUE),
      " assigned families"
    ),
    paste0(
      n_multihabitat_families,
      " families span more than one habitat; independent model df: habitat ",
      family_habitat_df, ", clade ", family_clade_df
    )
  ),
  allowed_interpretation = c(
    "All analyzed MAGs meet the stated inclusion threshold",
    "Passing 50/10 does not eliminate completeness confounding",
    "Use adjusted scores only as a sensitivity analysis under random gene loss",
    "Use GTDB family for host taxonomy; do not substitute CheckM lineage",
    "Report estimable within-family sensitivities, but do not claim that host or habitat confounding was eliminated"
  )
)

final_lock_status <- tibble(
  analysis_component = c(
    "Metabolic-potential figure",
    "Abundance/TPM figure",
    "Completeness sensitivity",
    "GTDB-family sensitivity"
  ),
  status = c(
    "verified", "not verified", "verified as sensitivity",
    "verified; residual confounding remains"
  ),
  permitted_statement = c(
    "Unadjusted genomic metabolic-potential association",
    "No abundance conclusion until the 56-sample TPM artifact is recovered and audited",
    "Locked Clade 4-IS result is stable at 50%; broad pathway patterns remain completeness-sensitive",
    "Exact family taxonomy is available, but sparse and nested family-habitat structure prevents a claim of fully adjusted habitat/clade effects"
  )
)

write.table(
  artifact_audit, file.path(out_dir, "artifact_49_sequence_sync_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  confounding_status, file.path(out_dir, "confounding_adjustment_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  final_lock_status, file.path(out_dir, "FINAL_LOCK_STATUS.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  clade4_is_locked_results,
  file.path(out_dir, "Clade4_IS_locked_arsenic_B12_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Verified current artifacts at 49/49 IDs: 4\n")
cat("Abundance artifact: not found / not verified\n")
cat("Metabolic matrix: 48 MAG x 40 features = 1920 unique rows\n")
cat("Unavailable-KO threshold flips: 41 at 50%; 10 at 75%\n")
cat("Dual-copy MAG exclusion: two n=1 AS cells become empty; Clade 4-IS unchanged\n")
cat("Completeness sensitivity: 83 flips at 50%; 35 at 75%; locked Clade 4-IS 50% result unchanged\n")
cat("GTDB taxonomy: exact 48/48 MAG match; family-adjusted model rank audited\n")
cat("Quality/host sensitivity is complete, but residual confounding remains\n")

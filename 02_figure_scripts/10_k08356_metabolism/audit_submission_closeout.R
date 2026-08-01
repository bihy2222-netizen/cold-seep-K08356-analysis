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

required <- c(
  classification_path, neighborhood_path, tree_svg_path, mapping_path,
  quality_path, score_path
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

quality <- read.delim(quality_path, check.names = FALSE)
if (nrow(quality) != 48L ||
    any(quality$completeness < 50) ||
    any(quality$contamination > 10)) {
  stop("MAG quality audit failed the expected 48-MAG 50/10 threshold.")
}

confounding_status <- tibble(
  item = c(
    "MAG quality threshold", "Completeness balance",
    "Completeness-adjusted model", "GTDB family-level host table",
    "Family/habitat/completeness-adjusted model"
  ),
  status = c(
    "complete", "imbalanced", "not performed", "unavailable",
    "not estimable in current package"
  ),
  evidence = c(
    "48/48 MAGs pass completeness >=50% and contamination <=10%",
    paste0(
      "Completeness range ", min(quality$completeness), "%-",
      max(quality$completeness),
      "%; habitat medians AS 60.79%, ES 71.22%, IS 88.25%, NS 79.02%"
    ),
    "No regression or completeness-normalized gene-set model is included",
    "No validated family table matched to the exact 48 MAG IDs was found in inspected inputs",
    "Validated family covariate is absent and habitat/clade cells are strongly unbalanced"
  ),
  allowed_interpretation = c(
    "All analyzed MAGs meet the stated inclusion threshold",
    "Passing 50/10 does not eliminate completeness confounding",
    "Do not state that completeness effects were excluded",
    "Do not substitute CheckM lineage for GTDB family",
    "Describe patterns as unadjusted genomic-potential associations only"
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

cat("Verified current artifacts at 49/49 IDs: 4\n")
cat("Abundance artifact: not found / not verified\n")
cat("Metabolic matrix: 48 MAG x 40 features = 1920 unique rows\n")
cat("Quality: 48/48 MAGs pass 50/10; confounding adjustment not performed\n")

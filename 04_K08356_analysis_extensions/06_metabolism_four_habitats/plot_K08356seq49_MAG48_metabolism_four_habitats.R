suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(script_arg)) {
  script_path <- sub("^--file=", "", script_arg[1])
  script_path <- gsub("~\\+~", " ", script_path)
  dirname(normalizePath(script_path))
} else {
  getwd()
}
first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) paths[1] else hit[1]
}
font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)

input_xlsx <- Sys.getenv(
  "METABOLIC_RESULT_XLSX",
  file.path(script_dir, "METABOLIC_group_derep49_result.xlsx")
)
classification_tsv <- Sys.getenv(
  "K08356_CLADE_TSV",
  first_existing(c(
    file.path(script_dir, "input", "four_clade_classification_49.tsv"),
    file.path(dirname(script_dir), "four_clade_classification_49.tsv")
  ))
)
manifest_tsv <- Sys.getenv(
  "MAG49_MANIFEST_TSV",
  first_existing(c(
    file.path(script_dir, "MAG49_input_manifest.tsv"),
    file.path(script_dir, "input", "MAG49_input_manifest.tsv"),
    file.path(dirname(dirname(script_dir)), "01_metagenomic_workflow",
              "04_metabolism", "METABOLIC_group_derep49", "manifests",
              "MAG49_input_manifest.tsv")
  ))
)
sample_group_csv <- Sys.getenv(
  "SAMPLE_GROUP_CSV",
  file.path(script_dir, "sample_group_corrected.csv")
)
kegg_identifier_dir <- Sys.getenv(
  "KEGG_IDENTIFIER_DIR",
  file.path(script_dir, "KEGG_identifier_result")
)
gene_set_definition_tsv <- Sys.getenv(
  "GENE_SET_DEFINITION_TSV",
  file.path(script_dir, "gene_set_definitions_previous_algorithm.tsv")
)
canonical_aioa_tblout <- Sys.getenv(
  "CANONICAL_AIOA_TBLOUT",
  file.path(script_dir, "quality", "raw_hmm",
            "canonical3_aioA_no_threshold.tblout")
)
metabolic_aioa_tblout <- Sys.getenv(
  "METABOLIC_AIOA_TBLOUT",
  file.path(script_dir, "quality", "raw_hmm",
            "aioA.hmm.total.hmmsearch_result.txt")
)
out_dir <- Sys.getenv(
  "METABOLIC_PLOT_OUT_DIR",
  file.path(script_dir, "results")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  input_xlsx,
  classification_tsv,
  manifest_tsv,
  sample_group_csv,
  gene_set_definition_tsv,
  canonical_aioa_tblout,
  metabolic_aioa_tblout
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "))
}
if (!dir.exists(kegg_identifier_dir)) {
  stop("Missing KEGG identifier result directory: ", kegg_identifier_dir)
}

clade_levels <- paste("Clade", 1:4)
clade_from_display <- c(
  "Canonical AioA" = "Clade 1",
  "Uncertain DMSOR" = "Clade 2",
  "IdrA phylogenetic" = "Clade 3",
  "DIRM-synteny IdrA" = "Clade 4"
)
clade_names <- c(
  "Clade 1" = "Canonical AioA-associated",
  "Clade 2" = "Unknown/uncertain DMSOR",
  "Clade 3" = "IdrA-associated without complete DIRM-like synteny",
  "Clade 4" = "IdrA-associated with DIRM-like synteny"
)

classification <- read.delim(
  classification_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
manifest <- read.delim(
  manifest_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sample_groups <- read.csv(
  sample_group_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  transmute(sample_id = sample.id, habitat_in_sample_metadata = group)

mapping <- classification %>%
  transmute(
    contig_gene = Sequence_ID,
    MAG = MAG_ID,
    habitat_in_classification = Habitat,
    clade_display,
    final_clade = unname(clade_from_display[clade_display]),
    final_clade_name = unname(clade_names[final_clade])
  ) %>%
  left_join(
    manifest %>%
      transmute(
        MAG = bin_id,
        sample_id = str_remove(bin_id, "_bin[^_]+$"),
        habitat_in_manifest = group
      ),
    by = "MAG"
  ) %>%
  left_join(sample_groups, by = "sample_id") %>%
  mutate(
    habitat = coalesce(habitat_in_sample_metadata, habitat_in_manifest),
    habitat_assignment_source = if_else(
      !is.na(habitat_in_sample_metadata),
      "corrected sample metadata",
      "group-derep manifest fallback"
    ),
    final_clade = factor(final_clade, levels = clade_levels),
    habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS")),
    classification_habitat_match = as.character(habitat) ==
      habitat_in_classification,
    manifest_habitat_match = is.na(habitat_in_sample_metadata) |
      habitat_in_sample_metadata == habitat_in_manifest
  )

expected_clade_counts <- c(
  "Clade 1" = 3L,
  "Clade 2" = 7L,
  "Clade 3" = 19L,
  "Clade 4" = 20L
)
observed_clade_counts <- table(mapping$final_clade)
if (nrow(mapping) != 49L || n_distinct(mapping$MAG) != 48L) {
  stop("Expected 49 K08356 sequence memberships from 48 unique MAGs.")
}
if (any(is.na(mapping$final_clade)) || any(is.na(mapping$habitat))) {
  stop("Missing clade or manifest-based habitat assignment.")
}
if (!all(observed_clade_counts[names(expected_clade_counts)] ==
         expected_clade_counts)) {
  stop("Expected clade sequence counts 3/7/19/20.")
}
if (!setequal(unique(mapping$MAG), manifest$bin_id)) {
  stop("The clade classification and 48-MAG input manifest do not match.")
}
if (any(!mapping$manifest_habitat_match)) {
  stop("Sample metadata and group-derep manifest disagree for one or more MAGs.")
}

normalize_metabolic_mag <- function(x, suffix) {
  x %>%
    str_remove(fixed(suffix)) %>%
    str_replace_all(fixed("."), "-")
}

gene_set_definitions <- read.delim(
  gene_set_definition_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    feature_name = if_else(
      feature_name == "Flagellar Assembly",
      "Flagellar Assembly (partial 36-KO reconstruction)",
      feature_name
    )
  ) %>%
  distinct(source_sheet, research_category, feature_name, gene_id)

if (nrow(gene_set_definitions) != 510L ||
    n_distinct(gene_set_definitions$feature_name) != 38L) {
  stop("Expected the previous algorithm's 510 rows / 38 gene sets.")
}

kegg_result_files <- list.files(
  kegg_identifier_dir,
  pattern = "\\.cds\\.result\\.txt$",
  full.names = TRUE
)
if (length(kegg_result_files) != 48L) {
  stop("Expected 48 METABOLIC KEGG result files; observed ",
       length(kegg_result_files), ".")
}

ko_hits <- bind_rows(lapply(kegg_result_files, function(path) {
  read.delim(
    path,
    header = FALSE,
    sep = "\t",
    fill = TRUE,
    quote = "",
    col.names = c("gene_id", "gene_count"),
    stringsAsFactors = FALSE
  ) %>%
    transmute(
      MAG = basename(path) %>% str_remove("\\.cds\\.result\\.txt$"),
      gene_id,
      gene_count = replace_na(suppressWarnings(as.numeric(gene_count)), 0)
    )
}))

definition_id_detail <- gene_set_definitions %>%
  mutate(
    identifier_type = if_else(
      str_detect(gene_id, "^K[0-9]{5}$"), "KO", "non-KO"
    ),
    represented_in_KO_result_catalog = gene_id %in% unique(ko_hits$gene_id),
    denominator_treatment = case_when(
      research_category == "Arsenic" ~
        "Legacy arsenic set excluded; replaced by independent HMM-marker rows",
      identifier_type == "KO" ~ "Included in predefined gene-set denominator",
      TRUE ~ "Invalid outside arsenic set"
    ),
    annotation_source = case_when(
      identifier_type == "KO" & represented_in_KO_result_catalog ~
        "METABOLIC KEGG identifier result catalog",
      identifier_type == "KO" & !represented_in_KO_result_catalog ~
        "Legacy predefined gene-set table; absent from current KO result catalog",
      gene_id == "EC:1.20.4.1" ~
        "IUBMB/ENZYME: ArsC-type arsenate reductase (glutathione/glutaredoxin); detoxification",
      gene_id == "EC:1.20.99.1" ~
        "IUBMB/ENZYME: ArrA-type arsenate reductase (donor); respiratory marker",
      TRUE ~ "Unresolved non-KO identifier"
    ),
    annotation_reference = case_when(
      gene_id == "EC:1.20.4.1" ~
        "https://enzyme.expasy.org/EC/1.20.4.1",
      gene_id == "EC:1.20.99.1" ~
        "https://enzyme.expasy.org/EC/1.20.99.1",
      TRUE ~ NA_character_
    )
  )

definition_id_summary <- bind_rows(
  tibble(
    metric = c(
      "definition_records_total",
      "definition_records_KO",
      "definition_records_non_KO",
      "unique_KO_identifiers",
      "unique_non_KO_identifiers",
      "KO_records_absent_from_KO_result_catalog",
      "non_KO_records_absent_from_KO_result_catalog"
    ),
    value = c(
      nrow(definition_id_detail),
      sum(definition_id_detail$identifier_type == "KO"),
      sum(definition_id_detail$identifier_type == "non-KO"),
      n_distinct(definition_id_detail$gene_id[
        definition_id_detail$identifier_type == "KO"
      ]),
      n_distinct(definition_id_detail$gene_id[
        definition_id_detail$identifier_type == "non-KO"
      ]),
      sum(definition_id_detail$identifier_type == "KO" &
            !definition_id_detail$represented_in_KO_result_catalog),
      sum(definition_id_detail$identifier_type == "non-KO" &
            !definition_id_detail$represented_in_KO_result_catalog)
    )
  )
)
if (sum(definition_id_detail$identifier_type == "KO") != 508L ||
    sum(definition_id_detail$identifier_type == "non-KO") != 2L) {
  stop("Expected 508 KO records and two legacy non-KO arsenic records.")
}

definition_denominator_impact <- definition_id_detail %>%
  group_by(source_sheet, research_category, feature_name) %>%
  summarise(
    definition_records = n(),
    records_represented_in_KO_result_catalog =
      sum(represented_in_KO_result_catalog),
    unavailable_definition_records =
      sum(!represented_in_KO_result_catalog),
    distinct_unavailable_identifiers = n_distinct(
      gene_id[!represented_in_KO_result_catalog]
    ),
    maximum_attainable_coverage_pct_in_current_KO_catalog =
      100 * records_represented_in_KO_result_catalog / definition_records,
    unavailable_identifiers = paste(
      sort(unique(gene_id[!represented_in_KO_result_catalog])),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  filter(unavailable_definition_records > 0)

selected_modules <- gene_set_definitions %>%
  filter(research_category != "Arsenic") %>%
  distinct(source_sheet, research_category, feature_name)

non_ko_gene_ids <- gene_set_definitions %>%
  filter(research_category != "Arsenic", !str_detect(gene_id, "^K[0-9]{5}$")) %>%
  pull(gene_id) %>%
  unique()
if (length(non_ko_gene_ids)) {
  stop("Non-KO IDs found outside the separately scored arsenic set: ",
       paste(non_ko_gene_ids, collapse = ", "))
}

gene_set_scores_detailed <- tidyr::crossing(
  MAG = unique(mapping$MAG),
  gene_set_definitions %>% filter(research_category != "Arsenic")
) %>%
  left_join(ko_hits, by = c("MAG", "gene_id")) %>%
  mutate(gene_count = replace_na(gene_count, 0)) %>%
  group_by(MAG, source_sheet, research_category, feature_name) %>%
  summarise(
    feature_component_count = n_distinct(gene_id),
    components_detected = n_distinct(gene_id[gene_count > 0]),
    reconstruction_score_pct = 100 * components_detected /
      feature_component_count,
    .groups = "drop"
  )

gene_set_scores <- gene_set_scores_detailed %>%
  transmute(
    MAG,
    research_category,
    feature_name,
    source_type = "Previous algorithm: distinct-gene set coverage",
    feature_component_count,
    components_detected,
    reconstruction_score_pct
  )

catalog_compatible_gene_set_scores <- tidyr::crossing(
  MAG = unique(mapping$MAG),
  definition_id_detail %>%
    filter(
      research_category != "Arsenic",
      represented_in_KO_result_catalog
    ) %>%
    select(source_sheet, research_category, feature_name, gene_id)
) %>%
  left_join(ko_hits, by = c("MAG", "gene_id")) %>%
  mutate(gene_count = replace_na(gene_count, 0)) %>%
  group_by(MAG, source_sheet, research_category, feature_name) %>%
  summarise(
    catalog_compatible_component_count = n_distinct(gene_id),
    catalog_compatible_components_detected =
      n_distinct(gene_id[gene_count > 0]),
    catalog_compatible_score_pct =
      100 * catalog_compatible_components_detected /
      catalog_compatible_component_count,
    .groups = "drop"
  )

affected_definition_features <- definition_denominator_impact %>%
  filter(research_category != "Arsenic") %>%
  select(
    source_sheet, research_category, feature_name,
    definition_records, unavailable_definition_records,
    distinct_unavailable_identifiers, unavailable_identifiers
  )

gene_set_denominator_sensitivity_MAG <- gene_set_scores_detailed %>%
  inner_join(
    affected_definition_features,
    by = c("source_sheet", "research_category", "feature_name")
  ) %>%
  rename(
    legacy_component_count = feature_component_count,
    legacy_components_detected = components_detected,
    legacy_score_pct = reconstruction_score_pct
  ) %>%
  left_join(
    catalog_compatible_gene_set_scores,
    by = c("MAG", "source_sheet", "research_category", "feature_name")
  ) %>%
  mutate(
    unavailable_definition_pct =
      100 * unavailable_definition_records / definition_records,
    score_delta_pct = catalog_compatible_score_pct - legacy_score_pct,
    score_changed = abs(score_delta_pct) > 1e-10,
    legacy_ge_50 = legacy_score_pct >= 50,
    catalog_compatible_ge_50 = catalog_compatible_score_pct >= 50,
    threshold_flip_50 = legacy_ge_50 != catalog_compatible_ge_50,
    legacy_ge_75 = legacy_score_pct >= 75,
    catalog_compatible_ge_75 = catalog_compatible_score_pct >= 75,
    threshold_flip_75 = legacy_ge_75 != catalog_compatible_ge_75
  )

gene_set_denominator_sensitivity_feature <-
  gene_set_denominator_sensitivity_MAG %>%
  group_by(source_sheet, research_category, feature_name) %>%
  summarise(
    definition_KO_total = first(definition_records),
    unavailable_KO_count = first(unavailable_definition_records),
    unavailable_KO_pct = first(unavailable_definition_pct),
    unavailable_KO_IDs = first(unavailable_identifiers),
    n_MAG_evaluated = n_distinct(MAG),
    n_MAG_score_changed = sum(score_changed),
    n_MAG_legacy_ge_50 = sum(legacy_ge_50),
    n_MAG_catalog_compatible_ge_50 = sum(catalog_compatible_ge_50),
    n_MAG_threshold_flip_50 = sum(threshold_flip_50),
    n_MAG_legacy_ge_75 = sum(legacy_ge_75),
    n_MAG_catalog_compatible_ge_75 = sum(catalog_compatible_ge_75),
    n_MAG_threshold_flip_75 = sum(threshold_flip_75),
    maximum_score_increase_pct = max(score_delta_pct),
    .groups = "drop"
  )

gene_set_denominator_threshold_flips <-
  gene_set_denominator_sensitivity_MAG %>%
  filter(threshold_flip_50 | threshold_flip_75)

hmm_hit <- read_excel(
  input_xlsx,
  sheet = "HMMHitNum",
  .name_repair = "minimal"
)
hmm_hit_suffix <- ".cds.Hits"
hmm_hit_columns <- names(hmm_hit)[endsWith(names(hmm_hit), hmm_hit_suffix)]
if (length(hmm_hit_columns) != 48L) {
  stop("Expected 48 MAG hit columns in HMMHitNum; observed ",
       length(hmm_hit_columns), ".")
}

arsenic_features <- c(
  "Respiratory arsenate reduction (arrA) marker carriage",
  "Arsenate detoxification/resistance (arsC) marker carriage",
  "Arsenite oxidation (arxA/aioA) marker carriage"
)
arsenic_hit_qc <- hmm_hit %>%
  filter(
    Category == "As cycling",
    Gene.abbreviation %in% c(
      "arrA", "arsC (grx)", "arsC (trx)", "arxA", "aioA"
    )
  ) %>%
  select(Function, Gene.abbreviation, Gene.name, Hmm.file,
         all_of(hmm_hit_columns)) %>%
  pivot_longer(
    all_of(hmm_hit_columns),
    names_to = "metabolic_column",
    values_to = "gene_hits"
  ) %>%
  mutate(
    MAG = normalize_metabolic_mag(metabolic_column, hmm_hit_suffix),
    gene_hits = as.character(gene_hits),
    feature_name = case_when(
      Gene.abbreviation == "arrA" ~
        "Respiratory arsenate reduction (arrA) marker carriage",
      Gene.abbreviation %in% c("arsC (grx)", "arsC (trx)") ~
        "Arsenate detoxification/resistance (arsC) marker carriage",
      Gene.abbreviation %in% c("arxA", "aioA") ~
        "Arsenite oxidation (arxA/aioA) marker carriage",
      TRUE ~ NA_character_
    ),
    marker_component = case_when(
      Gene.abbreviation == "arrA" ~ "dissimilatory arsenate reductase",
      Gene.abbreviation %in% c("arsC (grx)", "arsC (trx)") ~
        "detoxification arsenate reductase",
      Gene.abbreviation %in% c("arxA", "aioA") ~
        "arsenite oxidase large subunit",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(gene_hits), gene_hits != "", gene_hits != "None") %>%
  separate_longer_delim(gene_hits, delim = ",") %>%
  transmute(
    feature_name,
    MAG,
    marker = Gene.abbreviation,
    marker_component,
    hmm_file = Hmm.file,
    gene_id = str_trim(gene_hits),
    is_focal_K08356_sequence = gene_id %in% mapping$contig_gene,
    retained_after_focal_exclusion = !is_focal_K08356_sequence
  )

arsenic_nonfocal_counts <- arsenic_hit_qc %>%
  filter(retained_after_focal_exclusion) %>%
  group_by(feature_name, MAG) %>%
  summarise(
    components_detected = n_distinct(marker_component),
    nonfocal_marker_hit_count = n_distinct(gene_id),
    .groups = "drop"
  )

arsenic_scores <- tidyr::crossing(
  MAG = unique(mapping$MAG),
  feature_name = arsenic_features
) %>%
  left_join(arsenic_nonfocal_counts, by = c("MAG", "feature_name")) %>%
  mutate(
    nonfocal_marker_hit_count = replace_na(nonfocal_marker_hit_count, 0L),
    feature_component_count = 1L,
    components_detected = as.integer(nonfocal_marker_hit_count > 0),
    reconstruction_score_pct = 100 * components_detected
  ) %>%
  transmute(
    MAG,
    research_category = "Arsenic",
    feature_name,
    source_type = "Independent non-focal arsenic HMM-marker carriage",
    feature_component_count,
    components_detected,
    reconstruction_score_pct
  )

host_scores <- bind_rows(gene_set_scores, arsenic_scores) %>%
  mutate(
    score_ge_50 = reconstruction_score_pct >= 50,
    score_ge_75 = reconstruction_score_pct >= 75
  )
if (!setequal(unique(host_scores$MAG), unique(mapping$MAG))) {
  stop("METABOLIC workbook MAG names do not match the clade classification.")
}
if (nrow(host_scores) != 48L * 40L ||
    n_distinct(host_scores$feature_name) != 40L ||
    anyDuplicated(host_scores[c("MAG", "feature_name")])) {
  stop("Expected 48 MAG x 40 displayed features = 1920 unique rows.")
}

canonical_aioa_raw <- read.table(
  canonical_aioa_tblout,
  comment.char = "#",
  header = FALSE,
  fill = TRUE,
  quote = "",
  stringsAsFactors = FALSE
) %>%
  transmute(
    raw_target_id = V1,
    contig_gene = str_remove(V1, "^[A-Z]{2}\\|"),
    METABOLIC_aioA_full_sequence_E_value = V5,
    METABOLIC_aioA_full_sequence_score = as.numeric(V6),
    METABOLIC_aioA_best_domain_score = as.numeric(V9)
  )
metabolic_aioa_lines <- readLines(metabolic_aioa_tblout, warn = FALSE)
if (!any(str_detect(metabolic_aioa_lines, fixed("-T 800")))) {
  stop("The archived METABOLIC aioA output does not document -T 800.")
}
canonical_aioa_audit <- mapping %>%
  filter(final_clade == "Clade 1") %>%
  select(contig_gene, MAG, habitat, final_clade, final_clade_name) %>%
  left_join(
    classification %>%
      select(
        Sequence_ID,
        joint_screen_Combined_score = Combined_score,
        classification_AioA_score = AioA_score
      ),
    by = c("contig_gene" = "Sequence_ID")
  ) %>%
  left_join(canonical_aioa_raw, by = "contig_gene") %>%
  mutate(
    METABOLIC_reporting_threshold = 800,
    joint_screen_reporting_threshold = 640,
    joint_screen_Combined_640_pass =
      joint_screen_Combined_score >= joint_screen_reporting_threshold,
    METABOLIC_aioA_reported_in_original_tblout =
      METABOLIC_aioA_full_sequence_score >= METABOLIC_reporting_threshold,
    interpretation = case_when(
      METABOLIC_aioA_reported_in_original_tblout ~
        "Passed METABOLIC aioA -T 800; canonical placement independently supported by the joint screen and phylogeny",
      joint_screen_Combined_640_pass ~
        "Below METABOLIC aioA -T 800; retained by independent joint-screen score >=640 and phylogenetic placement, not by METABOLIC detection",
      TRUE ~
        "Below both reporting thresholds; requires independent phylogenetic support"
    )
  )
if (nrow(canonical_aioa_audit) != 3L ||
    sum(canonical_aioa_audit$METABOLIC_aioA_reported_in_original_tblout) != 2L ||
    sum(canonical_aioa_audit$joint_screen_Combined_640_pass) != 3L ||
    any(is.na(canonical_aioa_audit$METABOLIC_aioA_full_sequence_score))) {
  stop("Expected three joint-screen/phylogenetic canonical sequences and two METABOLIC aioA scores >=800.")
}

membership_scores <- host_scores %>%
  inner_join(
    mapping %>% select(MAG, contig_gene, habitat, final_clade,
                       final_clade_name),
    by = "MAG",
    relationship = "many-to-many"
  )

mag_clade_scores <- membership_scores %>%
  group_by(research_category, feature_name, source_type, habitat,
           final_clade, final_clade_name, MAG) %>%
  summarise(
    n_sequence_memberships = n_distinct(contig_gene),
    reconstruction_score_pct = max(reconstruction_score_pct),
    score_ge_50 = any(score_ge_50),
    score_ge_75 = any(score_ge_75),
    .groups = "drop"
  )
if (anyDuplicated(
  mag_clade_scores[c("feature_name", "habitat", "final_clade", "MAG")]
)) {
  stop("A MAG was counted more than once within a feature/habitat/clade cell.")
}

group_summary <- mag_clade_scores %>%
  group_by(research_category, feature_name, source_type, habitat,
           final_clade, final_clade_name) %>%
  summarise(
    n_sequence_memberships = sum(n_sequence_memberships),
    n_unique_MAG = n_distinct(MAG),
    mean_reconstruction_score_pct = mean(reconstruction_score_pct),
    prevalence_score_ge_50_pct = 100 * mean(score_ge_50),
    prevalence_score_ge_75_pct = 100 * mean(score_ge_75),
    inference_allowed = n_distinct(MAG) > 1,
    .groups = "drop"
  )

dual_copy_MAG <- "S1_9-12_bin1"
group_summary_excluding_dual_copy_MAG <- mag_clade_scores %>%
  filter(MAG != dual_copy_MAG) %>%
  group_by(research_category, feature_name, source_type, habitat,
           final_clade, final_clade_name) %>%
  summarise(
    n_unique_MAG_excluding_dual_copy = n_distinct(MAG),
    mean_score_pct_excluding_dual_copy = mean(reconstruction_score_pct),
    prevalence_ge_50_pct_excluding_dual_copy = 100 * mean(score_ge_50),
    prevalence_ge_75_pct_excluding_dual_copy = 100 * mean(score_ge_75),
    .groups = "drop"
  )

dual_copy_MAG_exclusion_sensitivity <- group_summary %>%
  select(
    research_category, feature_name, source_type, habitat, final_clade,
    final_clade_name, n_unique_MAG_all = n_unique_MAG,
    mean_score_pct_all = mean_reconstruction_score_pct,
    prevalence_ge_50_pct_all = prevalence_score_ge_50_pct,
    prevalence_ge_75_pct_all = prevalence_score_ge_75_pct
  ) %>%
  left_join(
    group_summary_excluding_dual_copy_MAG,
    by = c(
      "research_category", "feature_name", "source_type", "habitat",
      "final_clade", "final_clade_name"
    )
  ) %>%
  mutate(
    n_unique_MAG_excluding_dual_copy = replace_na(
      n_unique_MAG_excluding_dual_copy, 0L
    ),
    mean_score_delta_pct =
      mean_score_pct_excluding_dual_copy - mean_score_pct_all,
    prevalence_ge_50_delta_pct =
      prevalence_ge_50_pct_excluding_dual_copy - prevalence_ge_50_pct_all,
    prevalence_ge_75_delta_pct =
      prevalence_ge_75_pct_excluding_dual_copy - prevalence_ge_75_pct_all,
    exclusion_effect = case_when(
      n_unique_MAG_excluding_dual_copy == 0L ~
        "cell becomes empty; original n=1 was descriptive only",
      n_unique_MAG_excluding_dual_copy < n_unique_MAG_all ~
        "cell recalculated after dual-copy MAG exclusion",
      TRUE ~ "unchanged; dual-copy MAG absent from this cell"
    )
  )

dual_copy_MAG_exclusion_summary <-
  dual_copy_MAG_exclusion_sensitivity %>%
  summarise(
    excluded_MAG = dual_copy_MAG,
    excluded_sequence_memberships = sum(mapping$MAG == dual_copy_MAG),
    clades_occupied_by_excluded_MAG = paste(
      sort(unique(as.character(mapping$final_clade[mapping$MAG == dual_copy_MAG]))),
      collapse = ";"
    ),
    excluded_habitat = paste(
      sort(unique(as.character(mapping$habitat[mapping$MAG == dual_copy_MAG]))),
      collapse = ";"
    ),
    feature_clade_habitat_cells_evaluated = n(),
    cells_with_MAG_removed = sum(
      n_unique_MAG_excluding_dual_copy < n_unique_MAG_all
    ),
    cells_becoming_empty = sum(n_unique_MAG_excluding_dual_copy == 0L),
    populated_cells_with_numeric_change = sum(
      n_unique_MAG_excluding_dual_copy > 0L &
        (
          abs(replace_na(mean_score_delta_pct, 0)) > 1e-10 |
          abs(replace_na(prevalence_ge_50_delta_pct, 0)) > 1e-10 |
          abs(replace_na(prevalence_ge_75_delta_pct, 0)) > 1e-10
        )
    ),
    Clade4_IS_affected = any(
      final_clade == "Clade 4" & habitat == "IS" &
        n_unique_MAG_excluding_dual_copy < n_unique_MAG_all
    ),
    interpretation = paste(
      "The dual-copy MAG occupies Clade 1-AS and Clade 3-AS only.",
      "Those n=1 descriptive cells become empty; Clade 4-IS and all other cells are unchanged."
    )
  )

cell_counts <- mapping %>%
  group_by(habitat, final_clade) %>%
  summarise(
    n_sequence_memberships = n_distinct(contig_gene),
    n_unique_MAG = n_distinct(MAG),
    .groups = "drop"
  ) %>%
  complete(
    habitat = factor(c("IS", "AS", "ES", "NS"),
                     levels = c("IS", "AS", "ES", "NS")),
    final_clade = factor(clade_levels, levels = clade_levels),
    fill = list(n_sequence_memberships = 0L, n_unique_MAG = 0L)
  ) %>%
  mutate(
    cell_label = case_when(
      n_unique_MAG == 0 ~ "NA: no members",
      n_unique_MAG == 1 ~ paste0("n=", n_sequence_memberships,
                                " seq / ", n_unique_MAG, " MAG*"),
      TRUE ~ paste0("n=", n_sequence_memberships,
                    " seq / ", n_unique_MAG, " MAGs")
    ),
    x_numeric = as.numeric(habitat)
  )

category_levels <- c(
  "Carbon", "Nitrogen", "Sulfur", "Arsenic", "Vitamin B12", "Motility"
)
feature_order <- bind_rows(
  selected_modules,
  tibble(
    research_category = "Arsenic",
    feature_name = arsenic_features
  )
) %>%
  mutate(research_category = factor(research_category, levels = category_levels)) %>%
  arrange(research_category, feature_name) %>%
  transmute(plot_label = paste(research_category, str_wrap(feature_name, 54),
                               sep = " | ")) %>%
  pull(plot_label)

plot_data <- group_summary %>%
  mutate(
    habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS")),
    final_clade = factor(final_clade, levels = clade_levels),
    plot_label = paste(research_category, str_wrap(feature_name, 54), sep = " | "),
    plot_label = factor(plot_label, levels = rev(unique(feature_order))),
    x_numeric = as.numeric(habitat),
    y_numeric = as.numeric(plot_label)
  )

clade_labels <- setNames(
  paste0(
    clade_levels,
    " (",
    as.integer(observed_clade_counts[clade_levels]),
    " sequences)\n",
    str_wrap(clade_names[clade_levels], width = 31)
  ),
  clade_levels
)

n_plot_features <- length(feature_order)
absent_cells <- cell_counts %>%
  filter(n_unique_MAG == 0) %>%
  mutate(
    xmin = x_numeric - 0.48,
    xmax = x_numeric + 0.48,
    ymin = 0.5,
    ymax = n_plot_features + 0.5
  )
hatch_data <- tidyr::crossing(
  absent_cells,
  hatch_y = seq(0.5, n_plot_features + 0.5, by = 2.5)
) %>%
  mutate(hatch_yend = pmin(hatch_y + 2.0, n_plot_features + 0.5))

make_plot <- function(threshold) {
  prevalence_column <- paste0("prevalence_score_ge_", threshold, "_pct")
  caption_text <- paste(
    "Previous algorithm: each MAG's gene-set coverage is the percentage of distinct predefined genes detected; gene copy number does not increase coverage.",
    "For multi-gene rows, fill is mean gene-set coverage and size is the fraction of MAGs passing the threshold.",
    "For the three arsenic rows, each MAG is binary after excluding all 49 focal K08356 sequences: arrA denotes respiratory arsenate reduction, arsC detoxification/resistance, and arxA/aioA arsenite oxidation.",
    "The legacy denominator is retained exactly; 19 distinct KO identifiers are absent from the current KO result catalog and conservatively depress affected rows (see definition audit).",
    "The >=50% main display is a pre-specified partial gene-set reconstruction threshold; >=75% is a sensitivity analysis, and neither threshold alone establishes a complete pathway or mechanism.",
    "Canonical AioA placement is phylogenetic: 2/3 sequences exceeded METABOLIC aioA -T 800; the third scored 762.0 but passed the independent joint screen (969.5 >=640).",
    "Flagellar Assembly denotes partial 36-KO gene-set reconstruction, not complete motility.",
    paste0("Bubble size is true MAG prevalence at gene-set coverage >=", threshold,
           "%; fill is mean MAG gene-set coverage."),
    "Tiny open circles indicate members are present but the score is zero; gray hatched columns indicate no clade members.",
    "Cell labels show sequence memberships / unique MAGs; * denotes n=1 and descriptive-only evidence.",
    "Habitat uses corrected sample metadata with the group-derep manifest as fallback.",
    "IS, incipient seep; AS, active seep; ES, extinct seep; NS, non-seep background."
  ) %>% str_wrap(width = 190)

  ggplot(plot_data, aes(x = x_numeric, y = y_numeric)) +
    geom_rect(
      data = absent_cells,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "#ECEFF1",
      color = NA
    ) +
    geom_segment(
      data = hatch_data,
      aes(x = xmin, xend = xmax, y = hatch_y, yend = hatch_yend),
      inherit.aes = FALSE,
      color = "#CFD8DC",
      linewidth = 0.35
    ) +
    geom_point(
      shape = 21,
      size = 1.05,
      fill = "white",
      color = "#90A4AE",
      stroke = 0.32
    ) +
    geom_point(
      data = plot_data %>% filter(.data[[prevalence_column]] > 0),
      aes(
        size = .data[[prevalence_column]],
        fill = mean_reconstruction_score_pct
      ),
      shape = 21,
      color = "#263238",
      stroke = 0.28,
      alpha = 0.96
    ) +
    geom_label(
      data = cell_counts %>% filter(n_unique_MAG > 0),
      aes(x = x_numeric, y = n_plot_features + 1.05, label = cell_label),
      inherit.aes = FALSE,
      vjust = 1.15,
      size = 2.25,
      linewidth = 0.15,
      label.padding = unit(0.08, "lines"),
      fill = "white",
      color = "#37474F"
    ) +
    geom_text(
      data = absent_cells,
      aes(x = x_numeric, y = n_plot_features + 1.05, label = "NA: no members"),
      inherit.aes = FALSE,
      vjust = 1.25,
      size = 2.25,
      color = "#607D8B"
    ) +
    facet_grid(
      . ~ final_clade,
      labeller = as_labeller(clade_labels),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = 1:4,
      labels = c("IS", "AS", "ES", "NS"),
      limits = c(0.5, 4.5),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = seq_len(n_plot_features),
      labels = levels(plot_data$plot_label),
      limits = c(0.5, n_plot_features + 1.55),
      expand = expansion(mult = 0)
    ) +
    scale_size_area(
      max_size = 9,
      limits = c(0, 100),
      breaks = c(25, 50, 75, 100),
      name = paste0("MAG prevalence (%)\n(gene-set coverage >=", threshold, "%)")
    ) +
    scale_fill_gradientn(
      colors = c("#F7FCF0", "#7FCDBB", "#2C7FB8", "#253494"),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      name = "Mean gene-set\ncoverage (%)"
    ) +
    labs(
      title = "Metabolic patterns of K08356-containing MAGs across four habitats",
      subtitle = paste0(
        "Previous distinct-gene algorithm; K08356seq49_MAG48; MAG prevalence threshold ",
        threshold,
        "%"
      ),
      x = "Habitat",
      y = "Research-focused metabolic gene set or function",
      caption = caption_text
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_family = "serif", base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10, color = "#37474F"),
      plot.caption = element_text(size = 8.1, hjust = 0, color = "#37474F"),
      strip.background = element_rect(fill = "#D9E2E7", color = "#455A64"),
      strip.text = element_text(face = "bold", size = 9.4, lineheight = 1.05),
      axis.text.x = element_text(face = "bold", size = 9.5),
      axis.text.y = element_text(size = 7.4, color = "#263238"),
      panel.grid.major = element_line(color = "#E2E8EB", linewidth = 0.28),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box.spacing = unit(3, "mm")
    )
}

save_plot <- function(plot, file_stem) {
  ggsave(
    file.path(out_dir, paste0(file_stem, ".png")),
    plot,
    width = 17,
    height = 16.8,
    dpi = 360,
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(file_stem, ".pdf")),
    plot,
    width = 17,
    height = 16.8,
    device = cairo_pdf,
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(file_stem, ".svg")),
    plot,
    width = 17,
    height = 16.8,
    device = grDevices::svg,
    bg = "white"
  )
}

plot_50 <- make_plot(50)
plot_75 <- make_plot(75)
save_plot(
  plot_50,
  "K08356seq49_MAG48_previous_gene_set_coverage50"
)
save_plot(
  plot_75,
  "K08356seq49_MAG48_previous_gene_set_coverage75_sensitivity"
)

write.table(
  mapping,
  file.path(out_dir, "K08356seq49_MAG48_clade_habitat_mapping.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  mapping %>% select(contig_gene, MAG, habitat_in_classification,
                     habitat_in_manifest, habitat_in_sample_metadata,
                     habitat, habitat_assignment_source,
                     classification_habitat_match, manifest_habitat_match),
  file.path(out_dir, "habitat_assignment_QC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  host_scores,
  file.path(out_dir, "K08356seq49_MAG48_metabolic_feature_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  mag_clade_scores,
  file.path(out_dir, "K08356seq49_MAG48_MAG_clade_feature_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  group_summary,
  file.path(out_dir, "habitat_clade_metabolic_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  cell_counts,
  file.path(out_dir, "clade_habitat_cell_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  arsenic_hit_qc,
  file.path(out_dir, "arsenic_marker_hit_focal_exclusion_QC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  canonical_aioa_audit,
  file.path(out_dir, "canonical_AioA_METABOLIC_HMM_score_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  definition_id_summary,
  file.path(out_dir, "gene_set_definition_ID_audit_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  definition_id_detail %>%
    filter(identifier_type == "non-KO" |
             !represented_in_KO_result_catalog),
  file.path(out_dir, "gene_set_non_KO_and_unmatched_definitions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  definition_denominator_impact,
  file.path(out_dir, "gene_set_unavailable_definition_denominator_impact.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  gene_set_denominator_sensitivity_feature,
  file.path(
    out_dir,
    "gene_set_unavailable_KO_threshold_sensitivity_by_feature.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  gene_set_denominator_threshold_flips,
  file.path(
    out_dir,
    "gene_set_unavailable_KO_MAG_threshold_flips.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  dual_copy_MAG_exclusion_sensitivity,
  file.path(out_dir, "dual_copy_MAG_exclusion_sensitivity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  dual_copy_MAG_exclusion_summary,
  file.path(out_dir, "dual_copy_MAG_exclusion_sensitivity_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  plot_data,
  file.path(out_dir, "bubble_plot_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  bind_rows(
    selected_modules %>%
      mutate(source_type = "Previous algorithm: distinct-gene set coverage"),
    tibble(
      research_category = "Arsenic",
      feature_name = arsenic_features,
      source_type = "Independent non-focal arsenic HMM-marker carriage"
    )
  ),
  file.path(out_dir, "selected_metabolic_features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("K08356 sequence memberships:", nrow(mapping), "\n")
cat("Unique MAGs:", n_distinct(mapping$MAG), "\n")
cat("Clade counts:", paste(as.integer(observed_clade_counts), collapse = "/"), "\n")
cat("Classification habitat corrections:",
    sum(!mapping$classification_habitat_match), "\n")
cat("Focal arsenite-oxidation markers excluded:",
    sum(arsenic_hit_qc$feature_name ==
          "Arsenite oxidation (arxA/aioA) marker carriage" &
        arsenic_hit_qc$is_focal_K08356_sequence), "\n")
cat("Non-focal arsenite-oxidation markers retained:",
    sum(arsenic_hit_qc$feature_name ==
          "Arsenite oxidation (arxA/aioA) marker carriage" &
        arsenic_hit_qc$retained_after_focal_exclusion), "\n")
cat("Previous gene sets validated:",
    n_distinct(gene_set_definitions$feature_name), "\n")
cat("Selected metabolic features:", n_distinct(host_scores$feature_name), "\n")
cat("Unavailable-KO threshold flips at 50%:",
    sum(gene_set_denominator_sensitivity_feature$n_MAG_threshold_flip_50), "\n")
cat("Unavailable-KO threshold flips at 75%:",
    sum(gene_set_denominator_sensitivity_feature$n_MAG_threshold_flip_75), "\n")
cat("Clade 4-IS affected by dual-copy MAG exclusion:",
    dual_copy_MAG_exclusion_summary$Clade4_IS_affected, "\n")
cat("Output directory:", out_dir, "\n")

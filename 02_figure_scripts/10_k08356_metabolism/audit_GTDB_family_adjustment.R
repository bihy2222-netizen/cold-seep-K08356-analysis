suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
out_dir <- file.path(script_dir, "results")
taxonomy_dir <- file.path(script_dir, "taxonomy", "raw_gtdb")

mapping_path <- file.path(
  out_dir, "K08356seq49_MAG48_clade_habitat_mapping.tsv"
)
quality_path <- file.path(out_dir, "K08356seq49_MAG48_quality_audit.tsv")
score_path <- file.path(
  out_dir, "K08356seq49_MAG48_metabolic_feature_scores.tsv"
)
summary_paths <- list.files(
  taxonomy_dir,
  pattern = "gtdbtk.*(overall|existing|bac120|ar53).*summary\\.tsv$",
  full.names = TRUE
)

required <- c(mapping_path, quality_path, score_path)
if (any(!file.exists(required)) || !length(summary_paths)) {
  stop("Missing mapping, quality, score, or GTDB summary input.")
}

mapping <- read.delim(mapping_path, check.names = FALSE) %>%
  distinct(MAG, habitat, final_clade, final_clade_name)
quality <- read.delim(quality_path, check.names = FALSE) %>%
  transmute(MAG = MAG_ID, completeness, contamination)
scores <- read.delim(score_path, check.names = FALSE)
expected_MAGs <- sort(unique(mapping$MAG))

read_gtdb <- function(path) {
  read.delim(path, check.names = FALSE, colClasses = "character") %>%
    mutate(
      source_file = basename(path),
      source_priority = if_else(
        str_detect(source_file, fixed("R2111_N500_0-10_bin13")),
        1L,
        2L
      )
    )
}

gtdb_all <- bind_rows(lapply(summary_paths, read_gtdb)) %>%
  filter(user_genome %in% expected_MAGs) %>%
  arrange(user_genome, source_priority) %>%
  distinct(user_genome, .keep_all = TRUE)

missing_MAGs <- setdiff(expected_MAGs, gtdb_all$user_genome)
if (length(missing_MAGs) || nrow(gtdb_all) != 48L ||
    anyDuplicated(gtdb_all$user_genome)) {
  stop(
    "GTDB summaries do not resolve the exact 48 MAGs. Missing: ",
    paste(missing_MAGs, collapse = ", ")
  )
}

rank_names <- c(
  "domain", "phylum", "class", "order", "family", "genus", "species"
)
taxonomy <- gtdb_all %>%
  select(
    MAG = user_genome,
    classification,
    classification_method,
    note,
    warnings,
    msa_percent,
    source_file
  ) %>%
  separate(
    classification,
    into = rank_names,
    sep = ";",
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(
    across(
      all_of(rank_names),
      ~ str_replace(.x, "^[a-z]__", "")
    ),
    across(
      all_of(rank_names),
      ~ na_if(.x, "")
    ),
    family_model = if_else(is.na(family), "Unclassified_family", family)
  ) %>%
  left_join(quality, by = "MAG") %>%
  left_join(
    mapping %>%
      group_by(MAG) %>%
      summarise(
        habitat = paste(sort(unique(habitat)), collapse = ";"),
        clade_memberships = paste(sort(unique(final_clade)), collapse = ";"),
        n_clade_memberships = n_distinct(final_clade),
        .groups = "drop"
      ),
    by = "MAG"
  ) %>%
  arrange(MAG)

family_membership <- mapping %>%
  left_join(
    taxonomy %>% select(MAG, family, family_model),
    by = "MAG"
  )

family_distribution <- family_membership %>%
  group_by(family_model, habitat, final_clade, final_clade_name) %>%
  summarise(n_unique_MAG = n_distinct(MAG), .groups = "drop") %>%
  arrange(desc(n_unique_MAG), family_model, habitat, final_clade)

family_confounding <- family_membership %>%
  group_by(family_model) %>%
  summarise(
    n_unique_MAG = n_distinct(MAG),
    n_habitats = n_distinct(habitat),
    habitats = paste(sort(unique(habitat)), collapse = ";"),
    n_clades = n_distinct(final_clade),
    clades = paste(sort(unique(final_clade)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_unique_MAG), family_model)

model_metadata <- mapping %>%
  filter(MAG != "S1_9-12_bin1") %>%
  distinct(MAG, habitat, final_clade) %>%
  left_join(
    taxonomy %>% select(MAG, family_model),
    by = "MAG"
  ) %>%
  left_join(quality, by = "MAG") %>%
  mutate(
    habitat = factor(habitat),
    final_clade = factor(final_clade),
    family_model = factor(family_model)
  )

full_formula <- ~ completeness + habitat + final_clade + family_model
full_matrix <- model.matrix(full_formula, data = model_metadata)
full_rank <- qr(full_matrix)$rank

term_rank <- function(term) {
  reduced_formula <- switch(
    term,
    completeness = ~ habitat + final_clade + family_model,
    habitat = ~ completeness + final_clade + family_model,
    final_clade = ~ completeness + habitat + family_model,
    family_model = ~ completeness + habitat + final_clade
  )
  reduced_matrix <- model.matrix(reduced_formula, data = model_metadata)
  tibble(
    term,
    full_design_columns = ncol(full_matrix),
    full_design_rank = full_rank,
    reduced_design_rank = qr(reduced_matrix)$rank,
    independently_estimable_df = full_rank - qr(reduced_matrix)$rank
  )
}

model_estimability <- bind_rows(lapply(
  c("completeness", "habitat", "final_clade", "family_model"),
  term_rank
)) %>%
  mutate(
    n_MAG = nrow(model_metadata),
    residual_df = n_MAG - full_design_rank,
    interpretation = case_when(
      independently_estimable_df == 0L ~
        "Term is completely aliased with the other covariates",
      term == "family_model" ~
        "Some family structure is estimable, but sparse families limit inference",
      TRUE ~
        "Only the reported degrees of freedom are separable from the other covariates"
    )
  )

term_test <- function(df, term) {
  if (n_distinct(df$reconstruction_score_pct) < 2L) {
    return(tibble(term, df = NA_real_, F_value = NA_real_, p_value = NA_real_))
  }
  fit <- lm(
    reconstruction_score_pct ~
      completeness + habitat + final_clade + family_model,
    data = df
  )
  drop <- suppressWarnings(drop1(fit, test = "F"))
  if (!term %in% rownames(drop)) {
    return(tibble(term, df = NA_real_, F_value = NA_real_, p_value = NA_real_))
  }
  tibble(
    term,
    df = unname(drop[term, "Df"]),
    F_value = unname(drop[term, "F value"]),
    p_value = unname(drop[term, "Pr(>F)"])
  )
}

family_adjusted_tests <- scores %>%
  inner_join(model_metadata, by = "MAG") %>%
  group_by(research_category, feature_name, source_type) %>%
  group_modify(function(.x, .y) {
    bind_rows(lapply(
      c("completeness", "habitat", "final_clade", "family_model"),
      function(term) term_test(.x, term)
    ))
  }) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    n_MAG = nrow(model_metadata),
    model_note = paste(
      "Descriptive OLS after excluding dual-clade S1_9-12_bin1;",
      "interpret only terms with independently estimable degrees of freedom"
    )
  )

locked_features <- c(
  "Arsenate detoxification/resistance (arsC) marker carriage",
  "Respiratory arsenate reduction (arrA) marker carriage",
  "Arsenite oxidation (arxA/aioA) marker carriage",
  "Cobalamin biosynthesis, cobinamide => cobalamin"
)
rhodo_locked <- scores %>%
  inner_join(model_metadata, by = "MAG") %>%
  filter(
    habitat == "IS",
    family_model == "Rhodobacteraceae",
    final_clade %in% c("Clade 3", "Clade 4"),
    feature_name %in% locked_features
  )

compare_rhodo <- function(df) {
  clade_counts <- table(droplevels(df$final_clade))
  invariant <- n_distinct(df$reconstruction_score_pct) < 2L
  fisher_p <- NA_real_
  wilcoxon_p <- NA_real_
  clade_lm_p <- NA_real_
  if (!invariant && all(c("Clade 3", "Clade 4") %in% names(clade_counts))) {
    if (all(df$reconstruction_score_pct %in% c(0, 100))) {
      fisher_p <- fisher.test(
        table(df$final_clade, df$reconstruction_score_pct > 0)
      )$p.value
    } else {
      wilcoxon_p <- suppressWarnings(
        wilcox.test(reconstruction_score_pct ~ final_clade, data = df)$p.value
      )
    }
    fit <- lm(
      reconstruction_score_pct ~ completeness + final_clade,
      data = df
    )
    coefficient_row <- grep(
      "^final_clade",
      rownames(summary(fit)$coefficients),
      value = TRUE
    )
    if (length(coefficient_row)) {
      clade_lm_p <- summary(fit)$coefficients[
        coefficient_row[1], "Pr(>|t|)"
      ]
    }
  }
  tibble(
    n_Rhodobacteraceae_IS = n_distinct(df$MAG),
    n_Clade3 = sum(df$final_clade == "Clade 3"),
    n_Clade4 = sum(df$final_clade == "Clade 4"),
    mean_Clade3 = mean(
      df$reconstruction_score_pct[df$final_clade == "Clade 3"]
    ),
    mean_Clade4 = mean(
      df$reconstruction_score_pct[df$final_clade == "Clade 4"]
    ),
    fisher_p = fisher_p,
    wilcoxon_p = wilcoxon_p,
    completeness_adjusted_clade_lm_p = clade_lm_p,
    interpretation = if_else(
      invariant,
      "Outcome invariant; no clade test is possible",
      paste(
        "Within-family, within-habitat sensitivity only;",
        "small Clade 3 sample and non-independent pathway definitions remain"
      )
    )
  )
}

rhodo_locked_comparison <- rhodo_locked %>%
  group_by(research_category, feature_name, source_type) %>%
  group_modify(~ compare_rhodo(.x)) %>%
  ungroup() %>%
  mutate(
    fisher_BH = p.adjust(fisher_p, method = "BH"),
    wilcoxon_BH = p.adjust(wilcoxon_p, method = "BH"),
    completeness_adjusted_clade_lm_BH = p.adjust(
      completeness_adjusted_clade_lm_p,
      method = "BH"
    )
  )

write.table(
  taxonomy,
  file.path(out_dir, "K08356_MAG48_GTDB_r226_taxonomy.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  family_distribution,
  file.path(out_dir, "GTDB_family_by_habitat_clade.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  family_confounding,
  file.path(out_dir, "GTDB_family_confounding_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  model_estimability,
  file.path(out_dir, "GTDB_family_model_estimability.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  family_adjusted_tests,
  file.path(out_dir, "GTDB_family_adjusted_feature_tests.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  rhodo_locked_comparison,
  file.path(out_dir, "Rhodobacteraceae_IS_Clade3_vs_Clade4_locked_features.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Exact GTDB taxonomy rows:", nrow(taxonomy), "\n")
cat("Distinct assigned families:", n_distinct(taxonomy$family, na.rm = TRUE), "\n")
cat(
  "Families spanning more than one habitat:",
  sum(family_confounding$n_habitats > 1L), "\n"
)
cat(
  "Full family-adjusted design rank:", full_rank, "/", ncol(full_matrix),
  "; residual df:", nrow(model_metadata) - full_rank, "\n"
)
cat(
  "Rhodobacteraceae IS locked comparisons:",
  nrow(rhodo_locked_comparison), "\n"
)

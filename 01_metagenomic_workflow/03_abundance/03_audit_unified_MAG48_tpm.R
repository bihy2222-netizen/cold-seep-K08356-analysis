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
repo_root <- normalizePath(file.path(script_dir, "../.."))
run_dir <- Sys.getenv(
  "RUN_DIR",
  "/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/K08356_four_group_annotation/unified_MAG48_coverm_tpm_20260802"
)
mapping_path <- Sys.getenv(
  "K08356_MAPPING_TSV",
  file.path(
    repo_root,
    "02_figure_scripts/10_k08356_metabolism/results",
    "K08356seq49_MAG48_clade_habitat_mapping.tsv"
  )
)

sample_path <- file.path(run_dir, "manifests", "sample56_reads_habitat.tsv")
reference_path <- file.path(run_dir, "manifests", "MAG48_reference_manifest.tsv")
tpm_dir <- file.path(run_dir, "tpm_per_sample")
result_dir <- file.path(run_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(sample_path, reference_path, mapping_path)
if (any(!file.exists(required))) {
  stop("Missing audit input: ", paste(required[!file.exists(required)], collapse = ", "))
}

samples <- read.delim(sample_path, check.names = FALSE)
reference <- read.delim(reference_path, check.names = FALSE)
mapping <- read.delim(mapping_path, check.names = FALSE)

if (nrow(samples) != 56L || n_distinct(samples$sample) != 56L) {
  stop("Expected 56 unique samples.")
}
if (nrow(reference) != 48L || n_distinct(reference$MAG) != 48L) {
  stop("Expected 48 unique reference MAGs.")
}
if (!setequal(reference$MAG, unique(mapping$MAG))) {
  stop("Common reference does not match the exact current 48-MAG mapping.")
}

read_sample_tpm <- function(sample, habitat) {
  path <- file.path(tpm_dir, paste0(sample, "_MAG48_tpm.tsv"))
  if (!file.exists(path)) {
    stop("Missing TPM file: ", path)
  }
  raw <- read.delim(path, check.names = FALSE)
  if (ncol(raw) != 2L) {
    stop("Expected two CoverM columns in ", path)
  }
  names(raw) <- c("MAG", "TPM")
  raw <- raw %>%
    mutate(
      MAG = sub("\\.fasta$", "", basename(MAG)),
      TPM = as.numeric(TPM)
    )
  if (nrow(raw) != 48L || n_distinct(raw$MAG) != 48L ||
      !setequal(raw$MAG, reference$MAG) || anyNA(raw$TPM) ||
      any(raw$TPM < 0)) {
    stop("TPM file failed exact 48-MAG validation: ", path)
  }
  raw %>% mutate(sample = sample, habitat = habitat, .before = 1)
}

long_tpm <- bind_rows(Map(
  read_sample_tpm,
  samples$sample,
  samples$habitat
))
if (nrow(long_tpm) != 56L * 48L ||
    anyDuplicated(long_tpm[c("sample", "MAG")])) {
  stop("Expected 2,688 unique sample-MAG abundance rows.")
}

tpm_matrix <- long_tpm %>%
  select(sample, habitat, MAG, TPM) %>%
  pivot_wider(names_from = MAG, values_from = TPM) %>%
  arrange(factor(habitat, levels = c("IS", "AS", "ES", "NS")), sample)

MAG_membership <- mapping %>%
  distinct(MAG, final_clade, final_clade_name) %>%
  group_by(MAG) %>%
  mutate(
    n_clade_memberships = n_distinct(final_clade),
    multi_clade_host = n_clade_memberships > 1L
  ) %>%
  ungroup()

nonexclusive_clade_tpm <- long_tpm %>%
  inner_join(MAG_membership, by = "MAG", relationship = "many-to-many") %>%
  group_by(sample, habitat, final_clade, final_clade_name) %>%
  summarise(
    n_unique_MAG = n_distinct(MAG),
    sum_MAG_TPM = sum(TPM),
    mean_MAG_TPM = mean(TPM),
    .groups = "drop"
  ) %>%
  mutate(
    aggregation_note = paste(
      "Clade categories are non-exclusive because S1_9-12_bin1 belongs to",
      "Clade 1 and Clade 3; do not sum clades as total host abundance"
    )
  )

exclusive_membership <- MAG_membership %>%
  group_by(MAG) %>%
  summarise(
    host_category = if_else(
      any(multi_clade_host),
      "Multi-clade host",
      first(final_clade)
    ),
    .groups = "drop"
  )
exclusive_category_tpm <- long_tpm %>%
  inner_join(exclusive_membership, by = "MAG") %>%
  group_by(sample, habitat, host_category) %>%
  summarise(
    n_unique_MAG = n_distinct(MAG),
    sum_MAG_TPM = sum(TPM),
    mean_MAG_TPM = mean(TPM),
    .groups = "drop"
  )

exclude_dual_tpm <- long_tpm %>%
  filter(MAG != "S1_9-12_bin1") %>%
  inner_join(
    MAG_membership %>% filter(!multi_clade_host),
    by = "MAG"
  ) %>%
  group_by(sample, habitat, final_clade, final_clade_name) %>%
  summarise(
    n_unique_MAG = n_distinct(MAG),
    sum_MAG_TPM = sum(TPM),
    mean_MAG_TPM = mean(TPM),
    .groups = "drop"
  )

write.table(
  long_tpm,
  file.path(result_dir, "sample56_MAG48_TPM_long.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  tpm_matrix,
  file.path(result_dir, "sample56_by_MAG48_TPM_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  nonexclusive_clade_tpm,
  file.path(result_dir, "sample_habitat_clade_MAG_TPM_nonexclusive.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  exclusive_category_tpm,
  file.path(result_dir, "sample_habitat_host_category_MAG_TPM_exclusive.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  exclude_dual_tpm,
  file.path(result_dir, "sample_habitat_clade_MAG_TPM_excluding_dual_host.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Verified samples:", n_distinct(long_tpm$sample), "\n")
cat("Verified MAGs:", n_distinct(long_tpm$MAG), "\n")
cat("Unique sample-MAG rows:", nrow(long_tpm), "\n")
cat("Dual-clade MAG handled in non-exclusive, exclusive, and exclusion outputs.\n")

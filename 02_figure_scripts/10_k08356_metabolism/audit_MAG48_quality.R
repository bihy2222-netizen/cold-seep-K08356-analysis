#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) paths[1] else hit[1]
}

manifest_path <- Sys.getenv(
  "MAG49_MANIFEST_TSV",
  first_existing(c(
    file.path(script_dir, "MAG49_input_manifest.tsv"),
    file.path(script_dir, "input", "MAG49_input_manifest.tsv"),
    file.path(dirname(dirname(script_dir)), "01_metagenomic_workflow",
              "04_metabolism", "METABOLIC_group_derep49", "manifests",
              "MAG49_input_manifest.tsv")
  ))
)
classification_path <- Sys.getenv(
  "K08356_CLADE_TSV",
  first_existing(c(
    file.path(script_dir, "input", "four_clade_classification_49.tsv"),
    file.path(script_dir, "..", "four_clade_classification_49.tsv")
  ))
)
stats_dir <- Sys.getenv(
  "QUALITY_STATS_DIR",
  file.path(script_dir, "quality", "raw_stats")
)
out_dir <- file.path(script_dir, "results")

manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
classification <- read.delim(
  classification_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

manifest$sample_id <- sub("_bin[0-9]+$", "", manifest$bin_id)
manifest$stats_bin <- sub("^.*_bin([0-9]+)$", "bin.\\1", manifest$bin_id)
manifest$stats_path <- file.path(
  stats_dir,
  paste0(manifest$group, "__", manifest$sample_id, ".stats")
)

read_quality <- function(i) {
  row <- manifest[i, ]
  if (!file.exists(row$stats_path)) {
    stop("Missing quality file: ", row$stats_path)
  }
  stats <- read.delim(row$stats_path, check.names = FALSE, stringsAsFactors = FALSE)
  hit <- stats[stats$bin == row$stats_bin, , drop = FALSE]
  if (nrow(hit) != 1L) {
    stop("Expected one quality row for ", row$bin_id, "; found ", nrow(hit))
  }
  data.frame(
    MAG_ID = row$bin_id,
    Habitat = row$group,
    completeness = hit$completeness,
    contamination = hit$contamination,
    GC = hit$GC,
    checkm_lineage = hit$lineage,
    N50 = hit$N50,
    genome_size_bp = hit$size,
    binner = hit$binner,
    quality_source = basename(row$stats_path),
    stringsAsFactors = FALSE
  )
}

quality <- do.call(rbind, lapply(seq_len(nrow(manifest)), read_quality))
clade_by_mag <- aggregate(
  classification$clade_display,
  by = list(MAG_ID = classification$MAG_ID),
  FUN = function(x) paste(sort(unique(x)), collapse = "; ")
)
names(clade_by_mag)[2] <- "clade_memberships"
quality <- merge(quality, clade_by_mag, by = "MAG_ID", all.x = TRUE, sort = FALSE)
quality <- quality[match(manifest$bin_id, quality$MAG_ID), ]

if (nrow(quality) != 48L || anyDuplicated(quality$MAG_ID)) {
  stop("Quality audit did not resolve exactly 48 unique MAGs")
}

quality$meets_MQ_50_10 <- quality$completeness >= 50 & quality$contamination <= 10
write.table(
  quality,
  file.path(out_dir, "K08356seq49_MAG48_quality_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary <- aggregate(
  cbind(completeness, contamination) ~ Habitat,
  data = quality,
  FUN = function(x) round(c(min = min(x), median = median(x), max = max(x)), 2)
)
summary <- data.frame(
  Habitat = summary$Habitat,
  n_MAG = as.integer(table(quality$Habitat)[summary$Habitat]),
  completeness_min = summary$completeness[, "min"],
  completeness_median = summary$completeness[, "median"],
  completeness_max = summary$completeness[, "max"],
  contamination_min = summary$contamination[, "min"],
  contamination_median = summary$contamination[, "median"],
  contamination_max = summary$contamination[, "max"]
)
write.table(
  summary,
  file.path(out_dir, "K08356seq49_MAG48_quality_by_habitat.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mag_clade <- unique(classification[c("MAG_ID", "clade_display")])
quality_clade <- merge(
  quality[c("MAG_ID", "completeness", "contamination")],
  mag_clade,
  by = "MAG_ID",
  all.y = TRUE
)
clade_summary <- aggregate(
  cbind(completeness, contamination) ~ clade_display,
  data = quality_clade,
  FUN = function(x) round(c(min = min(x), median = median(x), max = max(x)), 2)
)
clade_summary <- data.frame(
  clade_display = clade_summary$clade_display,
  n_MAG = as.integer(table(quality_clade$clade_display)[clade_summary$clade_display]),
  completeness_min = clade_summary$completeness[, "min"],
  completeness_median = clade_summary$completeness[, "median"],
  completeness_max = clade_summary$completeness[, "max"],
  contamination_min = clade_summary$contamination[, "min"],
  contamination_median = clade_summary$contamination[, "median"],
  contamination_max = clade_summary$contamination[, "max"]
)
write.table(
  clade_summary,
  file.path(out_dir, "K08356seq49_MAG48_quality_by_clade.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Quality rows:", nrow(quality), "\n")
cat("Completeness range:", range(quality$completeness), "\n")
cat("Contamination range:", range(quality$contamination), "\n")
cat("All pass >=50 completeness and <=10 contamination:", all(quality$meets_MQ_50_10), "\n")

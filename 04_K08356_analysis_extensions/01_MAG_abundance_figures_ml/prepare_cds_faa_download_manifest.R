suppressPackageStartupMessages({
  library(dplyr)
})

branch_map_file <- "/Users/catherine/Downloads/MAG-bin-tpm/derepMAG_abundance_stats/k08356_branch_feature_map.tsv"
out_dir <- "/Users/catherine/Downloads/MAG-bin-tpm/K08356_MAG_gene_island_arrows_20260708"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

branch_map <- read.delim(branch_map_file, check.names = FALSE)

manifest <- branch_map %>%
  mutate(
    cds_faa = paste0(MAG, ".cds.faa"),
    hpc_source_path = paste0("/home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/CDS/", cds_faa),
    local_target_path = paste0(out_dir, "/cds_faa/", cds_faa)
  ) %>%
  distinct(MAG, cds_faa, hpc_source_path, local_target_path, Feature, Corrected_habitat, Original_MAG_ID) %>%
  arrange(Feature, MAG)

write.csv(manifest, file.path(out_dir, "required_cds_faa_manifest.csv"), row.names = FALSE)
unique_files <- manifest %>%
  distinct(cds_faa, hpc_source_path) %>%
  arrange(cds_faa)

writeLines(unique_files$cds_faa, file.path(out_dir, "required_cds_faa_filenames.txt"))
writeLines(unique_files$hpc_source_path, file.path(out_dir, "required_cds_faa_hpc_paths.txt"))

copy_script <- c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  "# Run this on the supercomputer. It creates a tarball in the current directory only,",
  "# and does not write anything into /home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/CDS.",
  "SRC='/home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/CDS'",
  "OUT='K08356_required_cds_faa.tar.gz'",
  "TMP_LIST='K08356_required_cds_faa_files.txt'",
  "",
  "cat > \"$TMP_LIST\" <<'EOF'",
  unique_files$cds_faa,
  "EOF",
  "",
  "tar -C \"$SRC\" -czf \"$OUT\" -T \"$TMP_LIST\"",
  "echo \"Wrote $OUT\""
)

writeLines(copy_script, file.path(out_dir, "make_K08356_required_cds_faa_tarball_on_hpc.sh"))

message("Wrote manifest to: ", out_dir)

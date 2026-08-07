suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

branch_map_file <- "/Users/catherine/Downloads/MAG-bin-tpm/derepMAG_abundance_stats/k08356_branch_feature_map.tsv"
out_dir <- "/Users/catherine/Downloads/MAG-bin-tpm/K08356_MAG_gene_island_arrows_20260708"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

remote_dir <- "/home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/CDS/cat_result/cat_others"

sample_from_protein <- function(x) {
  sub("_bin[0-9]+-k141_.*$", "", x)
}

branch_map <- read.delim(branch_map_file, check.names = FALSE) %>%
  mutate(
    Sample = sample_from_protein(Original_MAG_ID),
    remote_cds_faa = paste0(remote_dir, "/", Sample, ".cdhit.kofamscan.best.contig.cds.faa"),
    remote_gff = paste0(remote_dir, "/", Sample, ".cdhit.kofamscan.best.contig.gff"),
    local_cds_faa = paste0(out_dir, "/sample_cds_gff/", Sample, ".cdhit.kofamscan.best.contig.cds.faa"),
    local_gff = paste0(out_dir, "/sample_cds_gff/", Sample, ".cdhit.kofamscan.best.contig.gff")
  )

sample_manifest <- branch_map %>%
  distinct(Sample, remote_cds_faa, remote_gff, local_cds_faa, local_gff) %>%
  arrange(Sample)

target_manifest <- branch_map %>%
  select(MAG, Feature, Original_MAG_ID, Corrected_habitat, Corrected_sample_group,
         Sample, remote_cds_faa, remote_gff, local_cds_faa, local_gff) %>%
  arrange(Feature, Original_MAG_ID)

write.csv(sample_manifest, file.path(out_dir, "required_sample_cds_gff_manifest.csv"), row.names = FALSE)
write.csv(target_manifest, file.path(out_dir, "target_to_sample_cds_gff_manifest.csv"), row.names = FALSE)
writeLines(unique(c(sample_manifest$remote_cds_faa, sample_manifest$remote_gff)),
           file.path(out_dir, "required_sample_cds_gff_remote_paths.txt"))

message("Wrote sample-level CDS/GFF manifest to: ", out_dir)

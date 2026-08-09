suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

base_out <- "/Users/catherine/Downloads/MAG-bin-tpm/K08356_gene_and_MAG_dual_evidence_20260708"
fig_out <- file.path(base_out, "overall_derepMAG_final_figures")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

habitat_levels <- c("IS", "AS", "ES", "NS")

sample_to_habitat <- function(sample) {
  case_when(
    grepl("^SY", sample) ~ "IS",
    grepl("^SQ_", sample) ~ "AS",
    grepl("^(S1_|S2_|S4_)", sample) ~ "AS",
    grepl("^(C1_|C2_|C3_|S13_|S14_|S15_|ES_)", sample) ~ "ES",
    grepl("^(S3_|R2111_|NS_)", sample) ~ "NS",
    TRUE ~ NA_character_
  )
}

sample_group <- function(sample) {
  case_when(
    grepl("^SY365", sample) ~ "SY365",
    grepl("^SY366YB", sample) ~ "SY366YB",
    grepl("^SY366YW", sample) ~ "SY366YW",
    grepl("^SY368", sample) ~ "SY368",
    grepl("^SY456", sample) ~ "SY456",
    grepl("^SY457", sample) ~ "SY457",
    grepl("^SY459", sample) ~ "SY459",
    TRUE ~ sub("_.*$", "", sample)
  )
}

depth_start <- function(sample) {
  x <- sub("^.*_([0-9]+)-[0-9]+$", "\\1", sample)
  suppressWarnings(as.numeric(ifelse(x == sample, 0, x)))
}

branch_map <- read.delim("derepMAG_abundance_stats/k08356_branch_feature_map.tsv", check.names = FALSE) %>%
  filter(Feature == "canonical aioA-associated") %>%
  select(MAG, Original_MAG_ID, Feature, Corrected_habitat, Corrected_sample_group) %>%
  distinct()

mag_matrix <- read.delim("merged_tpm_matrix.txt", check.names = FALSE)
colnames(mag_matrix)[1] <- "MAG"
sample_cols <- setdiff(colnames(mag_matrix), "MAG")

canonical_long <- mag_matrix %>%
  filter(MAG %in% branch_map$MAG) %>%
  pivot_longer(all_of(sample_cols), names_to = "Sample", values_to = "TPM") %>%
  left_join(branch_map, by = "MAG") %>%
  mutate(
    Sample_habitat = factor(sample_to_habitat(Sample), levels = habitat_levels),
    MAG_source_habitat = Corrected_habitat,
    Sample_group = sample_group(Sample),
    depth = depth_start(Sample)
  ) %>%
  filter(!is.na(Sample_habitat))

canonical_is <- canonical_long %>%
  filter(Sample_habitat == "IS", TPM > 0) %>%
  arrange(Sample_group, depth, Sample, desc(TPM))

canonical_summary <- canonical_is %>%
  group_by(MAG, Original_MAG_ID, MAG_source_habitat, Corrected_sample_group) %>%
  summarise(
    n_IS_samples_nonzero = n_distinct(Sample),
    sum_IS_TPM = sum(TPM),
    max_IS_TPM = max(TPM),
    max_IS_sample = Sample[which.max(TPM)],
    .groups = "drop"
  ) %>%
  arrange(desc(sum_IS_TPM))

write_csv(canonical_is, file.path(fig_out, "canonical_aioA_IS_positive_TPM_by_MAG.csv"))
write_csv(canonical_summary, file.path(fig_out, "canonical_aioA_IS_positive_TPM_summary_by_MAG.csv"))

sample_order <- canonical_long %>%
  filter(Sample_habitat == "IS") %>%
  distinct(Sample, Sample_group, depth) %>%
  arrange(Sample_group, depth, Sample) %>%
  pull(Sample)

diag_fig <- canonical_long %>%
  filter(Sample_habitat == "IS", TPM > 0) %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    MAG_label = paste0(MAG, " (", MAG_source_habitat, "-origin)")
  ) %>%
  ggplot(aes(Sample, TPM, fill = MAG_label)) +
  geom_col(width = 0.82, color = "grey35", linewidth = 0.1) +
  scale_y_continuous(labels = scales::label_scientific()) +
  labs(
    title = "Diagnostic: canonical aioA-associated MAG TPM in IS samples",
    subtitle = "Orange in the stacked branch plot comes from cross-sample recruitment of canonical-associated host MAGs",
    x = "IS sample",
    y = "Host MAG TPM",
    fill = "canonical-associated derepMAG"
  ) +
  theme_classic(base_family = "Times", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

ggsave(file.path(fig_out, "canonical_aioA_IS_contribution_diagnostic.png"), diag_fig, width = 8.2, height = 6.2, dpi = 450)
ggsave(file.path(fig_out, "canonical_aioA_IS_contribution_diagnostic.pdf"), diag_fig, width = 8.2, height = 6.2)

note <- c(
  "# canonical aioA-associated stacked-bar diagnostic",
  "",
  "The branch stacked bar is not an origin-only occurrence plot. It shows read-recruited MAG TPM across all samples.",
  "",
  "For canonical aioA-associated K08356, the branch map contains three derepMAG entries:",
  paste0("- ", branch_map$MAG, " | source habitat: ", branch_map$Corrected_habitat, " | original gene id: ", branch_map$Original_MAG_ID),
  "",
  "In IS samples, nonzero canonical-associated TPM is contributed mainly by:",
  paste0("- ", canonical_summary$MAG, " | source habitat: ", canonical_summary$MAG_source_habitat,
         " | nonzero IS samples: ", canonical_summary$n_IS_samples_nonzero,
         " | sum IS TPM: ", round(canonical_summary$sum_IS_TPM, 3),
         " | max IS sample: ", canonical_summary$max_IS_sample,
         " | max TPM: ", round(canonical_summary$max_IS_TPM, 3)),
  "",
  "Interpretation: the IS-origin canonical MAG is present mainly at SY366YW-4-8, but the NS-origin R2111_S300 canonical derepMAG also recruits reads in multiple IS samples. Therefore orange segments in several IS bars are abundance/recruitment signals from the MAG TPM matrix, not additional IS-origin canonical bins.",
  "",
  "Recommended wording: `canonical aioA-associated host MAG abundance was not significantly different among habitats (Kruskal-Wallis p = 0.877); although only one canonical-associated MAG was recovered from IS, low-level read recruitment to canonical-associated derepMAG representatives was observed across multiple IS samples.`"
)
writeLines(note, file.path(fig_out, "README_canonical_aioA_stacked_bar_diagnostic.md"))

message("Wrote canonical diagnostic files to: ", fig_out)

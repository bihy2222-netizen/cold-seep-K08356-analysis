suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

base <- "/Users/catherine/Downloads/aioA 蛋白序列建树/tree11_no_uncertain_DryadAioAIdrAfinal"
outdir <- "/Users/catherine/Downloads/MAG-bin-tpm/k08356_gene_level_contig_style_report"
figdir <- file.path(outdir, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

data <- fromJSON(file.path(base, "tree11_integrated_context_data.json"), flatten = TRUE)
master <- as_tibble(data$master)

habitat_levels <- c("IS", "AS", "ES", "NS")
hab_cols <- c(IS = "#2ca25f", AS = "#8e63c7", ES = "#e39d25", NS = "#3b78c8")
clade_levels <- c("idrA", "canonical_aioA", "aioA_like", "aioA_IdrA_intermediate")
clade_labels <- c(
  idrA = "IdrA-associated",
  canonical_aioA = "canonical aioA-associated",
  aioA_like = "aioA-like-associated",
  aioA_IdrA_intermediate = "unknown AioA-like / uncertain DMSOR"
)

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 3))))
}

gene_long <- master %>%
  distinct(clade_short, clade, habitat, sample, sample_K08356_TPM) %>%
  mutate(
    clade_short = factor(clade_short, levels = clade_levels),
    clade_label = factor(clade_labels[as.character(clade_short)], levels = unname(clade_labels)),
    habitat = factor(habitat, levels = habitat_levels),
    TPM = as.numeric(sample_K08356_TPM),
    TPM = ifelse(is.na(TPM), 0, TPM),
    logTPM = log10(TPM + 1)
  ) %>%
  filter(!is.na(habitat), !is.na(clade_short))

sample_counts <- gene_long %>%
  distinct(sample, habitat) %>%
  count(habitat, name = "n_samples") %>%
  mutate(habitat_label = paste0(habitat, "\n(n=", n_samples, ")"))

gene_long <- gene_long %>%
  left_join(sample_counts %>% select(habitat, habitat_label), by = "habitat") %>%
  mutate(habitat_label = factor(habitat_label, levels = sample_counts$habitat_label[match(habitat_levels, sample_counts$habitat)]))

kw_by_clade <- gene_long %>%
  group_by(clade_short, clade_label) %>%
  group_modify(~ {
    if (nrow(.x) < 3 || length(unique(.x$habitat)) < 2) {
      return(tibble(statistic = NA_real_, df = NA_real_, p_value = NA_real_, n = nrow(.x), tested_value = "log10(TPM+1)"))
    }
    kt <- kruskal.test(logTPM ~ habitat, data = .x)
    tibble(statistic = unname(kt$statistic), df = unname(kt$parameter), p_value = kt$p.value, n = nrow(.x), tested_value = "log10(TPM+1)")
  }) %>%
  ungroup()

pairwise_rows <- list()
for (cl in levels(gene_long$clade_short)) {
  dat <- gene_long %>% filter(clade_short == cl)
  if (nrow(dat) >= 3 && length(unique(dat$habitat)) >= 2) {
    pw <- pairwise.wilcox.test(dat$logTPM, dat$habitat, p.adjust.method = "BH", exact = FALSE)
    mat <- pw$p.value
    if (!is.null(mat)) {
      for (r in rownames(mat)) {
        for (c in colnames(mat)) {
          p <- mat[r, c]
          if (!is.na(p)) {
            pairwise_rows[[length(pairwise_rows) + 1]] <- tibble(
              clade_short = cl,
              clade_label = clade_labels[[cl]],
              group1 = c,
              group2 = r,
              p_adj_BH = p,
              method = "Pairwise Wilcoxon rank-sum test on log10(TPM+1), BH-adjusted"
            )
          }
        }
      }
    }
  }
}
pairwise <- bind_rows(pairwise_rows)

summary_stats <- gene_long %>%
  group_by(clade_short, clade_label, habitat) %>%
  summarise(
    n_sample_clade_records = n(),
    n_detected_TPM_gt_0 = sum(TPM > 0),
    detection_rate = n_detected_TPM_gt_0 / n_sample_clade_records,
    mean_TPM = mean(TPM),
    median_TPM = median(TPM),
    max_TPM = max(TPM),
    mean_log10_TPM_plus_1 = mean(logTPM),
    median_log10_TPM_plus_1 = median(logTPM),
    .groups = "drop"
  )

kw_labels <- kw_by_clade %>%
  mutate(label = paste0("KW ", fmt_p(p_value)))

p_branch <- ggplot(gene_long, aes(x = habitat, y = logTPM, fill = habitat)) +
  geom_boxplot(width = 0.58, alpha = 0.72, color = "grey25", outlier.shape = NA, linewidth = 0.45) +
  geom_jitter(aes(color = habitat), width = 0.12, size = 1.9, alpha = 0.9, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.8, fill = "white", color = "black") +
  geom_text(
    data = kw_labels,
    aes(x = 2.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.35,
    size = 3.2
  ) +
  facet_wrap(~ clade_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = hab_cols, drop = FALSE) +
  scale_color_manual(values = hab_cols, drop = FALSE) +
  labs(
    title = "Gene/contig-level K08356 TPM by Tree11 branch",
    subtitle = "Abundance uses sample-level K08356 gene TPM, not host MAG TPM",
    x = "Habitat",
    y = "K08356 abundance, log10(TPM + 1)"
  ) +
  theme_classic(base_family = "Times", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    strip.background = element_rect(fill = "#EEF2F7", color = "#D8DEE9"),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "grey15"),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    plot.margin = margin(12, 14, 12, 12)
  )

ggsave(file.path(figdir, "K08356_gene_TPM_by_branch_boxplots.png"), p_branch, width = 8, height = 8, dpi = 450)
ggsave(file.path(figdir, "K08356_gene_TPM_by_branch_boxplots.pdf"), p_branch, width = 8, height = 8)
svg(file.path(figdir, "K08356_gene_TPM_by_branch_boxplots.svg"), width = 8, height = 8, family = "Times")
print(p_branch)
dev.off()

total_gene <- gene_long %>%
  group_by(sample, habitat) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(logTPM = log10(TPM + 1)) %>%
  left_join(sample_counts %>% select(habitat, habitat_label), by = "habitat") %>%
  mutate(habitat_label = factor(habitat_label, levels = sample_counts$habitat_label[match(habitat_levels, sample_counts$habitat)]))

kw_total <- kruskal.test(logTPM ~ habitat, data = total_gene)

p_total <- ggplot(total_gene, aes(x = habitat, y = logTPM, fill = habitat)) +
  geom_boxplot(width = 0.58, alpha = 0.72, color = "grey25", outlier.shape = NA, linewidth = 0.45) +
  geom_jitter(aes(color = habitat), width = 0.12, size = 2.2, alpha = 0.9, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.0, fill = "white", color = "black") +
  scale_fill_manual(values = hab_cols, drop = FALSE) +
  scale_color_manual(values = hab_cols, drop = FALSE) +
  labs(
    title = "Gene/contig-level total K08356 TPM across habitats",
    subtitle = paste0("Kruskal-Wallis on log10(TPM+1): ", fmt_p(kw_total$p.value)),
    x = "Habitat",
    y = "K08356 abundance, log10(TPM + 1)"
  ) +
  theme_classic(base_family = "Times", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(color = "grey15"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    plot.margin = margin(12, 14, 12, 12)
  )

ggsave(file.path(figdir, "K08356_gene_total_TPM_boxplot.png"), p_total, width = 6.4, height = 6.0, dpi = 450)
ggsave(file.path(figdir, "K08356_gene_total_TPM_boxplot.pdf"), p_total, width = 6.4, height = 6.0)
svg(file.path(figdir, "K08356_gene_total_TPM_boxplot.svg"), width = 6.4, height = 6.0, family = "Times")
print(p_total)
dev.off()

kw_total_tbl <- tibble(
  feature = "K08356_total_gene_TPM",
  statistic = unname(kw_total$statistic),
  df = unname(kw_total$parameter),
  p_value = kw_total$p.value,
  n = nrow(total_gene),
  tested_value = "log10(TPM+1)",
  note = "Total Tree11 branch-supported K08356 gene TPM summed per sample"
)

write_csv(gene_long, file.path(outdir, "K08356_gene_branch_sample_long.csv"))
write_csv(total_gene, file.path(outdir, "K08356_gene_total_sample_long.csv"))
write_csv(summary_stats, file.path(outdir, "K08356_gene_branch_habitat_summary.csv"))
write_csv(kw_by_clade, file.path(outdir, "K08356_gene_branch_Kruskal.csv"))
write_csv(pairwise, file.path(outdir, "K08356_gene_branch_pairwise_Wilcoxon_BH.csv"))
write_csv(kw_total_tbl, file.path(outdir, "K08356_gene_total_Kruskal.csv"))

method_note <- tibble(
  item = c("abundance_used", "normalization", "source", "important_correction", "old_MAG_host_report_warning"),
  note = c(
    "sample_K08356_TPM from tree11_integrated_context_data.json; this is gene/contig-level K08356 TPM by sample and Tree11 branch",
    "TPM was already normalized upstream in the contig/KO workflow; this script uses per-sample TPM and log10(TPM+1), with no habitat-sum normalization",
    file.path(base, "tree11_integrated_context_data.json"),
    "This report does not use merged_tpm_matrix.txt host MAG TPM",
    "The previous MAG report measured abundance of K08356-bearing host MAGs, not K08356 gene abundance; it should not be used to claim K08356 enzyme abundance."
  )
)
write_csv(method_note, file.path(outdir, "README_method_note.csv"))

cat(outdir, "\n")

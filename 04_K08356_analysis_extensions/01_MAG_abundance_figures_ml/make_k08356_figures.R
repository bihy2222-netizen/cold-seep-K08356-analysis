suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(vegan)
})

outdir <- "derepMAG_abundance_stats"
figdir <- file.path(outdir, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

habitat_levels <- c("IS", "AS", "ES", "NS")
habitat_colors <- c(IS = "#4C78A8", AS = "#59A14F", ES = "#F28E2B", NS = "#E15759")

abund <- read.delim(file.path(outdir, "feature_abundance_by_sample.tsv"), check.names = FALSE)
prop <- read.delim(file.path(outdir, "feature_proportion_by_sample.tsv"), check.names = FALSE)
kw <- read.delim(file.path(outdir, "feature_kruskal.tsv"), check.names = FALSE)
kw_prop <- read.delim(file.path(outdir, "feature_proportion_kruskal.tsv"), check.names = FALSE)

abund$Habitat <- factor(abund$Habitat, levels = habitat_levels)
prop$Habitat <- factor(prop$Habitat, levels = habitat_levels)

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.3f", p)))
}

clean_feature <- function(x) {
  x |>
    sub("_TPM$", "", x = _) |>
    sub("_log10_TPM_plus_1$", "", x = _) |>
    sub("_proportion_of_unique_K08356_MAG_total$", "", x = _)
}

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#EEF2F7", color = "#D8DEE9"),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(color = "#444444")
  )

total_feature <- "K08356_bearing_unique_MAG_total"
total_p <- kw$Kruskal_p[kw$Feature == paste0(total_feature, "_TPM")]
total_df <- abund |> filter(Feature == total_feature)

p_total <- ggplot(total_df, aes(Habitat, log10_abundance_plus_1, fill = Habitat)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.72) +
  geom_jitter(aes(color = Habitat), width = 0.12, size = 1.9, alpha = 0.85, show.legend = FALSE) +
  scale_fill_manual(values = habitat_colors, drop = FALSE) +
  scale_color_manual(values = habitat_colors, drop = FALSE) +
  labs(
    title = "K08356-bearing derepMAG total abundance",
    subtitle = paste0("Kruskal-Wallis P = ", fmt_p(total_p)),
    x = "Habitat",
    y = "log10(TPM + 1)"
  ) +
  theme_pub
ggsave(file.path(figdir, "K08356_total_abundance_boxplot.png"), p_total, width = 6.4, height = 6.0, dpi = 320)

branch_df <- abund |>
  filter(Feature != total_feature) |>
  mutate(
    Feature_clean = clean_feature(Feature),
    Feature_clean = factor(Feature_clean, levels = c(
      "IdrA-associated",
      "canonical aioA-associated",
      "aioA-like-associated",
      "unknown AioA-like / uncertain DMSOR"
    ))
  )

branch_p <- kw |>
  filter(grepl("_TPM$", Feature), !grepl("unique_MAG_total", Feature)) |>
  transmute(Feature_clean = clean_feature(Feature), label = paste0("KW P = ", fmt_p(Kruskal_p)))

p_branch <- ggplot(branch_df, aes(Habitat, log10_abundance_plus_1, fill = Habitat)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.72) +
  geom_jitter(aes(color = Habitat), width = 0.12, size = 1.5, alpha = 0.82, show.legend = FALSE) +
  geom_text(
    data = branch_p,
    aes(x = 2.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.35,
    size = 3.25
  ) +
  facet_wrap(~ Feature_clean, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = habitat_colors, drop = FALSE) +
  scale_color_manual(values = habitat_colors, drop = FALSE) +
  labs(
    title = "K08356 branch-associated derepMAG abundance",
    subtitle = "Each facet shows global Kruskal-Wallis P across four habitats",
    x = "Habitat",
    y = "log10(TPM + 1)"
  ) +
  theme_pub
ggsave(file.path(figdir, "K08356_branch_abundance_boxplots.png"), p_branch, width = 8.0, height = 8.0, dpi = 320)

prop_df <- prop |>
  mutate(
    Feature_clean = factor(Feature, levels = c(
      "IdrA-associated",
      "canonical aioA-associated",
      "aioA-like-associated",
      "unknown AioA-like / uncertain DMSOR"
    ))
  )

prop_p <- kw_prop |>
  transmute(Feature_clean = clean_feature(Feature), label = paste0("KW P = ", fmt_p(Kruskal_p)))

p_prop <- ggplot(prop_df, aes(Habitat, Abundance, fill = Habitat)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.72) +
  geom_jitter(aes(color = Habitat), width = 0.12, size = 1.5, alpha = 0.82, show.legend = FALSE) +
  geom_text(
    data = prop_p,
    aes(x = 2.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.35,
    size = 3.25
  ) +
  facet_wrap(~ Feature_clean, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = habitat_colors, drop = FALSE) +
  scale_color_manual(values = habitat_colors, drop = FALSE) +
  labs(
    title = "K08356 branch proportion within unique K08356-bearing MAG total",
    subtitle = "Each facet shows global Kruskal-Wallis P across four habitats",
    x = "Habitat",
    y = "Proportion"
  ) +
  theme_pub
ggsave(file.path(figdir, "K08356_branch_proportion_boxplots.png"), p_prop, width = 8.0, height = 8.0, dpi = 320)

wide <- branch_df |>
  select(Sample, Habitat, Feature_clean, Abundance) |>
  pivot_wider(names_from = Feature_clean, values_from = Abundance, values_fill = 0)

comm <- wide |>
  select(-Sample, -Habitat) |>
  as.data.frame()
rownames(comm) <- wide$Sample
keep <- rowSums(comm) > 0
comm <- comm[keep, , drop = FALSE]
meta <- wide[keep, c("Sample", "Habitat")]
meta$Habitat <- factor(meta$Habitat, levels = habitat_levels)

set.seed(20260707)
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100, autotransform = FALSE, trace = FALSE)
anos <- anosim(comm, meta$Habitat, distance = "bray", permutations = 999)

scores_df <- as.data.frame(scores(nmds, display = "sites"))
scores_df$Sample <- rownames(scores_df)
scores_df <- left_join(scores_df, meta, by = "Sample")

nmds_label <- sprintf("ANOSIM R = %.3f, P = %.3g; NMDS stress = %.3f", anos$statistic, anos$signif, nmds$stress)

p_nmds <- ggplot(scores_df, aes(NMDS1, NMDS2, color = Habitat)) +
  stat_ellipse(aes(group = Habitat), type = "norm", linetype = 2, linewidth = 0.45, alpha = 0.65, show.legend = FALSE) +
  geom_point(size = 2.8, alpha = 0.9) +
  scale_color_manual(values = habitat_colors, drop = FALSE) +
  labs(
    title = "NMDS of K08356 branch-associated MAG abundance composition",
    subtitle = nmds_label,
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_pub
ggsave(file.path(figdir, "K08356_branch_NMDS_ANOSIM.png"), p_nmds, width = 6.4, height = 6.0, dpi = 320)

write.table(
  data.frame(
    Metric = c("ANOSIM_R", "ANOSIM_P", "NMDS_stress", "Samples_used", "Samples_removed_zero_total"),
    Value = c(anos$statistic, anos$signif, nmds$stress, nrow(comm), sum(!keep))
  ),
  file = file.path(outdir, "nmds_anosim_stats.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

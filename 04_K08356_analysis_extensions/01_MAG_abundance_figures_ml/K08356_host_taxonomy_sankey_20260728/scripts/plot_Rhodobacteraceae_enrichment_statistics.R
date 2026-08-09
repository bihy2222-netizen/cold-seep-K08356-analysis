suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
})

root_dir <- "/Users/catherine/Downloads/MAG-bin-tpm/K08356_host_taxonomy_sankey_20260728"
fig_dir <- file.path(root_dir, "figures")
res_dir <- file.path(root_dir, "results")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- "/Users/catherine/Downloads/营养盐热图绘制/results_revision/02_K08356_host_origin/K08356_sequence_MAG_taxonomy_master.tsv"

clade_cols <- c(
  "Clade 1" = "#80CDB7",
  "Clade 2" = "#BDBDBD",
  "Clade 3" = "#BFA7DF",
  "Clade 4" = "#8E68B1"
)

df <- read_tsv(input_path, show_col_types = FALSE) %>%
  mutate(
    Clade = factor(Clade, levels = c("Clade 1", "Clade 2", "Clade 3", "Clade 4")),
    Family_group = if_else(Family == "Rhodobacteraceae", "Rhodobacteraceae", "Other families")
  )

tab <- table(df$Clade, df$Family_group)
fisher_res <- fisher.test(tab)
chi_res <- suppressWarnings(chisq.test(tab, correct = FALSE))
cramers_v <- sqrt(unname(chi_res$statistic) /
                    (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)))

plot_counts <- as.data.frame(tab) %>%
  rename(Clade = Var1, Family_group = Var2, n_loci = Freq) %>%
  group_by(Clade) %>%
  mutate(total_loci = sum(n_loci),
         proportion = n_loci / total_loci,
         label = paste0(n_loci, "/", total_loci)) %>%
  ungroup()

rhodo_props <- plot_counts %>%
  filter(Family_group == "Rhodobacteraceae") %>%
  mutate(
    percent_label = percent(proportion, accuracy = 1),
    label_y = if_else(proportion == 0, 0.085, pmin(proportion + 0.055, 1.03))
  )

residual_df <- as.data.frame(as.table(chi_res$stdres)) %>%
  rename(Clade = Var1, Family_group = Var2, standardized_residual = Freq) %>%
  mutate(
    Clade = factor(Clade, levels = c("Clade 1", "Clade 2", "Clade 3", "Clade 4")),
    residual_label = sprintf("%.2f", standardized_residual)
  )

summary_tbl <- tibble(
  test = "Fisher exact test, Clade x Rhodobacteraceae-vs-other families",
  p_value = fisher_res$p.value,
  cramers_v = cramers_v,
  clade4_rhodobacteraceae_loci = tab["Clade 4", "Rhodobacteraceae"],
  clade4_total_loci = sum(tab["Clade 4", ]),
  clade4_rhodobacteraceae_percent = tab["Clade 4", "Rhodobacteraceae"] / sum(tab["Clade 4", ]),
  clade4_rhodobacteraceae_standardized_residual = chi_res$stdres["Clade 4", "Rhodobacteraceae"]
)

write_csv(plot_counts, file.path(res_dir, "Rhodobacteraceae_enrichment_plot_counts.csv"))
write_csv(residual_df, file.path(res_dir, "Rhodobacteraceae_enrichment_residual_heatmap_data.csv"))
write_csv(summary_tbl, file.path(res_dir, "Rhodobacteraceae_enrichment_plot_statistics.csv"))

stat_label <- paste0(
  "Fisher exact P = ", formatC(fisher_res$p.value, format = "e", digits = 2),
  "; Cramer's V = ", sprintf("%.3f", cramers_v),
  "; Clade 4 r = ", sprintf("%.2f", chi_res$stdres["Clade 4", "Rhodobacteraceae"])
)

p_bar <- ggplot(plot_counts, aes(x = Clade, y = proportion, fill = Family_group)) +
  geom_col(width = 0.68, color = "#30363D", linewidth = 0.25) +
  geom_text(
    data = rhodo_props,
    aes(x = Clade, y = label_y, label = paste0(label, "\n", percent_label)),
    inherit.aes = FALSE,
    family = "Times New Roman",
    fontface = "bold",
    size = 3.4,
    color = "#1F2933"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.12),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c("Rhodobacteraceae" = "#8E68B1", "Other families" = "#D8D8D8")) +
  labs(
    title = "A  Rhodobacteraceae-assigned loci by clade",
    subtitle = "Bars show locus-level proportions within each clade",
    x = NULL,
    y = "Proportion of K08356 loci",
    fill = NULL
  ) +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35)
  )

p_resid <- ggplot(residual_df, aes(x = Family_group, y = Clade, fill = standardized_residual)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = residual_label), family = "Times New Roman", fontface = "bold", size = 3.5) +
  scale_fill_gradient2(
    low = "#3B73B9", mid = "#F7F7F7", high = "#B64B6B",
    midpoint = 0, limits = c(-5, 5), oob = squish,
    name = "Std. residual"
  ) +
  labs(
    title = "B  Standardized residuals",
    subtitle = "Positive residuals indicate more loci than expected",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35)
  )

p_combined <- (p_bar | p_resid) +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    title = "Clade 4 enrichment of Rhodobacteraceae-assigned K08356 loci",
    subtitle = stat_label,
    caption = "Counts represent 48 K08356 candidate loci from 47 unique MAGs; this plot does not represent taxon abundance in the complete MAG background.",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(family = "Times New Roman", hjust = 0.5, size = 10.5),
      plot.caption = element_text(family = "Times New Roman", size = 8.5, color = "#4B5563", hjust = 0)
    )
  )

save_plot <- function(plot, filename, width = 10.8, height = 5.8) {
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  grDevices::svg(file.path(fig_dir, paste0(filename, ".svg")), width = width, height = height, family = "Times New Roman")
  print(plot)
  grDevices::dev.off()
  ggsave(file.path(fig_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 600)
  ggsave(file.path(fig_dir, paste0(filename, ".eps")), plot, width = width, height = height, device = cairo_ps)
}

save_plot(p_combined, "Rhodobacteraceae_enrichment_by_clade_statistics")

cat("Rhodobacteraceae enrichment figure done\n")
cat("Fisher exact P:", fisher_res$p.value, "\n")
cat("Cramer's V:", cramers_v, "\n")
cat("Clade 4 Rhodobacteraceae:", tab["Clade 4", "Rhodobacteraceae"], "/", sum(tab["Clade 4", ]), "\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggalluvial)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

root_dir <- "/Users/catherine/Downloads/MAG-bin-tpm/K08356_host_taxonomy_sankey_20260728"
input_path <- "/Users/catherine/Downloads/营养盐热图绘制/results_revision/02_K08356_host_origin/K08356_sequence_MAG_taxonomy_master.tsv"
fig_dir <- file.path(root_dir, "figures")
res_dir <- file.path(root_dir, "results")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

levels_tax <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
levels_clade_tax <- c("Clade", levels_tax)
prefix_map <- c(Phylum = "p__", Class = "c__", Order = "o__", Family = "f__", Genus = "g__", Species = "s__")

clean_taxon <- function(x, level) {
  prefix <- prefix_map[[level]]
  x <- ifelse(is.na(x), "", trimws(as.character(x)))
  x <- str_remove(x, "^[a-z]__")
  x <- ifelse(x == "" | x == "NA" | x == "Unclassified genus", "unclassified", x)
  paste0(prefix, x)
}

parse_species <- function(classification) {
  classification <- ifelse(is.na(classification), "", classification)
  hit <- str_match(classification, "(^|;)s__([^;]*)")[, 3]
  clean_taxon(hit, "Species")
}

display_taxon <- function(x) {
  x <- as.character(x)
  x <- str_remove(x, "^[a-z]__")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace(x, "^unclassified$", "Unclassified")
  x <- str_replace(x, "^rare species \\(<2 seq\\)$", "Rare species")
  x
}

truncate_label <- function(x, max_chars = 26) {
  if_else(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

df <- read_tsv(input_path, show_col_types = FALSE) %>%
  mutate(
    Clade = factor(Clade, levels = c("Clade 1", "Clade 2", "Clade 3", "Clade 4")),
    Clade_name_full = Clade_name,
    Phylum = clean_taxon(Phylum, "Phylum"),
    Class = clean_taxon(Class, "Class"),
    Order = clean_taxon(Order, "Order"),
    Family = clean_taxon(Family, "Family"),
    Genus = clean_taxon(Genus, "Genus"),
    Species = parse_species(Classification),
    Clade_name_short = case_when(
      str_detect(Clade_name, regex("canonical", ignore_case = TRUE)) ~ "canonical aioA-associated",
      str_detect(Clade_name, regex("unknown|uncertain", ignore_case = TRUE)) ~ "unknown AioA-like / uncertain DMSOR",
      str_detect(Clade_name, regex("IdrA", ignore_case = TRUE)) ~ "IdrA-associated",
      str_detect(Clade_name, regex("aioA-like", ignore_case = TRUE)) ~ "aioA-like-associated",
      TRUE ~ Clade_name
    )
  )

path_counts_full <- df %>%
  group_by(across(all_of(levels_tax))) %>%
  summarise(
    n_sequences = n(),
    n_unique_MAGs = n_distinct(MAG_ID),
    clades = paste(sort(unique(Clade_name_short)), collapse = "; "),
    habitats = paste(sort(unique(Source_habitat)), collapse = "; "),
    MAG_IDs = paste(sort(unique(MAG_ID)), collapse = "; "),
    sequence_IDs = paste(K08356_sequence_ID, collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sequences), Phylum, Class, Order, Family, Genus, Species)

species_counts <- path_counts_full %>%
  group_by(Species) %>%
  summarise(total_sequences = sum(n_sequences), .groups = "drop")

rare_species <- species_counts %>%
  filter(total_sequences < 2, Species != "s__unclassified") %>%
  pull(Species)

path_counts_plot <- path_counts_full %>%
  mutate(Species_plot = if_else(Species %in% rare_species, "s__rare species (<2 seq)", Species)) %>%
  group_by(Phylum, Class, Order, Family, Genus, Species = Species_plot) %>%
  summarise(
    n_sequences = sum(n_sequences),
    n_unique_MAGs = sum(n_unique_MAGs),
    clades = paste(sort(unique(unlist(str_split(clades, "; ")))), collapse = "; "),
    habitats = paste(sort(unique(unlist(str_split(habitats, "; ")))), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sequences), Phylum, Class, Order, Family, Genus, Species)

path_counts_clade_full <- df %>%
  group_by(across(all_of(levels_clade_tax))) %>%
  summarise(
    n_sequences = n(),
    n_unique_MAGs = n_distinct(MAG_ID),
    clade_name = paste(sort(unique(Clade_name_full)), collapse = "; "),
    habitats = paste(sort(unique(Source_habitat)), collapse = "; "),
    MAG_IDs = paste(sort(unique(MAG_ID)), collapse = "; "),
    sequence_IDs = paste(K08356_sequence_ID, collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Clade, desc(n_sequences), Phylum, Class, Order, Family, Genus, Species)

species_counts_by_clade <- path_counts_clade_full %>%
  group_by(Clade, Species) %>%
  summarise(total_sequences = sum(n_sequences), .groups = "drop")

rare_species_by_clade <- species_counts_by_clade %>%
  filter(total_sequences < 2, Species != "s__unclassified") %>%
  transmute(Clade, Species, rare_species = TRUE)

path_counts_clade_plot <- path_counts_clade_full %>%
  left_join(rare_species_by_clade, by = c("Clade", "Species")) %>%
  mutate(Species_plot = if_else(!is.na(rare_species), "s__rare species (<2 seq)", Species)) %>%
  select(-rare_species) %>%
  group_by(Clade, Phylum, Class, Order, Family, Genus, Species = Species_plot) %>%
  summarise(
    n_sequences = sum(n_sequences),
    n_unique_MAGs = sum(n_unique_MAGs),
    clade_name = paste(sort(unique(unlist(str_split(clade_name, "; ")))), collapse = "; "),
    habitats = paste(sort(unique(unlist(str_split(habitats, "; ")))), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Clade, desc(n_sequences), Phylum, Class, Order, Family, Genus, Species)

write_csv(df, file.path(res_dir, "K08356_host_taxonomy_source_table_phylum_to_species.csv"))
write_csv(path_counts_full, file.path(res_dir, "K08356_host_taxonomy_path_counts_full_species.csv"))
write_csv(path_counts_plot, file.path(res_dir, "K08356_host_taxonomy_path_counts_plot_species_rare_collapsed.csv"))
write_csv(path_counts_clade_full, file.path(res_dir, "K08356_host_taxonomy_clade_path_counts_full_species.csv"))
write_csv(path_counts_clade_plot, file.path(res_dir, "K08356_host_taxonomy_clade_path_counts_plot_species_rare_collapsed.csv"))

family_tab <- table(df$Clade, df$Family)
set.seed(20260730)
family_fisher_sim <- fisher.test(family_tab, simulate.p.value = TRUE, B = 100000)
family_chisq <- suppressWarnings(chisq.test(family_tab, correct = FALSE))
family_cramers_v <- sqrt(unname(family_chisq$statistic) /
                           (sum(family_tab) * min(nrow(family_tab) - 1, ncol(family_tab) - 1)))
family_residuals <- as.data.frame(as.table(family_chisq$stdres)) %>%
  rename(Clade = Var1, Family = Var2, standardized_residual = Freq) %>%
  arrange(desc(abs(standardized_residual)))

rhodo_family <- "f__Rhodobacteraceae"
rhodo_tab <- table(df$Clade, if_else(df$Family == rhodo_family, "Rhodobacteraceae", "Other families"))
rhodo_fisher <- fisher.test(rhodo_tab)
rhodo_chisq <- suppressWarnings(chisq.test(rhodo_tab, correct = FALSE))
rhodo_cramers_v <- sqrt(unname(rhodo_chisq$statistic) /
                          (sum(rhodo_tab) * min(nrow(rhodo_tab) - 1, ncol(rhodo_tab) - 1)))
rhodo_residuals <- as.data.frame(as.table(rhodo_chisq$stdres)) %>%
  rename(Clade = Var1, Family_group = Var2, standardized_residual = Freq) %>%
  arrange(desc(abs(standardized_residual)))

family_stats_summary <- tibble::tibble(
  test = c("Clade x Family Fisher-Freeman-Halton Monte Carlo",
           "Clade x Rhodobacteraceae-vs-other Fisher exact"),
  table_scope = c("all GTDB family categories", "Rhodobacteraceae vs all other families"),
  p_value = c(family_fisher_sim$p.value, rhodo_fisher$p.value),
  p_value_note = c("Monte Carlo simulation, B=100000; exact R calculation overflows for sparse 4 x family table",
                   "Exact Fisher test"),
  cramers_v = c(family_cramers_v, rhodo_cramers_v),
  key_result = c(
    paste0("Clade 4 x Rhodobacteraceae standardized residual = ",
           round(family_chisq$stdres["Clade 4", rhodo_family], 3)),
    paste0("Clade 4: ", rhodo_tab["Clade 4", "Rhodobacteraceae"], "/",
           sum(rhodo_tab["Clade 4", ]), " loci in Rhodobacteraceae; standardized residual = ",
           round(rhodo_chisq$stdres["Clade 4", "Rhodobacteraceae"], 3))
  )
)
write_csv(as.data.frame.matrix(family_tab) %>% tibble::rownames_to_column("Clade"),
          file.path(res_dir, "Clade_by_family_contingency_table.csv"))
write_csv(as.data.frame.matrix(rhodo_tab) %>% tibble::rownames_to_column("Clade"),
          file.path(res_dir, "Clade_by_Rhodobacteraceae_vs_other_table.csv"))
write_csv(family_stats_summary, file.path(res_dir, "Clade_family_enrichment_statistics_summary.csv"))
write_csv(family_residuals, file.path(res_dir, "Clade_family_standardized_residuals.csv"))
write_csv(rhodo_residuals, file.path(res_dir, "Clade_Rhodobacteraceae_standardized_residuals.csv"))

palette <- c(
  "p__Pseudomonadota" = "#1F78B4",
  "p__Actinomycetota" = "#33A02C",
  "p__Chloroflexota" = "#FF7F00",
  "p__Thermoproteota" = "#6A3D9A",
  "p__Nitrospinota" = "#B15928",
  "p__Campylobacterota" = "#E31A1C",
  "p__Desulfobacterota" = "#A6CEE3",
  "p__unclassified" = "#999999"
)
missing_phyla <- setdiff(unique(path_counts_plot$Phylum), names(palette))
if (length(missing_phyla) > 0) {
  extra <- grDevices::hcl.colors(length(missing_phyla), palette = "Dark 3")
  names(extra) <- missing_phyla
  palette <- c(palette, extra)
}

plot_df <- path_counts_plot %>%
  mutate(path_id = row_number(), Phylum_fill = Phylum) %>%
  pivot_longer(cols = all_of(levels_tax), names_to = "Rank", values_to = "Taxon") %>%
  mutate(
    Rank = factor(Rank, levels = levels_tax),
    Taxon_label = Taxon,
    Taxon_label = if_else(nchar(Taxon_label) > 34, paste0(substr(Taxon_label, 1, 31), "..."), Taxon_label)
  )

axis_labels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
names(axis_labels) <- levels_tax

p <- ggplot(
  plot_df,
  aes(x = Rank, stratum = Taxon_label, alluvium = path_id, y = n_sequences, fill = Phylum_fill)
) +
  geom_flow(alpha = 0.42, color = "grey72", linewidth = 0.18, curve_type = "quintic") +
  geom_stratum(width = 0.24, color = "grey30", linewidth = 0.25, alpha = 0.96) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.7,
    family = "Times New Roman",
    color = "black",
    min.y = 0.8
  ) +
  scale_x_discrete(labels = axis_labels, expand = c(0.045, 0.045)) +
  scale_fill_manual(values = palette, drop = FALSE, na.translate = FALSE) +
  labs(
    title = "K08356-bearing MAG host taxonomy Sankey",
    subtitle = "Taxonomic flow from phylum to species; weights are K08356 sequence counts, rare species collapsed for readability",
    x = NULL,
    y = "Number of K08356 sequences",
    fill = "Phylum"
  ) +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
    plot.margin = margin(12, 18, 12, 18)
  )

save_plot <- function(plot, filename, width = 13.5, height = 8.0) {
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  grDevices::svg(file.path(fig_dir, paste0(filename, ".svg")), width = width, height = height, family = "Times New Roman")
  print(plot)
  grDevices::dev.off()
  ggsave(file.path(fig_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 600)
  ggsave(file.path(fig_dir, paste0(filename, ".eps")), plot, width = width, height = height, device = cairo_ps)
}

save_plot(p, "K08356_host_taxonomy_sankey_phylum_to_species")

clade_cols <- c(
  "Clade 1" = "#80CDB7",
  "Clade 2" = "#BDBDBD",
  "Clade 3" = "#BFA7DF",
  "Clade 4" = "#8E68B1"
)

clade_labels <- df %>%
  distinct(Clade, Clade_name_full) %>%
  left_join(df %>% count(Clade, name = "n"), by = "Clade") %>%
  arrange(Clade) %>%
  mutate(label = paste0(as.character(Clade), " (n=", n, "): ", Clade_name_full))

plot_clade_df <- path_counts_clade_plot %>%
  mutate(path_id = row_number(), Clade_fill = as.character(Clade)) %>%
  pivot_longer(cols = all_of(levels_clade_tax), names_to = "Rank", values_to = "Taxon") %>%
  mutate(
    Rank = factor(Rank, levels = levels_clade_tax),
    Taxon_label = if_else(Rank == "Clade", as.character(Taxon), display_taxon(Taxon)),
    Taxon_label = truncate_label(Taxon_label, 28)
  )

node_labels_clade <- plot_clade_df %>%
  group_by(Rank, Taxon_label) %>%
  summarise(node_n = sum(n_sequences), .groups = "drop") %>%
  mutate(Node_label_with_n = paste0(Taxon_label, "\n(n=", node_n, ")"))

plot_clade_df <- plot_clade_df %>%
  left_join(node_labels_clade, by = c("Rank", "Taxon_label"))

axis_labels_clade <- c("Clade", "Phylum", "Class", "Order", "Family", "Genus", "Species")
names(axis_labels_clade) <- levels_clade_tax

p_clade <- ggplot(
  plot_clade_df,
  aes(x = Rank, stratum = Taxon_label, alluvium = path_id, y = n_sequences)
) +
  geom_flow(aes(fill = Clade_fill), alpha = 0.56, color = "white", linewidth = 0.16,
            curve_type = "quintic") +
  geom_stratum(width = 0.28, fill = "#F8FAFC", color = "#5B6773", linewidth = 0.38) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.45,
    family = "Times New Roman",
    fontface = "bold",
    color = "#1F2933",
    min.y = 0.85
  ) +
  scale_x_discrete(labels = axis_labels_clade, expand = c(0.035, 0.035)) +
  scale_fill_manual(
    values = clade_cols,
    breaks = clade_labels$Clade,
    labels = clade_labels$label,
    name = NULL,
    drop = FALSE
  ) +
  labs(
    title = "Host taxonomic distribution of K08356 loci across four phylogenetic clades",
    subtitle = "Taxonomic flows from clade to species-level assignment; ribbon widths represent K08356 locus counts (48 loci from 47 MAGs)",
    x = NULL,
    y = "Number of K08356 sequences"
  ) +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 10.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.58, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
    plot.margin = margin(10, 16, 10, 16)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

save_plot(p_clade, "K08356_host_taxonomy_sankey_clade_to_species", width = 16.5, height = 8.8)

p_clade_n <- ggplot(
  plot_clade_df,
  aes(x = Rank, stratum = Node_label_with_n, alluvium = path_id, y = n_sequences)
) +
  geom_flow(aes(fill = Clade_fill), alpha = 0.58, color = "white", linewidth = 0.17,
            curve_type = "quintic") +
  geom_stratum(width = 0.30, fill = "#F8FAFC", color = "#4B5563", linewidth = 0.42) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 2.35,
    family = "Times New Roman",
    fontface = "bold",
    color = "#111827",
    lineheight = 0.82,
    min.y = 1.15
  ) +
  scale_x_discrete(labels = axis_labels_clade, expand = c(0.03, 0.03)) +
  scale_fill_manual(
    values = clade_cols,
    breaks = clade_labels$Clade,
    labels = clade_labels$label,
    name = NULL,
    drop = FALSE
  ) +
  labs(
    title = "Host taxonomic distribution of K08356 loci across four phylogenetic clades",
    subtitle = "Taxonomic flows from clade to species-level assignment; ribbon widths represent K08356 locus counts (48 loci from 47 MAGs)",
    x = NULL,
    y = "Number of K08356 sequences"
  ) +
  theme_classic(base_family = "Times New Roman", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 10.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 8.8),
    legend.key.width = unit(0.55, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.38),
    plot.margin = margin(10, 16, 10, 16)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

save_plot(p_clade_n, "K08356_host_taxonomy_sankey_clade_to_species_with_n", width = 17.5, height = 9.4)

readme <- c(
  "# K08356 host taxonomy Sankey",
  "",
  paste0("Input: ", input_path),
  "",
  "Taxonomic ranks shown: Phylum -> Class -> Order -> Family -> Genus -> Species.",
  "Weights are K08356 sequence counts from the 48-row K08356 MAG host taxonomy master table.",
  "Species was parsed from the GTDB-style Classification field. Empty species annotations were written as s__unclassified.",
  "For readability, species observed in fewer than 2 sequences were collapsed to s__rare species (<2 seq) in the plotted figure.",
  "The full uncollapsed path count table is retained in results/K08356_host_taxonomy_path_counts_full_species.csv.",
  "",
  "Statistical note:",
  paste0("- Clade x Rhodobacteraceae-vs-other Fisher exact test: P = ",
         signif(rhodo_fisher$p.value, 4), "; Cramer's V = ", round(rhodo_cramers_v, 3), "."),
  paste0("- Clade 4 contained ", rhodo_tab["Clade 4", "Rhodobacteraceae"], " of ",
         sum(rhodo_tab["Clade 4", ]), " K08356 loci in Rhodobacteraceae-assigned MAGs; standardized residual r = ",
         round(rhodo_chisq$stdres["Clade 4", "Rhodobacteraceae"], 2), "."),
  "- These statistics describe enrichment within the 48 candidate-locus dataset and do not represent taxon abundance in the complete 808-MAG background.",
  "",
  "Outputs:",
  "- figures/K08356_host_taxonomy_sankey_phylum_to_species.pdf",
  "- figures/K08356_host_taxonomy_sankey_phylum_to_species.svg",
  "- figures/K08356_host_taxonomy_sankey_phylum_to_species.png",
  "- figures/K08356_host_taxonomy_sankey_phylum_to_species.eps",
  "- figures/K08356_host_taxonomy_sankey_clade_to_species.pdf",
  "- figures/K08356_host_taxonomy_sankey_clade_to_species.svg",
  "- figures/K08356_host_taxonomy_sankey_clade_to_species.png",
  "- figures/K08356_host_taxonomy_sankey_clade_to_species.eps",
  "- figures/K08356_host_taxonomy_sankey_clade_to_species_with_n.pdf",
  "- figures/K08356_host_taxonomy_sankey_clade_to_species_with_n.svg",
  "- figures/K08356_host_taxonomy_sankey_clade_to_species_with_n.png",
  "- figures/K08356_host_taxonomy_sankey_clade_to_species_with_n.eps"
)
writeLines(readme, file.path(root_dir, "README.md"))

cat("Sankey done\n")
cat("Input rows:", nrow(df), "\n")
cat("Full paths:", nrow(path_counts_full), "\n")
cat("Plot paths:", nrow(path_counts_plot), "\n")

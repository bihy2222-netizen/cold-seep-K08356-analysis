suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readxl)
  library(stringr)
})

out_dir <- "/Users/catherine/Downloads/MAG-bin-tpm/K08356_MAG_gene_island_arrows_20260708"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

branch_map_file <- "/Users/catherine/Downloads/MAG-bin-tpm/derepMAG_abundance_stats/k08356_branch_feature_map.tsv"
hitdata_file <- "/Users/catherine/Downloads/hitdata.txt"
sy_excel <- "/Users/catherine/Downloads/箭头图/箭头图原始文件.xlsx"

branch_order <- c(
  "IdrA-associated",
  "canonical aioA-associated",
  "aioA-like-associated",
  "unknown AioA-like / uncertain DMSOR"
)

branch_colors <- c(
  "IdrA-associated" = "#5B8DB8",
  "canonical aioA-associated" = "#2CA25F",
  "aioA-like-associated" = "#009E73",
  "unknown AioA-like / uncertain DMSOR" = "#CC79A7"
)

gene_colors <- c(
  "Target K08356" = "#65A88F",
  "Rieske/AioB" = "#B9A7D8",
  "Arsenic resistance" = "#A987B5",
  "Efflux/transporter" = "#B9B9B9",
  "Nitrogen/reduction" = "#8FB7D9",
  "Electron transfer" = "#E7B85A",
  "Other" = "#BDBDBD"
)

class_gene <- function(gene) {
  case_when(
    gene %in% c("AioA", "IdrA", "AioA-like") ~ "Target K08356",
    gene %in% c("AioB", "Rieske", "Rieske subunit") ~ "Rieske/AioB",
    grepl("^Ars|ACR3", gene) ~ "Arsenic resistance",
    grepl("Pst|Tar|transporter|Virulence", gene, ignore.case = TRUE) ~ "Efflux/transporter",
    gene %in% c("NirD", "NirB", "SseA") ~ "Nitrogen/reduction",
    gene %in% c("Fer2", "GlpF", "YjbI") ~ "Electron transfer",
    TRUE ~ "Other"
  )
}

arrow_poly <- function(df, head_frac = 0.18) {
  pieces <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    xmin <- r$Start
    xmax <- r$End
    y <- r$Molecule_y
    h <- 0.32
    len <- max(xmax - xmin, 1)
    head <- min(len * head_frac, len * 0.45, 450)
    if (r$Direction >= 0) {
      xs <- c(xmin, xmax - head, xmax - head, xmax, xmax - head, xmax - head, xmin)
    } else {
      xs <- c(xmax, xmin + head, xmin + head, xmin, xmin + head, xmin + head, xmax)
    }
    ys <- c(y - h, y - h, y - h * 1.6, y, y + h * 1.6, y + h, y + h)
    pieces[[i]] <- data.frame(
      row_id = r$row_id,
      Molecule = r$Molecule,
      Gene = r$Gene,
      Gene_class = r$Gene_class,
      x = xs,
      y = ys
    )
  }
  bind_rows(pieces)
}

parse_pos <- function(x) {
  nums <- str_extract_all(gsub(",", "", x), "\\d+")[[1]]
  as.numeric(nums[1:2])
}

## Three available canonical aioA-associated gene neighborhoods.
sy_raw <- read_excel(sy_excel) %>%
  rename(Position = `位置 (bp)`, Gene = `基因`, Direction = `链`, Module = `功能模块`)
sy_genes <- sy_raw %>%
  rowwise() %>%
  mutate(Start = parse_pos(Position)[1], End = parse_pos(Position)[2]) %>%
  ungroup() %>%
  mutate(Molecule = "SY366YW-4-8_bin7-k141_14765\n(canonical aioA-associated; IS)")

s1_genes <- tibble(
  Molecule = "S1_9-12_bin1-k141_637215\n(canonical aioA-associated; AS)",
  Gene = c("Rieske subunit", "AioA", "SseA"),
  Start = c(1070, 1613, 4246),
  End = c(1588, 4084, 5181),
  Direction = c(1, 1, 1),
  Module = c("Rieske/AioB", "Target", "Reduction")
)

r2111_genes <- tibble(
  Molecule = "R2111_S300_0-10_bin19-k141_825759\n(canonical aioA-associated; NS)",
  Gene = c("Acs", "CaiD", "GlpF", "YjbI", "AioA", "Rieske", "NirD", "Fer2", "NirB",
           "HTH_metalloreg", "PRK07534", "Virulence_fact"),
  Start = c(2, 1181, 2845, 3805, 5705, 8198, 8821, 9301, 9624, 12597, 13612, 14795),
  End = c(730, 2848, 3687, 4716, 8194, 8704, 9150, 9621, 10835, 13577, 14604, 15088),
  Direction = c(1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1),
  Module = c("Other", "Other", "Electron transfer", "Electron transfer", "Target",
             "Rieske/AioB", "Reduction", "Electron transfer", "Reduction",
             "Regulation", "Other", "Transporter")
)

gene_islands <- bind_rows(sy_genes, s1_genes, r2111_genes) %>%
  mutate(
    Start = as.numeric(Start),
    End = as.numeric(End),
    Gene_class = class_gene(Gene),
    Molecule_f = factor(Molecule, levels = unique(Molecule)),
    Molecule_y = as.numeric(Molecule_f) * 2,
    row_id = row_number(),
    Mid = (Start + End) / 2
  )

write.csv(gene_islands, file.path(out_dir, "available_canonical_gene_island_coordinates.csv"), row.names = FALSE)

poly_df <- arrow_poly(gene_islands)

island_plot <- ggplot() +
  geom_segment(
    data = gene_islands %>% group_by(Molecule, Molecule_f, Molecule_y) %>% summarise(xmin = min(Start), xmax = max(End), .groups = "drop"),
    aes(x = xmin, xend = xmax, y = Molecule_y, yend = Molecule_y),
    color = "grey70",
    linewidth = 0.35
  ) +
  geom_polygon(
    data = poly_df,
    aes(x = x, y = y, group = row_id, fill = Gene_class),
    color = "grey20",
    linewidth = 0.25
  ) +
  geom_text(
    data = gene_islands,
    aes(x = Mid, y = Molecule_y + 0.78, label = Gene),
    size = 2.85,
    family = "Times"
  ) +
  scale_fill_manual(values = gene_colors, drop = FALSE) +
  scale_y_continuous(
    breaks = seq_along(levels(gene_islands$Molecule_f)) * 2,
    labels = levels(gene_islands$Molecule_f),
    expand = expansion(add = 1.2)
  ) +
  scale_x_continuous(labels = function(x) paste0(x / 1000, "k")) +
  labs(
    title = "Available MAG gene neighborhoods around canonical aioA-associated K08356",
    subtitle = "Three neighborhoods can be drawn from available CDS coordinate files; colors follow Fig. 5-style functional highlighting",
    x = "Genomic position (bp)",
    y = NULL,
    fill = "Gene/function class"
  ) +
  theme_classic(base_family = "Times", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed")
  )

ggsave(file.path(out_dir, "available_canonical_MAG_gene_island_arrows.png"), island_plot, width = 12, height = 6.6, dpi = 450, bg = "white")
ggsave(file.path(out_dir, "available_canonical_MAG_gene_island_arrows.pdf"), island_plot, width = 12, height = 6.6, bg = "white")

## Four-branch K08356 protein domain architecture from CD-search hitdata.
branch_map <- read.delim(branch_map_file, check.names = FALSE) %>%
  mutate(Feature = factor(Feature, levels = branch_order))

hit_lines <- readLines(hitdata_file)
hit_table <- read.delim(text = paste(hit_lines[!grepl("^#", hit_lines) & nzchar(hit_lines)], collapse = "\n"),
                        check.names = FALSE, stringsAsFactors = FALSE)

domains <- hit_table %>%
  mutate(
    Original_MAG_ID = sub("^Q#[0-9]+ - >", "", Query),
    Domain = .data[["Short name"]],
    From = as.numeric(From),
    To = as.numeric(To)
  ) %>%
  left_join(branch_map %>% select(MAG, Feature, Original_MAG_ID, Corrected_habitat), by = "Original_MAG_ID") %>%
  filter(!is.na(Feature)) %>%
  mutate(
    Protein_label = paste0(Original_MAG_ID, " | ", Corrected_habitat),
    Domain_class = case_when(
      grepl("arsenite_ox_L", Domain) ~ "Molybdopterin oxidoreductase",
      grepl("MopB_CT", Domain) ~ "MopB_CT",
      grepl("Molybdopterin-Binding", Domain) ~ "Molybdopterin-binding",
      TRUE ~ Domain
    )
  )

protein_order <- domains %>%
  group_by(Feature, Protein_label) %>%
  summarise(start = min(From), end = max(To), .groups = "drop") %>%
  arrange(Feature, Protein_label) %>%
  pull(Protein_label)

domains <- domains %>%
  mutate(Protein_f = factor(Protein_label, levels = rev(unique(protein_order))))

write.csv(domains, file.path(out_dir, "K08356_four_branch_CDsearch_domain_hits.csv"), row.names = FALSE)

domain_colors <- c(
  "Molybdopterin oxidoreductase" = "#7B5AA6",
  "MopB_CT" = "#D95F02",
  "Molybdopterin-binding" = "#1B9E77"
)

domain_plot <- ggplot(domains, aes(y = Protein_f)) +
  geom_segment(
    data = domains %>% group_by(Feature, Protein_f) %>% summarise(xmin = min(From), xmax = max(To), .groups = "drop"),
    aes(x = xmin, xend = xmax, y = Protein_f, yend = Protein_f),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.35
  ) +
  geom_segment(
    aes(x = From, xend = To, yend = Protein_f, color = Domain_class),
    linewidth = 4.2,
    lineend = "butt"
  ) +
  facet_grid(Feature ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = domain_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.04))) +
  labs(
    title = "K08356 protein domain architecture across four MAG-associated branches",
    subtitle = "CD-search hits from hitdata.txt; this is protein-domain evidence, not a gene-neighborhood plot",
    x = "Amino-acid position",
    y = NULL,
    color = "CD-search domain"
  ) +
  theme_classic(base_family = "Times", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 5.8),
    axis.title.x = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.spacing.y = unit(0.35, "lines")
  )

ggsave(file.path(out_dir, "K08356_four_branch_protein_domain_architecture.png"), domain_plot, width = 10, height = 13.5, dpi = 450, bg = "white")
ggsave(file.path(out_dir, "K08356_four_branch_protein_domain_architecture.pdf"), domain_plot, width = 10, height = 13.5, bg = "white")

readme <- c(
  "# K08356 MAG gene-island arrow figure notes",
  "",
  "Generated outputs:",
  "- `available_canonical_MAG_gene_island_arrows.png/pdf`: true gene-neighborhood arrow plot for the three MAGs with available CDS coordinate information.",
  "- `K08356_four_branch_protein_domain_architecture.png/pdf`: CD-search protein-domain architecture for all 48 K08356 proteins across the four branch groups.",
  "",
  "Important limitation:",
  "Only three `.cds.faa` files are currently available on disk: `S1_9-12_bin1.cds.faa`, `SY366YW-4-8_bin7.cds.faa`, and `R2111_S300_0-10_bin19.cds.faa`. These cover the canonical aioA-associated examples drawn previously.",
  "The supplied `/Users/catherine/Downloads/hitdata.txt` contains CD-search domain hits for the target K08356 proteins, but it does not contain upstream/downstream ORF coordinates or annotations. Therefore it can support the four-branch protein-domain architecture plot, but cannot by itself support a true four-branch gene-island plot.",
  "",
  "To draw true gene-island arrows for IdrA-associated, aioA-like-associated, and unknown AioA-like / uncertain DMSOR groups, provide the corresponding MAG `.cds.faa`, `.gff`, `.gbk`, or Prodigal output files for those MAGs/contigs.",
  "",
  "Reference style used:",
  "The gene-neighborhood panel follows the logic of IdrA ISMEJ 2021 Fig. 5: target molybdopterin oxidoreductase neighborhood, functional coloring, and adjacent gene arrows with strand direction."
)
writeLines(readme, file.path(out_dir, "README_K08356_MAG_gene_island_arrows.md"))

message("Wrote K08356 MAG arrow/domain figures to: ", out_dir)

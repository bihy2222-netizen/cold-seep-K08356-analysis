suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(vegan)
  library(cluster)
})

root_dir <- "/Users/catherine/Downloads/MAG-bin-tpm/机器学习/03_protein_sequence_space_validation_v1"
fasta_path <- "/Users/catherine/Downloads/aioA 蛋白序列建树/tree11_no_uncertain_DryadAioAIdrAfinal/all_proteins.faa"
clade_detail_path <- "/Users/catherine/Downloads/aioA 蛋白序列建树/tree11_no_uncertain_DryadAioAIdrAfinal/tree11_clade_detail.csv"
mag_label_path <- "/Users/catherine/Downloads/MAG-bin-tpm/机器学习/00_input/K08356_branch_label_table.csv"

dir.create(file.path(root_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "Rlib"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(root_dir, "Rlib"), .libPaths()))

set.seed(20260708)

read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header_idx <- grep("^>", lines)
  records <- vector("list", length(header_idx))
  for (i in seq_along(header_idx)) {
    start <- header_idx[i]
    end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
    header <- sub("^>", "", lines[start])
    seq <- paste(lines[(start + 1):end], collapse = "")
    seq <- toupper(gsub("[^A-Z]", "", seq))
    records[[i]] <- data.frame(protein_id = header, sequence = seq, stringsAsFactors = FALSE)
  }
  bind_rows(records)
}

infer_reference_branch <- function(id) {
  case_when(
    grepl("^canonical_AioA_|^Dryad_canonical_AioA|AioA_or_AioA-like", id) ~ "canonical aioA-associated",
    grepl("^AioA_like_", id) ~ "aioA-like-associated",
    grepl("^IdrA_|Dryad_IdrA|Fig5_IdrA", id) ~ "IdrA-associated",
    grepl("Fig5_unknown|Unknown", id) ~ "unknown AioA-like / uncertain DMSOR",
    grepl("^MicrobiologySpectrum_AioA\\|", id) ~ "canonical aioA-associated",
    grepl("^MicrobiologySpectrum_IdrA\\|", id) ~ "IdrA-associated",
    grepl("^ISMEJFig5\\|Arsenite_oxidase_clade\\|", id) ~ "canonical aioA-associated",
    grepl("^ISMEJFig5\\|Iodate_reductase_clade\\|", id) ~ "IdrA-associated",
    grepl("^ISMEJFig5\\|Unknown_clade\\|", id) ~ "unknown AioA-like / uncertain DMSOR",
    grepl("AioA_Like|AioA-like|aioA-like", id) ~ "aioA-like-associated",
    TRUE ~ "other DMSOR-family reference-associated"
  )
}

infer_reference_set <- function(id) {
  case_when(
    grepl("^ColdSeep_contig_K08356_", id) ~ "Cold seep contig target",
    grepl("^canonical_AioA_|^AioA_like_|^Dryad_canonical_AioA", id) ~ "AioA/AioA-like reference",
    grepl("^IdrA_|^Dryad_IdrA", id) ~ "IdrA reference",
    grepl("^Fig5_|ISMEJFig5", id) ~ "ISMEJ Fig5 reference",
    grepl("^NarH", id) ~ "NarH outgroup",
    grepl("^MicrobiologySpectrum_", id) ~ "MicrobiologySpectrum reference",
    grepl("^ISMEJFig5\\|", id) ~ "ISMEJ Fig5 reference",
    TRUE ~ "legacy DMSOR-family reference"
  )
}

aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
di <- as.vector(outer(aa, aa, paste0))

calc_features <- function(seq) {
  chars <- strsplit(seq, "")[[1]]
  valid <- chars[chars %in% aa]
  len <- length(valid)
  aa_counts <- table(factor(valid, levels = aa))
  aa_freq <- as.numeric(aa_counts) / max(len, 1)
  names(aa_freq) <- paste0("AA_", aa)

  if (len >= 2) {
    pairs <- paste0(valid[-len], valid[-1])
    di_counts <- table(factor(pairs, levels = di))
    di_freq <- as.numeric(di_counts) / (len - 1)
  } else {
    di_freq <- rep(0, length(di))
  }
  names(di_freq) <- paste0("DI_", di)

  c(length_aa = len, aa_freq, di_freq)
}

clean_branch <- function(x) {
  case_when(
    x == "IdrA-associated clade" ~ "IdrA-associated",
    x == "Fig5 unknown AioA-like / uncertain DMSOR clade" ~ "unknown AioA-like / uncertain DMSOR",
    x == "canonical aioA" ~ "canonical aioA-associated",
    x == "unknown_uncertain" ~ "unknown AioA-like / uncertain DMSOR",
    TRUE ~ x
  )
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

normalize_coldseep_id <- function(id) {
  id %>%
    sub("^ColdSeep_contig_K08356_", "", .) %>%
    sub("_unknown$", "", .)
}

clade_tbl <- read_csv(clade_detail_path, show_col_types = FALSE) %>%
  mutate(branch = clean_branch(clade),
         normalized_id = sequence_id,
         source_type = "contig-only",
         is_reference = FALSE) %>%
  select(normalized_id, branch, habitat, sample_group, display, source_type, is_reference)

mag_label_tbl <- read_csv(mag_label_path, show_col_types = FALSE) %>%
  mutate(branch_mag = clean_branch(branch),
         mag_contig_normalized_id = paste0(sample_id, "-", target_orf_id)) %>%
  select(mag_contig_normalized_id, MAG_id, target_gene_id, branch_mag)

seq_tbl <- read_fasta(fasta_path) %>%
  mutate(length_aa = nchar(sequence),
         normalized_id = normalize_coldseep_id(protein_id),
         is_coldseep_target = grepl("^ColdSeep_contig_K08356_", protein_id))

metadata <- seq_tbl %>%
  left_join(clade_tbl, by = "normalized_id") %>%
  left_join(mag_label_tbl, by = c("normalized_id" = "mag_contig_normalized_id")) %>%
  mutate(
    branch = if_else(is.na(branch), infer_reference_branch(protein_id), branch),
    branch = clean_branch(branch),
    source_type = if_else(is.na(source_type), infer_reference_set(protein_id), source_type),
    source_type = if_else(!is.na(target_gene_id), "MAG-derived contig", source_type),
    is_reference = if_else(is.na(is_reference), TRUE, is_reference),
    habitat = if_else(is.na(habitat), "reference", habitat),
    sample_id = if_else(is_coldseep_target, sub("-k141_.*$", "", normalized_id), NA_character_),
    MAG_id = if_else(is.na(MAG_id), NA_character_, MAG_id),
    display_label = case_when(
      is_coldseep_target ~ normalized_id,
      is_reference & grepl("^ISMEJFig5\\|", protein_id) ~ sub("^ISMEJFig5\\|", "", protein_id),
      is_reference & grepl("^MicrobiologySpectrum_", protein_id) ~ sub("^MicrobiologySpectrum_", "", protein_id),
      TRUE ~ protein_id
    )
  )

feature_mat <- t(vapply(seq_tbl$sequence, calc_features, numeric(1 + 20 + 400)))
feature_df <- as.data.frame(feature_mat)
feature_df$protein_id <- seq_tbl$protein_id
feature_df <- feature_df %>% relocate(protein_id)

analysis_df <- metadata %>%
  filter(length_aa >= 80, branch != "other DMSOR-family reference-associated")

target_analysis_df <- analysis_df %>%
  filter(is_coldseep_target, branch %in% c("canonical aioA-associated", "IdrA-associated",
                                           "aioA-like-associated",
                                           "unknown AioA-like / uncertain DMSOR"))

feature_analysis <- feature_df %>%
  filter(protein_id %in% analysis_df$protein_id)

x <- feature_analysis %>% select(-protein_id) %>% as.matrix()
rownames(x) <- feature_analysis$protein_id
x_scaled <- scale(x)
x_scaled[is.na(x_scaled)] <- 0

pca <- prcomp(x_scaled, center = FALSE, scale. = FALSE)
pca_scores <- as.data.frame(pca$x[, 1:min(10, ncol(pca$x)), drop = FALSE])
pca_scores$protein_id <- rownames(pca$x)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)

dist_mat <- dist(x_scaled, method = "euclidean")
pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)
pcoa_scores <- data.frame(
  protein_id = rownames(x_scaled),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
)
pcoa_var <- pcoa$eig / sum(abs(pcoa$eig))

has_umap <- requireNamespace("umap", quietly = TRUE)
if (has_umap) {
  umap_cfg <- umap::umap.defaults
  umap_cfg$n_neighbors <- min(15, nrow(x_scaled) - 1)
  umap_cfg$min_dist <- 0.2
  umap_cfg$metric <- "euclidean"
  umap_res <- umap::umap(x_scaled, config = umap_cfg)
  umap_scores <- data.frame(
    protein_id = rownames(x_scaled),
    UMAP1 = umap_res$layout[, 1],
    UMAP2 = umap_res$layout[, 2],
    stringsAsFactors = FALSE
  )
} else {
  umap_scores <- data.frame(
    protein_id = character(),
    UMAP1 = numeric(),
    UMAP2 = numeric(),
    note = character()
  )
}

plot_meta <- analysis_df %>%
  select(protein_id, normalized_id, branch, source_type, habitat, sample_id, sample_group, MAG_id,
         target_gene_id, is_reference, is_coldseep_target, display_label, length_aa)

pca_plot_df <- pca_scores %>% left_join(plot_meta, by = "protein_id")
pcoa_plot_df <- pcoa_scores %>% left_join(plot_meta, by = "protein_id")
umap_plot_df <- umap_scores %>% left_join(plot_meta, by = "protein_id")

branch_levels <- c("canonical aioA-associated", "IdrA-associated",
                   "aioA-like-associated", "unknown AioA-like / uncertain DMSOR")
branch_colors <- c(
  "canonical aioA-associated" = "#D95F02",
  "IdrA-associated" = "#1F78B4",
  "aioA-like-associated" = "#009E73",
  "unknown AioA-like / uncertain DMSOR" = "#CC79A7"
)

shape_values <- c(
  "contig-only" = 21,
  "MAG-derived contig" = 22,
  "AioA/AioA-like reference" = 24,
  "IdrA reference" = 25,
  "ISMEJ Fig5 reference" = 23,
  "NarH outgroup" = 4,
  "legacy DMSOR-family reference" = 22
)

theme_pub <- function() {
  theme_classic(base_family = "Times New Roman", base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 9.5),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
    )
}

plot_embedding <- function(df, xvar, yvar, title, subtitle, xlab, ylab, label_refs = FALSE) {
  label_df <- df %>%
    filter(is_reference, source_type == "ISMEJ Fig5 reference") %>%
    group_by(branch) %>%
    slice_head(n = 2) %>%
    ungroup()

  p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]],
                      color = branch, shape = source_type, fill = branch)) +
    geom_point(data = df %>% filter(is_reference), size = 1.7, alpha = 0.28, stroke = 0.25) +
    geom_point(data = df %>% filter(!is_reference), size = 2.6, alpha = 0.9, stroke = 0.45) +
    stat_ellipse(data = df %>% filter(!is_reference), aes(group = branch), linewidth = 0.45,
                 linetype = "dashed", alpha = 0.7, show.legend = FALSE) +
    scale_color_manual(values = branch_colors, breaks = branch_levels, drop = FALSE) +
    scale_fill_manual(values = branch_colors, breaks = branch_levels, drop = FALSE) +
    scale_shape_manual(values = shape_values) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab,
         color = "Branch", fill = "Branch", shape = "Source") +
    theme_pub()
  if (label_refs && nrow(label_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      aes(label = display_label),
      size = 2.3,
      max.overlaps = 40,
      min.segment.length = 0,
      segment.size = 0.2,
      show.legend = FALSE
    )
  }
  p
}

p_pca <- plot_embedding(
  pca_plot_df, "PC1", "PC2",
  "K08356-related protein sequence-space PCA",
  "230 ColdSeep contig targets plus contextual references",
  sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
  sprintf("PC2 (%.1f%%)", 100 * pca_var[2])
)

p_pcoa <- plot_embedding(
  pcoa_plot_df, "PCoA1", "PCoA2",
  "K08356-related protein sequence-space PCoA",
  "Scaled sequence-feature distance; target ellipses shown",
  sprintf("PCoA1 (%.1f%%)", 100 * pcoa_var[1]),
  sprintf("PCoA2 (%.1f%%)", 100 * pcoa_var[2])
)

if (has_umap && nrow(umap_plot_df) > 0) {
  p_umap <- plot_embedding(
    umap_plot_df, "UMAP1", "UMAP2",
    "K08356-related protein sequence-space UMAP",
    "UMAP on scaled sequence-feature embedding",
    "UMAP1", "UMAP2"
  )
} else {
  p_umap <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = "UMAP package not available to Rscript.\nPCoA/MDS was generated as the non-linear-distance fallback.",
             family = "Times New Roman", size = 5) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = "UMAP not generated in this run") +
    theme_void(base_family = "Times New Roman")
}

target_ids <- target_analysis_df$protein_id
x_target <- x_scaled[rownames(x_scaled) %in% target_ids, , drop = FALSE]
target_meta <- plot_meta %>% filter(protein_id %in% rownames(x_target))
target_dist <- dist(x_target, method = "euclidean")

centroids <- as.data.frame(x_target) %>%
  mutate(protein_id = rownames(x_target)) %>%
  left_join(target_meta %>% select(protein_id, branch), by = "protein_id") %>%
  group_by(branch) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

centroid_mat <- centroids %>% select(-branch) %>% as.matrix()
rownames(centroid_mat) <- centroids$branch
centroid_dist <- as.matrix(dist(centroid_mat, method = "euclidean"))
centroid_long <- as.data.frame(as.table(centroid_dist)) %>%
  rename(branch_1 = Var1, branch_2 = Var2, centroid_distance = Freq)

p_centroid <- ggplot(centroid_long, aes(branch_1, branch_2, fill = centroid_distance)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", centroid_distance)), family = "Times New Roman", size = 3.5) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC") +
  labs(title = "Branch centroid distance in sequence feature space",
       x = NULL, y = NULL, fill = "Distance") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

sil <- silhouette(as.integer(factor(target_meta$branch[match(rownames(x_target), target_meta$protein_id)])),
                  target_dist)
sil_df <- as.data.frame(sil[, 1:3])
sil_df$protein_id <- rownames(x_target)
sil_df <- sil_df %>%
  rename(cluster = cluster, neighbor = neighbor, silhouette_width = sil_width) %>%
  left_join(target_meta, by = "protein_id")

sil_summary <- sil_df %>%
  group_by(branch) %>%
  summarise(n = n(),
            mean_silhouette = mean(silhouette_width),
            median_silhouette = median(silhouette_width),
            .groups = "drop") %>%
  arrange(desc(mean_silhouette))

set.seed(20260708)
permanova <- adonis2(dist_mat ~ branch, data = plot_meta, permutations = 999)
permanova_targets <- adonis2(target_dist ~ branch, data = target_meta, permutations = 999)
permanova_df <- as.data.frame(permanova)
permanova_df$term <- rownames(permanova_df)
permanova_df <- permanova_df %>% relocate(term) %>% mutate(analysis_scope = "targets_plus_references")
permanova_targets_df <- as.data.frame(permanova_targets)
permanova_targets_df$term <- rownames(permanova_targets_df)
permanova_targets_df <- permanova_targets_df %>% relocate(term) %>% mutate(analysis_scope = "coldseep_targets_only")
permanova_df <- bind_rows(permanova_targets_df, permanova_df)

branch_counts <- metadata %>%
  count(branch, source_type, name = "n") %>%
  arrange(branch, source_type)

matching_check <- data.frame(
  metric = c("fasta_sequences_total", "coldseep_contig_targets_in_fasta",
             "coldseep_targets_with_tree11_clade_label", "tree11_clade_targets_in_table",
             "MAG_targets_in_label_table", "MAG_targets_found_in_contig_fasta",
             "analysis_sequences_after_filter", "umap_available_to_Rscript"),
  value = c(nrow(seq_tbl), sum(seq_tbl$is_coldseep_target),
            sum(seq_tbl$normalized_id %in% clade_tbl$normalized_id), nrow(clade_tbl),
            nrow(mag_label_tbl), sum(mag_label_tbl$mag_contig_normalized_id %in% seq_tbl$normalized_id),
            nrow(analysis_df), as.character(has_umap))
)

write_csv(metadata %>% select(-sequence), file.path(root_dir, "inputs", "protein_embedding_metadata.csv"))
write_csv(feature_df, file.path(root_dir, "results", "protein_sequence_feature_embedding_matrix.csv"))
write_csv(pca_plot_df, file.path(root_dir, "results", "protein_sequence_PCA_scores.csv"))
write_csv(pcoa_plot_df, file.path(root_dir, "results", "protein_sequence_PCoA_scores.csv"))
write_csv(umap_plot_df, file.path(root_dir, "results", "protein_sequence_UMAP_scores.csv"))
write_csv(centroid_long, file.path(root_dir, "results", "branch_centroid_distance_long.csv"))
write_csv(as.data.frame(centroid_dist) %>% mutate(branch = rownames(centroid_dist)) %>% relocate(branch),
          file.path(root_dir, "results", "branch_centroid_distance_matrix.csv"))
write_csv(sil_df, file.path(root_dir, "results", "silhouette_by_protein.csv"))
write_csv(sil_summary, file.path(root_dir, "results", "silhouette_summary_by_branch.csv"))
write_csv(permanova_df, file.path(root_dir, "results", "PERMANOVA_embedding_distance_by_branch.csv"))
write_csv(branch_counts, file.path(root_dir, "results", "sequence_counts_by_branch_source.csv"))
write_csv(matching_check, file.path(root_dir, "results", "input_matching_check.csv"))

save_plot <- function(plot, name, width = 7.2, height = 6.2) {
  ggsave(file.path(root_dir, "figures", paste0(name, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  grDevices::svg(file.path(root_dir, "figures", paste0(name, ".svg")), width = width, height = height,
                 family = "Times New Roman")
  print(plot)
  grDevices::dev.off()
  ggsave(file.path(root_dir, "figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 450)
}

save_plot(p_pca, "sequence_feature_PCA_K08356_branches")
save_plot(p_pcoa, "sequence_feature_PCoA_K08356_branches")
save_plot(p_umap, "sequence_feature_UMAP_K08356_branches")
save_plot(p_centroid, "branch_centroid_distance_heatmap", width = 6.8, height = 5.8)

readme <- c(
  "# Protein Sequence-Space Validation v1",
  "",
  "This run validates K08356-related DMSOR branch separation using a reproducible local sequence-feature embedding.",
  "",
  "## Inputs",
  paste0("- FASTA: ", fasta_path),
  paste0("- tree11 contig clade table: ", clade_detail_path),
  paste0("- MAG branch label table for MAG-derived flags: ", mag_label_path),
  "",
  "## Method",
  "- Protein sequences were represented by amino-acid composition plus dipeptide frequency features.",
  "- Features were scaled before PCA, PCoA/MDS, centroid distance, silhouette, and PERMANOVA analyses.",
  "- The main quantitative summaries use the 230 ColdSeep contig K08356 targets with tree11 clade labels.",
  "- Reference sequences are retained in PCA/PCoA/UMAP plots as contextual anchors.",
  "- This is a local sequence-space validation, not an ESM-2 protein language model run.",
  "- If the R `umap` package is available to Rscript, UMAP is generated automatically. Otherwise the UMAP figure states that PCoA/MDS is the fallback.",
  "",
  "## Interpretation boundary",
  "- This analysis can support whether branches occupy distinguishable protein sequence-feature space.",
  "- It cannot directly prove substrate specificity, catalytic activity, or enzyme function.",
  "- The current run uses the tree11 contig-level K08356 target set and clade labels.",
  "- The 48 MAG-derived K08356 labels are used only to flag which contig targets are also MAG-derived.",
  "- ESM-2/ProtT5 can be added later as an upgraded protein language model embedding, but was not run in this local v1.",
  "",
  "## Key outputs",
  "- protein_sequence_feature_embedding_matrix.csv",
  "- protein_sequence_PCA_scores.csv",
  "- protein_sequence_UMAP_scores.csv if available",
  "- branch_centroid_distance_matrix.csv",
  "- silhouette_summary_by_branch.csv",
  "- PERMANOVA_embedding_distance_by_branch.csv",
  "",
  paste0("UMAP available to Rscript in this run: ", has_umap)
)
writeLines(readme, file.path(root_dir, "README_interpretation.md"))

cat("Protein sequence-space validation complete.\n")
cat("FASTA sequences:", nrow(seq_tbl), "\n")
cat("ColdSeep contig targets:", sum(seq_tbl$is_coldseep_target), "\n")
cat("Tree11 clade targets matched:", sum(seq_tbl$normalized_id %in% clade_tbl$normalized_id), "of", nrow(clade_tbl), "\n")
cat("MAG targets matched to contig fasta:", sum(mag_label_tbl$mag_contig_normalized_id %in% seq_tbl$normalized_id), "of", nrow(mag_label_tbl), "\n")
cat("UMAP available:", has_umap, "\n")

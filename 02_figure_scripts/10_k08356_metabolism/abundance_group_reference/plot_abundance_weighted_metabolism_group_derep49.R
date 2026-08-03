suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(script_arg)) {
  script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_arg[1]))
  dirname(normalizePath(script_path))
} else {
  getwd()
}
analysis_dir <- if (file.exists(file.path(
  script_dir, "results", "K08356seq49_MAG48_metabolic_feature_scores.tsv"
))) script_dir else dirname(script_dir)

abundance_path <- Sys.getenv(
  "GROUP_TPM_LONG",
  file.path(
    analysis_dir, "abundance_group_reference_20260802", "raw",
    "IS_AS_ES_NS_MAG_TPM_long.tsv"
  )
)
run_summary_path <- Sys.getenv(
  "GROUP_TPM_SUMMARY",
  file.path(
    analysis_dir, "abundance_group_reference_20260802", "raw",
    "IS_AS_ES_NS_MAG_TPM_summary.tsv"
  )
)
score_path <- Sys.getenv(
  "METABOLIC_SCORE_TSV",
  file.path(
    analysis_dir, "results", "K08356seq49_MAG48_metabolic_feature_scores.tsv"
  )
)
mapping_path <- Sys.getenv(
  "K08356_MAPPING_TSV",
  file.path(
    analysis_dir, "results", "K08356seq49_MAG48_clade_habitat_mapping.tsv"
  )
)
sample_metadata_path <- Sys.getenv(
  "SAMPLE_METADATA_CSV",
  file.path(
    dirname(dirname(dirname(analysis_dir))),
    "TPM_combo_figure_20260728", "sample_metadata_used.csv"
  )
)
out_dir <- Sys.getenv(
  "ABUNDANCE_METABOLISM_OUT_DIR",
  file.path(analysis_dir, "abundance_group_reference_20260802", "results")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)

required <- c(
  abundance_path, run_summary_path, score_path, mapping_path,
  sample_metadata_path
)
if (any(!file.exists(required))) {
  stop("Missing input: ", paste(required[!file.exists(required)], collapse = ", "))
}

abundance_all <- read.delim(abundance_path, check.names = FALSE)
run_summary <- read.delim(run_summary_path, check.names = FALSE)
scores <- read.delim(score_path, check.names = FALSE)
mapping <- read.delim(mapping_path, check.names = FALSE)
sample_metadata <- read.csv(sample_metadata_path, check.names = FALSE) %>%
  transmute(
    Sample = sample,
    sample_habitat = habitat,
    depth_mid
  )

habitat_levels <- c("IS", "AS", "ES", "NS")
clade_levels <- paste("Clade", 1:4)
category_levels <- c(
  "Carbon", "Nitrogen", "Sulfur", "Arsenic", "Vitamin B12", "Motility"
)

membership <- mapping %>%
  transmute(
    Sequence_ID = contig_gene,
    MAG,
    source_habitat = habitat,
    final_clade = as.character(final_clade),
    final_clade_name
  ) %>%
  distinct()
target_MAGs <- sort(unique(membership$MAG))

if (nrow(membership) != 49L || length(target_MAGs) != 48L ||
    nrow(scores) != 1920L || n_distinct(scores$MAG) != 48L ||
    n_distinct(scores$feature_name) != 40L ||
    !setequal(target_MAGs, scores$MAG)) {
  stop("Expected the locked 49-sequence/48-MAG/40-feature analysis set.")
}

abundance <- abundance_all %>%
  filter(Genome %in% target_MAGs) %>%
  transmute(
    reference_group = Group,
    Sample,
    MAG = Genome,
    TPM = as.numeric(TPM),
    Source_file
  ) %>%
  left_join(sample_metadata, by = "Sample")

MAG_reference_group <- abundance %>%
  distinct(MAG, reference_group) %>%
  left_join(
    membership %>% distinct(MAG, source_habitat),
    by = "MAG"
  )
if (nrow(abundance) != 56L * 48L ||
    n_distinct(abundance$Sample) != 56L ||
    n_distinct(abundance$MAG) != 48L ||
    anyDuplicated(abundance[c("Sample", "MAG")]) ||
    anyNA(abundance$TPM) || any(abundance$TPM < 0) ||
    anyNA(abundance$sample_habitat) ||
    any(table(abundance$MAG) != 56L) ||
    nrow(MAG_reference_group) != 48L ||
    any(MAG_reference_group$reference_group !=
          MAG_reference_group$source_habitat)) {
  stop("Filtered group-reference abundance failed exact 56 x 48 validation.")
}

expected_reference_sizes <- c(IS = 315L, AS = 128L, ES = 252L, NS = 156L)
observed_reference_sizes <- setNames(
  run_summary$Genome_rows_per_file,
  run_summary$Group
)
if (!identical(
  as.integer(observed_reference_sizes[names(expected_reference_sizes)]),
  as.integer(expected_reference_sizes)
)) {
  stop("Unexpected group-reference genome counts.")
}

clade_counts <- membership %>%
  group_by(final_clade, final_clade_name) %>%
  summarise(
    n_sequences = n_distinct(Sequence_ID),
    n_unique_MAG = n_distinct(MAG),
    .groups = "drop"
  ) %>%
  arrange(factor(final_clade, levels = clade_levels))
if (!identical(clade_counts$n_sequences, c(3L, 7L, 19L, 20L))) {
  stop("Expected clade sequence counts 3/7/19/20.")
}

score_membership <- scores %>%
  inner_join(
    membership %>%
      select(MAG, final_clade, final_clade_name) %>%
      distinct(),
    by = "MAG",
    relationship = "many-to-many"
  )

calculate_summary <- function(exclude_dual_MAG = FALSE) {
  current_abundance <- abundance
  current_scores <- score_membership
  if (exclude_dual_MAG) {
    current_abundance <- current_abundance %>%
      filter(MAG != "S1_9-12_bin1")
    current_scores <- current_scores %>%
      filter(MAG != "S1_9-12_bin1")
  }

  sample_feature <- current_abundance %>%
    inner_join(
      current_scores,
      by = "MAG",
      relationship = "many-to-many"
    ) %>%
    group_by(
      Sample, sample_habitat, final_clade, final_clade_name,
      research_category, feature_name, source_type
    ) %>%
    summarise(
      n_unique_MAG = n_distinct(MAG),
      total_host_TPM = sum(TPM),
      functional_host_TPM_50 = sum(TPM[score_ge_50]),
      score_TPM_numerator = sum(TPM * reconstruction_score_pct),
      any_host_detected = total_host_TPM > 0,
      functional_host_detected_50 = functional_host_TPM_50 > 0,
      .groups = "drop"
    ) %>%
    mutate(
      sample_TPM_weighted_score_pct = if_else(
        any_host_detected,
        score_TPM_numerator / total_host_TPM,
        NA_real_
      )
    )

  habitat_summary <- sample_feature %>%
    group_by(
      sample_habitat, final_clade, final_clade_name,
      research_category, feature_name, source_type
    ) %>%
    summarise(
      n_samples = n_distinct(Sample),
      n_unique_MAG = first(n_unique_MAG),
      total_host_TPM_all_samples = sum(total_host_TPM),
      mean_total_host_TPM_per_sample = mean(total_host_TPM),
      median_total_host_TPM_per_sample = median(total_host_TPM),
      mean_functional_host_TPM_50_per_sample = mean(functional_host_TPM_50),
      n_samples_host_detected = sum(any_host_detected),
      sample_host_detection_prevalence_pct = 100 * mean(any_host_detected),
      sample_functional_detection_prevalence_50_pct =
        100 * mean(functional_host_detected_50),
      abundance_weighted_score_pct = if_else(
        any(any_host_detected),
        mean(sample_TPM_weighted_score_pct, na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      dual_MAG_excluded = exclude_dual_MAG,
      reference_bias_note = paste(
        "Descriptive only: target MAGs were quantified within four different",
        "background references (IS 315, AS 128, ES 252, NS 156 MAGs)"
      ),
      score_aggregation_note = paste(
        "Each sample is equally weighted: calculate sum(TPM*score)/sum(TPM)",
        "within detected hosts per sample, then average non-missing sample scores",
        "within habitat; host-nondetection samples are NA for fill"
      )
    )

  list(sample_feature = sample_feature, habitat_summary = habitat_summary)
}

full <- calculate_summary(FALSE)
exclude_dual <- calculate_summary(TRUE)

dual_sensitivity <- full$habitat_summary %>%
  select(
    sample_habitat, final_clade, research_category, feature_name,
    full_weighted_score_pct = abundance_weighted_score_pct,
    full_detection_prevalence_50_pct =
      sample_functional_detection_prevalence_50_pct
  ) %>%
  left_join(
    exclude_dual$habitat_summary %>%
      select(
        sample_habitat, final_clade, research_category, feature_name,
        excluded_weighted_score_pct = abundance_weighted_score_pct,
        excluded_detection_prevalence_50_pct =
          sample_functional_detection_prevalence_50_pct
      ),
    by = c(
      "sample_habitat", "final_clade", "research_category", "feature_name"
    )
  ) %>%
  mutate(
    weighted_score_delta_pct =
      excluded_weighted_score_pct - full_weighted_score_pct,
    detection_prevalence_delta_pct =
      excluded_detection_prevalence_50_pct -
      full_detection_prevalence_50_pct
  )

clade_reference_audit <- membership %>%
  distinct(MAG, final_clade) %>%
  left_join(
    MAG_reference_group %>% select(MAG, reference_group),
    by = "MAG"
  ) %>%
  group_by(final_clade) %>%
  summarise(
    n_unique_MAG = n_distinct(MAG),
    n_reference_groups = n_distinct(reference_group),
    reference_groups = paste(sort(unique(reference_group)), collapse = ";"),
    reference_consistent_for_cross_sample_habitat_comparison =
      n_reference_groups == 1L,
    .groups = "drop"
  )

clade4_sample_abundance <- abundance %>%
  inner_join(
    membership %>%
      filter(final_clade == "Clade 4") %>%
      distinct(MAG, final_clade),
    by = "MAG"
  ) %>%
  group_by(Sample, sample_habitat, final_clade) %>%
  summarise(
    n_Clade4_MAG = n_distinct(MAG),
    Clade4_host_MAG_TPM = sum(TPM),
    .groups = "drop"
  ) %>%
  mutate(
    sample_habitat = factor(sample_habitat, levels = habitat_levels),
    log10_TPM_plus_1 = log10(Clade4_host_MAG_TPM + 1)
  )
if (nrow(clade4_sample_abundance) != 56L ||
    any(clade4_sample_abundance$n_Clade4_MAG != 20L) ||
    clade_reference_audit$n_reference_groups[
      clade_reference_audit$final_clade == "Clade 4"
    ] != 1L) {
  stop("Clade 4 must contain 20 MAGs quantified against one reference.")
}

clade4_KW <- suppressWarnings(
  kruskal.test(log10_TPM_plus_1 ~ sample_habitat,
               data = clade4_sample_abundance)
)
habitat_pairs <- combn(habitat_levels, 2, simplify = FALSE)
clade4_pairwise <- bind_rows(lapply(habitat_pairs, function(pair) {
  x <- clade4_sample_abundance$log10_TPM_plus_1[
    clade4_sample_abundance$sample_habitat == pair[1]
  ]
  y <- clade4_sample_abundance$log10_TPM_plus_1[
    clade4_sample_abundance$sample_habitat == pair[2]
  ]
  identical_distributions <- length(unique(c(x, y))) < 2L
  p_value <- if (identical_distributions) {
    1
  } else {
    suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
  }
  tibble(
    habitat_1 = pair[1],
    habitat_2 = pair[2],
    Wilcoxon_p = p_value,
    test_note = if_else(
      identical_distributions,
      "Both groups have identical values; p set to 1",
      "Two-sided Wilcoxon rank-sum test"
    )
  )
})) %>%
  mutate(Wilcoxon_BH = p.adjust(Wilcoxon_p, method = "BH"))
clade4_stats <- tibble(
  test = "Kruskal-Wallis on log10(TPM+1)",
  statistic = unname(clade4_KW$statistic),
  df = unname(clade4_KW$parameter),
  p_value = clade4_KW$p.value,
  scope = paste(
    "20 Clade 4 MAGs; all samples mapped against the same IS group reference;",
    "station/depth dependence not modeled"
  )
)

clade4_detection_by_habitat <- clade4_sample_abundance %>%
  mutate(detected = Clade4_host_MAG_TPM > 0) %>%
  group_by(sample_habitat) %>%
  summarise(
    n_samples = n(),
    n_detected = sum(detected),
    n_not_detected = sum(!detected),
    detection_prevalence_pct = 100 * mean(detected),
    .groups = "drop"
  )

clade4_IS_nonIS <- clade4_sample_abundance %>%
  mutate(
    habitat_binary = if_else(sample_habitat == "IS", "IS", "non-IS"),
    detected = Clade4_host_MAG_TPM > 0
  ) %>%
  count(habitat_binary, detected) %>%
  complete(
    habitat_binary = c("IS", "non-IS"),
    detected = c(FALSE, TRUE),
    fill = list(n = 0L)
  )

count_value <- function(habitat, detected_value) {
  clade4_IS_nonIS$n[
    clade4_IS_nonIS$habitat_binary == habitat &
      clade4_IS_nonIS$detected == detected_value
  ]
}
fisher_matrix <- matrix(
  c(
    count_value("IS", TRUE), count_value("IS", FALSE),
    count_value("non-IS", TRUE), count_value("non-IS", FALSE)
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    habitat = c("IS", "non-IS"),
    status = c("detected", "not_detected")
  )
)
if (!identical(as.integer(fisher_matrix), c(20L, 1L, 1L, 34L))) {
  stop("Unexpected Clade 4 IS/non-IS detection table.")
}
fisher_p <- fisher.test(fisher_matrix)$p.value

# Conditional maximum-likelihood odds ratio and central exact interval.
conditional_support <- seq.int(
  max(0, sum(fisher_matrix[1, ]) - sum(fisher_matrix[, 2])),
  min(sum(fisher_matrix[1, ]), sum(fisher_matrix[, 1]))
)
conditional_log_base <-
  lchoose(sum(fisher_matrix[, 1]), conditional_support) +
  lchoose(
    sum(fisher_matrix[, 2]),
    sum(fisher_matrix[1, ]) - conditional_support
  )
conditional_probabilities <- function(log_theta) {
  log_weight <- conditional_log_base + conditional_support * log_theta
  weight <- exp(log_weight - max(log_weight))
  weight / sum(weight)
}
conditional_mean <- function(log_theta) {
  sum(conditional_support * conditional_probabilities(log_theta))
}
solve_increasing <- function(fn, target) {
  exp(uniroot(function(log_theta) fn(log_theta) - target,
              interval = c(-50, 50))$root)
}
observed_a <- fisher_matrix["IS", "detected"]
conditional_OR <- solve_increasing(conditional_mean, observed_a)
conditional_CI_lower <- solve_increasing(
  function(log_theta) {
    probabilities <- conditional_probabilities(log_theta)
    sum(probabilities[conditional_support >= observed_a])
  },
  0.025
)
conditional_CI_upper <- solve_increasing(
  function(log_theta) {
    probabilities <- conditional_probabilities(log_theta)
    sum(probabilities[conditional_support > observed_a])
  },
  0.975
)

clade4_fisher <- tibble(
  IS_detected = fisher_matrix["IS", "detected"],
  IS_not_detected = fisher_matrix["IS", "not_detected"],
  nonIS_detected = fisher_matrix["non-IS", "detected"],
  nonIS_not_detected = fisher_matrix["non-IS", "not_detected"],
  fisher_exact_two_sided_p = fisher_p,
  conditional_odds_ratio = conditional_OR,
  conditional_OR_95CI_lower = conditional_CI_lower,
  conditional_OR_95CI_upper = conditional_CI_upper,
  interpretation = paste(
    "Strongly IS-associated, not IS-exclusive; MAG detection is TPM >0 after",
    "CoverM filters (>=10% covered fraction, >=95% read identity,",
    ">=75% aligned-read fraction)"
  )
)

feature_order <- scores %>%
  distinct(research_category, feature_name) %>%
  mutate(
    research_category = factor(research_category, levels = category_levels)
  ) %>%
  arrange(research_category, feature_name) %>%
  transmute(
    plot_label = paste(
      research_category,
      str_wrap(feature_name, width = 49),
      sep = " | "
    )
  ) %>%
  pull(plot_label)

sample_counts <- sample_metadata %>%
  count(sample_habitat, name = "n_samples") %>%
  mutate(
    sample_habitat = factor(sample_habitat, levels = habitat_levels),
    axis_label = paste0(sample_habitat, "\n(n=", n_samples, ")")
  )
habitat_axis_labels <- setNames(
  sample_counts$axis_label,
  as.character(sample_counts$sample_habitat)
)

clade_labels <- setNames(
  paste0(
    clade_counts$final_clade,
    " (", clade_counts$n_sequences, " seq / ",
    clade_counts$n_unique_MAG, " MAGs)\n",
    str_wrap(clade_counts$final_clade_name, 31)
  ),
  clade_counts$final_clade
)

plot_data <- full$habitat_summary %>%
  mutate(
    sample_habitat = factor(sample_habitat, levels = habitat_levels),
    final_clade = factor(final_clade, levels = clade_levels),
    plot_label = paste(
      research_category,
      str_wrap(feature_name, width = 49),
      sep = " | "
    ),
    plot_label = factor(plot_label, levels = rev(feature_order)),
    zero_host_detection = sample_host_detection_prevalence_pct == 0
  )

caption_text <- paste(
  "Fill: first calculate sum(TPM*score)/sum(TPM) within each sample, then equally average detected-host samples within habitat; host-nondetection samples are NA/blank.",
  "Size: sample detection prevalence of host MAGs with >=50% partial reconstruction.",
  "Host detection is TPM >0 after CoverM filters: >=10% covered fraction, >=95% read identity, >=75% aligned-read fraction, and 0.1/0.9 end trimming.",
  "Arsenic rows use independent non-focal binary markers; other rows use the predefined legacy gene-set coverage.",
  "S1_9-12_bin1 contributes to Clades 1 and 3, so clade categories are non-exclusive; an exclusion sensitivity table is provided.",
  "Important: Clades 1-3 mix target MAGs quantified in different background references (315/128/252/156 MAGs) and are descriptive.",
  "All 20 Clade 4 MAGs use the same IS reference across 56 samples, so its companion habitat comparison is reference-consistent but still lacks station/depth adjustment.",
  "Do not use this plot alone to claim a habitat-driven mechanism or spatial heterogeneity.",
  "IS, incipient seep; AS, active seep; ES, extinct seep; NS, non-seep background."
)

p <- ggplot(
  plot_data,
  aes(
    x = sample_habitat,
    y = plot_label,
    size = sample_functional_detection_prevalence_50_pct,
    fill = abundance_weighted_score_pct
  )
) +
  geom_point(
    shape = 21,
    color = "#263238",
    stroke = 0.28,
    alpha = 0.94
  ) +
  geom_point(
    data = plot_data %>% filter(zero_host_detection),
    shape = 21,
    size = 1.25,
    fill = "white",
    color = "#8a8a8a",
    stroke = 0.35,
    inherit.aes = FALSE,
    aes(x = sample_habitat, y = plot_label)
  ) +
  facet_grid(
    . ~ final_clade,
    labeller = labeller(final_clade = clade_labels),
    scales = "fixed",
    space = "fixed"
  ) +
  scale_x_discrete(labels = habitat_axis_labels, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    na.value = "white",
    name = "Mean within-sample\nTPM-weighted score (%)"
  ) +
  scale_size_continuous(
    range = c(1.2, 6.2),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = paste0(
      "Sample detection prevalence\nof host MAGs with ",
      "\u226550% partial reconstruction (%)"
    )
  ) +
  labs(
    title = "Abundance-weighted metabolic potential of group-dereplicated K08356 MAGs",
    subtitle = paste(
      "49 K08356 sequences / 48 unique MAGs / 56 samples;",
      "clades 3/7/19/20; descriptive group-reference TPM analysis"
    ),
    x = "Sample habitat",
    y = "Metabolic feature",
    caption = str_wrap(caption_text, width = 185)
  ) +
  theme_bw(base_family = "Arial", base_size = 8.2) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 8.5, color = "#37474f"),
    plot.caption = element_text(size = 6.1, lineheight = 1.05, hjust = 0),
    strip.background = element_rect(fill = "#eceff1", color = "#59636a"),
    strip.text = element_text(face = "bold", size = 7.1, lineheight = 0.95),
    panel.grid.major = element_line(color = "#e4e8eb", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.06, "lines"),
    axis.text.x = element_text(face = "bold", size = 7.1),
    axis.text.y = element_text(size = 5.25, lineheight = 0.88),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    plot.margin = margin(8, 8, 8, 8)
  )

base_name <- "group_derep49_MAG48_abundance_weighted_metabolism"
ggsave(
  file.path(out_dir, paste0(base_name, ".png")),
  p, width = 15.8, height = 11.6, units = "in", dpi = 320, bg = "white"
)
ggsave(
  file.path(out_dir, paste0(base_name, ".pdf")),
  p, width = 15.8, height = 11.6, units = "in", device = cairo_pdf,
  bg = "white"
)

habitat_colors <- c(
  IS = "#3aa88f", AS = "#ed9b52", ES = "#758fb8", NS = "#df705b"
)
p_clade4 <- ggplot(
  clade4_sample_abundance,
  aes(x = sample_habitat, y = log10_TPM_plus_1, fill = sample_habitat)
) +
  geom_boxplot(
    width = 0.58, outlier.shape = NA, alpha = 0.78,
    color = "#37474f", linewidth = 0.45
  ) +
  geom_jitter(
    aes(color = sample_habitat),
    width = 0.12, size = 2.0, alpha = 0.82, show.legend = FALSE
  ) +
  scale_fill_manual(values = habitat_colors, drop = FALSE) +
  scale_color_manual(values = habitat_colors, drop = FALSE) +
  annotate(
    "label", x = 4.45, y = Inf,
    label = paste0(
      "Kruskal-Wallis p ",
      if_else(clade4_KW$p.value < 0.001, "< 0.001",
              paste0("= ", formatC(clade4_KW$p.value, digits = 3,
                                    format = "f"))),
      "\nIS vs non-IS Fisher p = ",
      format(fisher_p, scientific = TRUE, digits = 3),
      "\nConditional OR = ", formatC(conditional_OR, digits = 1, format = "f"),
      " (95% CI ", formatC(conditional_CI_lower, digits = 1, format = "f"),
      "-", formatC(conditional_CI_upper, digits = 1, format = "f"), ")"
    ),
    hjust = 1, vjust = 1.3, size = 3.2, fill = "white"
  ) +
  labs(
    title = "Clade 4 host-MAG abundance across four sample habitats",
    subtitle = paste(
      "20 Clade 4 MAGs quantified in all 56 samples against one IS reference;",
      "reference-consistent sample comparison"
    ),
    x = "Sample habitat",
    y = expression(log[10]("summed host-MAG TPM" + 1)),
    caption = str_wrap(
      paste(
        "Points are samples. This is a sample-level association, not evidence of activity or causality.",
        "Detection is TPM >0 after >=10% covered fraction, >=95% read identity and >=75% aligned-read fraction.",
        "Station/core and depth are not adjusted in this provisional nonparametric test."
      ),
      width = 118
    )
  ) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "#455a64"),
    plot.caption = element_text(size = 7, lineheight = 1.05, hjust = 0),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none",
    plot.margin = margin(10, 14, 10, 10)
  )

ggsave(
  file.path(out_dir, "Clade4_same_reference_host_MAG_TPM_by_habitat.png"),
  p_clade4, width = 8.6, height = 5.7, units = "in", dpi = 320, bg = "white"
)
ggsave(
  file.path(out_dir, "Clade4_same_reference_host_MAG_TPM_by_habitat.pdf"),
  p_clade4, width = 8.6, height = 5.7, units = "in", device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(out_dir, "Clade4_same_reference_host_MAG_TPM_by_habitat.svg"),
  p_clade4, width = 8.6, height = 5.7, units = "in", device = svg,
  bg = "white"
)
ggsave(
  file.path(out_dir, paste0(base_name, ".svg")),
  p, width = 15.8, height = 11.6, units = "in", device = svg,
  bg = "white"
)

write.table(
  abundance,
  file.path(out_dir, "validated_sample56_MAG48_group_reference_TPM.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  full$sample_feature,
  file.path(out_dir, "sample_clade_feature_abundance_long.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  full$habitat_summary,
  file.path(out_dir, "abundance_weighted_clade_habitat_feature_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  dual_sensitivity,
  file.path(out_dir, "dual_clade_MAG_exclusion_abundance_sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  MAG_reference_group,
  file.path(out_dir, "MAG48_reference_group_ID_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  clade_reference_audit,
  file.path(out_dir, "clade_reference_consistency_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  clade4_sample_abundance,
  file.path(out_dir, "Clade4_same_reference_sample_TPM.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  clade4_stats,
  file.path(out_dir, "Clade4_same_reference_Kruskal_Wallis.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  clade4_pairwise,
  file.path(out_dir, "Clade4_same_reference_pairwise_Wilcoxon_BH.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  clade4_detection_by_habitat,
  file.path(out_dir, "Clade4_detection_prevalence_by_habitat.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  clade4_fisher,
  file.path(out_dir, "Clade4_IS_vs_nonIS_Fisher_exact.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

qc <- c(
  "analysis_status\tdescriptive_group_specific_reference_TPM",
  paste0("input_rows_all_MAGs\t", nrow(abundance_all)),
  paste0("validated_target_rows\t", nrow(abundance)),
  paste0("samples\t", n_distinct(abundance$Sample)),
  paste0("target_MAGs\t", n_distinct(abundance$MAG)),
  paste0("K08356_sequence_memberships\t", nrow(membership)),
  "clade_sequence_counts\t3/7/19/20",
  paste0(
    "group_reference_sizes\t",
    paste(names(observed_reference_sizes), observed_reference_sizes,
          sep = ":", collapse = ";")
  ),
  "same_target_MAG_set_across_samples\tTRUE",
  "same_background_reference_across_source_groups\tFALSE",
  "dual_clade_host\tS1_9-12_bin1",
  "reference_consistent_clade\tClade 4",
  "host_detection_rule\tTPM >0 after CoverM >=10% covered fraction, >=95% read identity, >=75% aligned-read fraction, trim 0.1/0.9",
  paste0("Clade4_IS_vs_nonIS_Fisher_p\t", fisher_p),
  paste0("Clade4_conditional_OR\t", conditional_OR),
  paste0("Clade4_conditional_OR_95CI\t", conditional_CI_lower, ";", conditional_CI_upper),
  "permitted_interpretation\tdescriptive abundance-weighted genomic potential; provisional Clade 4 sample-level habitat association",
  "prohibited_interpretation\tfully adjusted habitat mechanism or spatial heterogeneity"
)
writeLines(qc, file.path(out_dir, "ABUNDANCE_WEIGHTED_METABOLISM_QC.txt"))

cat(paste(qc, collapse = "\n"), "\n")
cat("Output directory:", out_dir, "\n")

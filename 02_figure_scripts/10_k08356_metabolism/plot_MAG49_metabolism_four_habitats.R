suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(script_arg)) {
  script_path <- sub("^--file=", "", script_arg[1])
  script_path <- gsub("~\\+~", " ", script_path)
  dirname(normalizePath(script_path))
} else {
  getwd()
}
font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)

input_xlsx <- Sys.getenv(
  "METABOLIC_RESULT_XLSX",
  file.path(script_dir, "METABOLIC_group_derep49_result.xlsx")
)
classification_tsv <- Sys.getenv(
  "K08356_CLADE_TSV",
  file.path(dirname(script_dir), "four_clade_classification_49.tsv")
)
manifest_tsv <- Sys.getenv(
  "MAG49_MANIFEST_TSV",
  file.path(script_dir, "MAG49_input_manifest.tsv")
)
out_dir <- Sys.getenv(
  "METABOLIC_PLOT_OUT_DIR",
  file.path(script_dir, "results")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(input_xlsx, classification_tsv, manifest_tsv)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "))
}

clade_levels <- paste("Clade", 1:4)
clade_from_display <- c(
  "Canonical AioA" = "Clade 1",
  "Uncertain DMSOR" = "Clade 2",
  "IdrA phylogenetic" = "Clade 3",
  "DIRM-synteny IdrA" = "Clade 4"
)
clade_names <- c(
  "Clade 1" = "Canonical AioA-associated",
  "Clade 2" = "Unknown/uncertain DMSOR",
  "Clade 3" = "IdrA-associated without complete DIRM-like synteny",
  "Clade 4" = "IdrA-associated with DIRM-like synteny"
)

classification <- read.delim(
  classification_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
manifest <- read.delim(
  manifest_tsv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mapping <- classification %>%
  transmute(
    contig_gene = Sequence_ID,
    MAG = MAG_ID,
    habitat_in_classification = Habitat,
    clade_display,
    final_clade = unname(clade_from_display[clade_display]),
    final_clade_name = unname(clade_names[final_clade])
  ) %>%
  left_join(
    manifest %>% select(MAG = bin_id, habitat = group),
    by = "MAG"
  ) %>%
  mutate(
    final_clade = factor(final_clade, levels = clade_levels),
    habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS")),
    habitat_source_match = as.character(habitat) == habitat_in_classification
  )

expected_clade_counts <- c(
  "Clade 1" = 3L,
  "Clade 2" = 7L,
  "Clade 3" = 19L,
  "Clade 4" = 20L
)
observed_clade_counts <- table(mapping$final_clade)
if (nrow(mapping) != 49L || n_distinct(mapping$MAG) != 48L) {
  stop("Expected 49 K08356 sequence memberships from 48 unique MAGs.")
}
if (any(is.na(mapping$final_clade)) || any(is.na(mapping$habitat))) {
  stop("Missing clade or manifest-based habitat assignment.")
}
if (!all(observed_clade_counts[names(expected_clade_counts)] ==
         expected_clade_counts)) {
  stop("Expected clade sequence counts 3/7/19/20.")
}
if (!setequal(unique(mapping$MAG), manifest$bin_id)) {
  stop("The clade classification and 48-MAG input manifest do not match.")
}

normalize_metabolic_mag <- function(x, suffix) {
  x %>%
    str_remove(fixed(suffix)) %>%
    str_replace_all(fixed("."), "-")
}

selected_modules <- tibble::tribble(
  ~research_category, ~feature_name,
  "Carbon", "Acetyl-CoA pathway, CO2 => acetyl-CoA",
  "Carbon", "Citrate cycle (TCA cycle, Krebs cycle)",
  "Carbon", "Citrate cycle, first carbon oxidation, oxaloacetate => 2-oxoglutarate",
  "Carbon", "Citrate cycle, second carbon oxidation, 2-oxoglutarate => oxaloacetate",
  "Carbon", "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate",
  "Carbon", "Glycolysis, core module involving three-carbon compounds",
  "Carbon", "Incomplete reductive citrate cycle, acetyl-CoA => oxoglutarate",
  "Carbon", "Methane oxidation, methanotroph, methane => formaldehyde",
  "Carbon", "Methanogenesis, CO2 => methane",
  "Carbon", "Methanogenesis, acetate => methane",
  "Carbon", "Methanogenesis, methanol => methane",
  "Carbon", "Methanogenesis, methylamine/dimethylamine/trimethylamine => methane",
  "Carbon", "Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)",
  "Carbon", "Reductive citrate cycle (Arnon-Buchanan cycle)",
  "Nitrogen", "Assimilatory nitrate reduction, nitrate => ammonia",
  "Nitrogen", "Complete nitrification, comammox, ammonia => nitrite => nitrate",
  "Nitrogen", "Denitrification, nitrate => nitrogen",
  "Nitrogen", "Dissimilatory nitrate reduction, nitrate => ammonia",
  "Nitrogen", "Nitrification, ammonia => nitrite",
  "Nitrogen", "Nitrogen fixation, nitrogen => ammonia",
  "Nitrogen", "Urea cycle",
  "Sulfur", "Assimilatory sulfate reduction, sulfate => H2S",
  "Sulfur", "Dissimilatory sulfate reduction, sulfate => H2S",
  "Sulfur", "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate",
  "Vitamin B12", "Cobalamin biosynthesis, aerobic, uroporphyrinogen III => precorrin 2 => cobyrinate a,c-diamide",
  "Vitamin B12", "Cobalamin biosynthesis, anaerobic, uroporphyrinogen III => sirohydrochlorin => cobyrinate a,c-diamide",
  "Vitamin B12", "Cobalamin biosynthesis, cobyrinate a,c-diamide => cobalamin"
)

module_step <- read_excel(
  input_xlsx,
  sheet = "KEGGModuleStepHit",
  .name_repair = "minimal"
)
module_suffix <- ".cds.Module.step.presence"
module_mag_columns <- names(module_step)[endsWith(names(module_step), module_suffix)]
if (length(module_mag_columns) != 48L) {
  stop("Expected 48 MAG columns in KEGGModuleStepHit; observed ",
       length(module_mag_columns), ".")
}

available_modules <- unique(module_step$Module)
missing_modules <- setdiff(selected_modules$feature_name, available_modules)
if (length(missing_modules)) {
  stop("Selected KEGG modules absent from METABOLIC output: ",
       paste(missing_modules, collapse = "; "))
}

module_scores <- module_step %>%
  filter(Module %in% selected_modules$feature_name) %>%
  select(Module.step, Module, Module.Category, all_of(module_mag_columns)) %>%
  pivot_longer(
    all_of(module_mag_columns),
    names_to = "metabolic_column",
    values_to = "step_status"
  ) %>%
  mutate(
    MAG = normalize_metabolic_mag(metabolic_column, module_suffix),
    step_present = str_to_lower(as.character(step_status)) == "present"
  ) %>%
  group_by(Module, Module.Category, MAG) %>%
  summarise(
    feature_component_count = n_distinct(Module.step),
    components_detected = n_distinct(Module.step[step_present]),
    reconstruction_score_pct = 100 * components_detected /
      feature_component_count,
    .groups = "drop"
  ) %>%
  inner_join(
    selected_modules,
    by = c("Module" = "feature_name")
  ) %>%
  transmute(
    MAG,
    research_category,
    feature_name = Module,
    source_type = "KEGG module-step coverage",
    feature_component_count,
    components_detected,
    reconstruction_score_pct
  )

function_hit <- read_excel(
  input_xlsx,
  sheet = "FunctionHit",
  .name_repair = "minimal"
)
function_suffix <- ".cds.Function.presence"
function_mag_columns <- names(function_hit)[endsWith(
  names(function_hit),
  function_suffix
)]
if (length(function_mag_columns) != 48L) {
  stop("Expected 48 MAG columns in FunctionHit; observed ",
       length(function_mag_columns), ".")
}

arsenic_features <- c("Arsenate reduction", "Arsenite oxidation")
arsenic_scores <- function_hit %>%
  filter(Category == "As cycling", Function %in% arsenic_features) %>%
  select(Category, Function, all_of(function_mag_columns)) %>%
  pivot_longer(
    all_of(function_mag_columns),
    names_to = "metabolic_column",
    values_to = "function_status"
  ) %>%
  mutate(
    MAG = normalize_metabolic_mag(metabolic_column, function_suffix),
    function_present = str_to_lower(as.character(function_status)) == "present"
  ) %>%
  group_by(Function, MAG) %>%
  summarise(
    feature_component_count = n(),
    components_detected = sum(function_present),
    reconstruction_score_pct = 100 * mean(function_present),
    .groups = "drop"
  ) %>%
  transmute(
    MAG,
    research_category = "Arsenic",
    feature_name = Function,
    source_type = "METABOLIC FunctionHit presence",
    feature_component_count,
    components_detected,
    reconstruction_score_pct
  )

host_scores <- bind_rows(module_scores, arsenic_scores) %>%
  mutate(score_ge_50 = reconstruction_score_pct >= 50)
if (!setequal(unique(host_scores$MAG), unique(mapping$MAG))) {
  stop("METABOLIC workbook MAG names do not match the clade classification.")
}

membership_scores <- host_scores %>%
  inner_join(
    mapping %>% select(MAG, contig_gene, habitat, final_clade,
                       final_clade_name),
    by = "MAG",
    relationship = "many-to-many"
  )

group_summary <- membership_scores %>%
  group_by(research_category, feature_name, source_type, habitat,
           final_clade, final_clade_name) %>%
  summarise(
    n_sequence_memberships = n(),
    n_unique_MAG = n_distinct(MAG),
    mean_reconstruction_score_pct = mean(reconstruction_score_pct),
    prevalence_score_ge_50_pct = 100 * mean(score_ge_50),
    inference_allowed = n_distinct(MAG) > 1,
    .groups = "drop"
  )

category_levels <- c("Carbon", "Nitrogen", "Sulfur", "Arsenic", "Vitamin B12")
feature_order <- bind_rows(
  selected_modules,
  tibble(
    research_category = "Arsenic",
    feature_name = arsenic_features
  )
) %>%
  mutate(research_category = factor(research_category, levels = category_levels)) %>%
  arrange(research_category, feature_name) %>%
  transmute(plot_label = paste(research_category, str_wrap(feature_name, 54),
                               sep = " | ")) %>%
  pull(plot_label)

plot_data <- group_summary %>%
  mutate(
    habitat = factor(habitat, levels = c("IS", "AS", "ES", "NS")),
    final_clade = factor(final_clade, levels = clade_levels),
    plot_label = paste(research_category, str_wrap(feature_name, 54), sep = " | "),
    plot_label = factor(plot_label, levels = rev(unique(feature_order)))
  )

clade_labels <- setNames(
  paste0(
    clade_levels,
    " (",
    as.integer(observed_clade_counts[clade_levels]),
    " sequences)\n",
    str_wrap(clade_names[clade_levels], width = 31)
  ),
  clade_levels
)

caption_text <- paste(
  "KEGG modules are scored by the percentage of METABOLIC module steps detected;",
  "arsenic functions are scored from METABOLIC FunctionHit presence.",
  "Bubble size is the proportion of K08356 sequence memberships whose MAG score is at least 50%;",
  "fill is the mean score. Habitat assignments come from the group-derep input manifest.",
  "IS, incipient seep; AS, active seep; ES, extinct seep; NS, non-seep background.",
  "A MAG carrying two K08356 homologs contributes once to each assigned clade."
) %>% str_wrap(width = 190)

p <- ggplot(
  plot_data,
  aes(
    x = habitat,
    y = plot_label,
    size = prevalence_score_ge_50_pct,
    fill = mean_reconstruction_score_pct
  )
) +
  geom_point(shape = 21, color = "#263238", stroke = 0.28, alpha = 0.96) +
  facet_grid(
    . ~ final_clade,
    labeller = as_labeller(clade_labels),
    drop = FALSE
  ) +
  scale_size_area(
    max_size = 9,
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100),
    name = "MAG prevalence (%)\n(reconstruction score >=50%)"
  ) +
  scale_fill_gradientn(
    colors = c("#F7FCF0", "#7FCDBB", "#2C7FB8", "#253494"),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Mean reconstruction\nscore (%)"
  ) +
  labs(
    title = "Metabolic patterns of K08356-containing MAGs across four habitats",
    subtitle = "METABOLIC-G v4.0; 49 K08356 sequences from 48 unique group-dereplicated MAGs",
    x = "Habitat",
    y = "Research-focused metabolic module or function",
    caption = caption_text
  ) +
  theme_bw(base_family = "serif", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10, color = "#37474F"),
    plot.caption = element_text(size = 8.3, hjust = 0, color = "#37474F"),
    strip.background = element_rect(fill = "#D9E2E7", color = "#455A64"),
    strip.text = element_text(face = "bold", size = 9.4, lineheight = 1.05),
    axis.text.x = element_text(face = "bold", size = 9.5),
    axis.text.y = element_text(size = 7.4, color = "#263238"),
    panel.grid.major = element_line(color = "#E2E8EB", linewidth = 0.28),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box.spacing = unit(3, "mm")
  )

file_stem <- "K08356_MAG49_METABOLIC_four_habitat_four_clade"
ggsave(
  file.path(out_dir, paste0(file_stem, ".png")),
  p,
  width = 17,
  height = 12.5,
  dpi = 360,
  bg = "white"
)
ggsave(
  file.path(out_dir, paste0(file_stem, ".pdf")),
  p,
  width = 17,
  height = 12.5,
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(out_dir, paste0(file_stem, ".svg")),
  p,
  width = 17,
  height = 12.5,
  device = grDevices::svg,
  bg = "white"
)

write.table(
  mapping,
  file.path(out_dir, "MAG49_clade_habitat_mapping.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  mapping %>% select(contig_gene, MAG, habitat_in_classification,
                     habitat, habitat_source_match),
  file.path(out_dir, "habitat_assignment_QC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  host_scores,
  file.path(out_dir, "MAG48_metabolic_feature_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  group_summary,
  file.path(out_dir, "habitat_clade_metabolic_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  plot_data,
  file.path(out_dir, "bubble_plot_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  bind_rows(
    selected_modules %>% mutate(source_type = "KEGG module-step coverage"),
    tibble(
      research_category = "Arsenic",
      feature_name = arsenic_features,
      source_type = "METABOLIC FunctionHit presence"
    )
  ),
  file.path(out_dir, "selected_metabolic_features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("K08356 sequence memberships:", nrow(mapping), "\n")
cat("Unique MAGs:", n_distinct(mapping$MAG), "\n")
cat("Clade counts:", paste(as.integer(observed_clade_counts), collapse = "/"), "\n")
cat("Habitat corrections from manifest:", sum(!mapping$habitat_source_match), "\n")
cat("Selected metabolic features:", n_distinct(host_scores$feature_name), "\n")
cat("Output directory:", out_dir, "\n")

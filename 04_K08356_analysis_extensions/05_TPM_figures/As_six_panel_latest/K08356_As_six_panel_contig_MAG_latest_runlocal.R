suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(grid)
})

# =============================================================================
# K08356 / arsenic-cycle six-panel figure
#
# A: Contig-level arsenic-gene relative abundance by habitat
# C: Contig-level arsenic-gene TPM heatmap (56 samples x arsenic genes)
# E: Contig total arsenic-gene abundance vs four K08356 host-MAG clades
# B: MAG-abundance-weighted arsenic-gene relative abundance by habitat
# D: MAG-abundance-weighted arsenic-gene heatmap (56 samples x arsenic genes)
# F: MAG-weighted total arsenic-gene abundance vs four K08356 host-MAG clades
#
# Statistical replicate = metagenomic sample.
# Candidate proteins = 48; unique K08356-bearing host MAGs = 47.
# MAG TPM is read-recruitment abundance, not the number of recovered MAGs.
# =============================================================================

# ------------------------------- 1. Paths ------------------------------------
outdir <- "/Users/catherine/Downloads/aioA 蛋白序列建树/TPM_combo_figure_20260728/As_six_panel_latest"

contig_xlsx <- "/Users/catherine/Downloads/Contig 热图柱状图 new/all_koid_tpm-final.xlsx"
contig_sheet <- "all_koid_tpm"

gene_list_xlsx <- "/Users/catherine/Downloads/Contig 热图柱状图 new/精简版_CNSAs_contig热图基因清单.xlsx"
gene_list_sheet <- "精简热图基因清单"

mag_tpm_xlsx <- "/Users/catherine/Downloads/MAG-bin-tpm/副本merged_tpm_matrix.xlsx"
mag_tpm_sheet <- "调整顺序后的bin的tpm"

sample_map_tsv <- "/Users/catherine/Downloads/MAG-bin-tpm/derepMAG_abundance_stats/sample_habitat_map.tsv"

strict_qc_tsv <- paste0(
  "/Users/catherine/Downloads/aioA 蛋白序列建树/MAG 补充 idra-tree/",
  "joint_reference_search_48/updated_final_evidence_table_joint_reference_strict_QC.tsv"
)

# Required MAG arsenic-gene table.
# Accepted wide format:
#   MAG_ID | K08356 | K00537 | K03741 | ...
# Accepted long format:
#   MAG_ID | KO | gene_coverage
#
# gene_coverage can be:
#   0/1 presence, gene copy number, or a 0-1 gene-set coverage value.
# Do not use MAG TPM in this table; MAG TPM is read from mag_tpm_xlsx.
mag_as_gene_file <- paste0(
  "/Users/catherine/Downloads/aioA 蛋白序列建树/TPM_combo_figure_20260728/As_six_panel_latest/",
  "host_MAG_arsenic_gene_coverage_47_MAGs.csv"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ----------------------------- 2. Plot settings -------------------------------
habitat_levels <- c("IS", "AS", "ES", "NS")
habitat_cols <- c(
  IS="#66C2A5",
  AS="#F4A261",
  ES="#8DA0CB",
  NS="#E76F51"
)

clade_levels <- c(
  "canonical AioA-associated",
  "unknown/uncertain DMSOR",
  "IdrA-associated without complete DIRM-like synteny",
  "IdrA-associated with DIRM-like synteny"
)

clade_cols <- c(
  "canonical AioA-associated"="#009E73",
  "unknown/uncertain DMSOR"="#8C8C8C",
  "IdrA-associated without complete DIRM-like synteny"="#9E77C8",
  "IdrA-associated with DIRM-like synteny"="#5B2A86"
)

heat_cols <- c("#2166AC", "#F7F7F7", "#F4A582", "#B2182B")

base_theme <- theme_classic(base_family="Times", base_size=9) +
  theme(
    plot.title=element_text(face="bold", size=10.5),
    axis.title=element_text(size=8.5),
    axis.text=element_text(size=7),
    legend.title=element_text(size=7.5),
    legend.text=element_text(size=6.8),
    strip.background=element_rect(fill="grey94", color="grey70", linewidth=0.25),
    strip.text=element_text(face="bold", size=7.2),
    plot.margin=margin(3, 4, 3, 4)
  )

# ------------------------------- 3. Helpers -----------------------------------
fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

sig_code <- function(q) {
  case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ ""
  )
}

parse_depth_mid <- function(x) {
  m <- str_match(x, "[-_](\\d+)-(\\d+)$")
  (as.numeric(m[,2]) + as.numeric(m[,3])) / 2
}

read_table_auto <- function(path, sheet=NULL) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call.=FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    read_excel(path, sheet=sheet)
  } else if (ext == "csv") {
    read_csv(path, show_col_types=FALSE)
  } else {
    read.delim(path, check.names=FALSE)
  }
}

map_clade <- function(x) {
  case_when(
    x == "canonical AioA-associated" ~
      "canonical AioA-associated",
    x == "synteny-supported IdrA/DIRM-like neighborhood-associated lineage" ~
      "IdrA-associated with DIRM-like synteny",
    x == "IdrA-associated phylogenetic lineage" ~
      "IdrA-associated without complete DIRM-like synteny",
    x == "unknown/uncertain DMSOR" ~
      "unknown/uncertain DMSOR",
    x == "unresolved DMSOR-related" ~
      "unknown/uncertain DMSOR",
    TRUE ~ as.character(x)
  )
}

standardize_habitat <- function(x) {
  x <- as.character(x)
  x[x %in% c("BG", "Background", "background")] <- "NS"
  factor(x, levels=habitat_levels)
}

gene_level_tests <- function(df, value_col) {
  value_col <- rlang::ensym(value_col)

  kw <- df %>%
    group_by(gene_label) %>%
    summarise(
      p=suppressWarnings(
        kruskal.test(!!value_col ~ habitat)$p.value
      ),
      .groups="drop"
    ) %>%
    mutate(
      BH_FDR=p.adjust(p, method="BH"),
      significance=sig_code(BH_FDR)
    )

  pairwise <- df %>%
    group_by(gene_label) %>%
    group_modify(~{
      pw <- suppressWarnings(
        pairwise.wilcox.test(
          x=dplyr::pull(.x, !!value_col),
          g=.x$habitat,
          p.adjust.method="BH",
          exact=FALSE
        )
      )

      as.data.frame(as.table(pw$p.value), stringsAsFactors=FALSE) %>%
        filter(!is.na(Freq)) %>%
        transmute(
          habitat_1=as.character(Var1),
          habitat_2=as.character(Var2),
          BH_FDR=as.numeric(Freq)
        )
    }) %>%
    ungroup()

  list(kw=kw, pairwise=pairwise)
}

association_stats <- function(df, x_col, y_col) {
  x_col <- rlang::ensym(x_col)
  y_col <- rlang::ensym(y_col)

  ans <- df %>%
    group_by(MAG_subclass, habitat) %>%
    summarise(
      n=sum(is.finite(!!x_col) & is.finite(!!y_col)),
      rho={
        x <- !!x_col
        y <- !!y_col
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
          NA_real_
        } else {
          unname(suppressWarnings(
            cor.test(x[ok], y[ok], method="spearman", exact=FALSE)$estimate
          ))
        }
      },
      p={
        x <- !!x_col
        y <- !!y_col
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
          NA_real_
        } else {
          suppressWarnings(
            cor.test(x[ok], y[ok], method="spearman", exact=FALSE)$p.value
          )
        }
      },
      .groups="drop"
    ) %>%
    mutate(BH_FDR=p.adjust(p, method="BH"))

  label_df <- ans %>%
    arrange(MAG_subclass, habitat) %>%
    mutate(
      line=paste0(
        habitat, ": rho=",
        ifelse(is.na(rho), "NA", sprintf("%.2f", rho)),
        ", q=", vapply(BH_FDR, fmt_p, character(1))
      )
    ) %>%
    group_by(MAG_subclass) %>%
    summarise(label=paste(line, collapse="\n"), .groups="drop")

  list(table=ans, labels=label_df)
}

# ------------------------------- 4. Metadata ----------------------------------
sample_meta <- read.delim(sample_map_tsv, check.names=FALSE)

if (!all(c("Sample", "Habitat") %in% names(sample_meta))) {
  stop("sample_habitat_map.tsv must contain Sample and Habitat columns.", call.=FALSE)
}

sample_meta <- sample_meta %>%
  transmute(
    sample=as.character(Sample),
    habitat=standardize_habitat(Habitat),
    depth_mid=parse_depth_mid(sample)
  ) %>%
  filter(!is.na(habitat)) %>%
  distinct(sample, .keep_all=TRUE)

# Preferred ordering: habitat, depth layer, sample name.
sample_order <- sample_meta %>%
  arrange(habitat, depth_mid, sample) %>%
  pull(sample)

if (nrow(sample_meta) != 56) {
  stop(
    "Expected 56 metagenomic samples, but sample map contains ",
    nrow(sample_meta), ". Check sample_habitat_map.tsv.",
    call.=FALSE
  )
}

# ------------------------- 5. Arsenic-gene definition -------------------------
gene_key_raw <- read_excel(gene_list_xlsx, sheet=gene_list_sheet)

gene_key <- gene_key_raw %>%
  mutate(
    KO=str_extract(as.character(`KO/注释标识`), "K\\d{5}"),
    gene=as.character(`推荐基因`),
    module_cn=as.character(`模块`),
    keep_flag=as.character(`主图建议`)
  ) %>%
  filter(
    str_detect(keep_flag, "保留"),
    str_detect(module_cn, "砷"),
    !is.na(KO),
    !is.na(gene)
  ) %>%
  transmute(
    KO,
    gene,
    pathway=case_when(
      KO == "K08356" | str_detect(gene, regex("aioA|K08356", TRUE)) ~
        "DMSOR-family K08356",
      str_detect(gene, regex("arrA|arsC", TRUE)) ~
        "As(V) reduction",
      str_detect(gene, regex("arsM", TRUE)) ~
        "Methylation",
      str_detect(gene, regex("acr3|arsB|arsA", TRUE)) ~
        "As transport",
      str_detect(gene, regex("arsR|arsH", TRUE)) ~
        "Regulation / detoxification",
      TRUE ~ "Other arsenic functions"
    )
  ) %>%
  distinct(KO, .keep_all=TRUE)

# Always retain K08356 if it is present in the contig table.
if (!"K08356" %in% gene_key$KO) {
  gene_key <- bind_rows(
    gene_key,
    tibble(
      KO="K08356",
      gene="K08356 (DMSOR-family AioA/IdrA-related)",
      pathway="DMSOR-family K08356"
    )
  )
}

pathway_levels <- c(
  "DMSOR-family K08356",
  "As(V) reduction",
  "Methylation",
  "As transport",
  "Regulation / detoxification",
  "Other arsenic functions"
)

gene_key <- gene_key %>%
  mutate(
    pathway=factor(pathway, levels=pathway_levels),
    gene_label=make.unique(gene)
  ) %>%
  arrange(pathway, gene_label)

gene_levels <- gene_key$gene_label

# ------------------------- 6. Contig-level abundance --------------------------
contig <- read_excel(contig_xlsx, sheet=contig_sheet)
names(contig)[1] <- "KO"
contig <- contig %>%
  mutate(KO=as.character(KO)) %>%
  group_by(KO) %>%
  summarise(
    across(where(is.numeric), ~sum(.x, na.rm=TRUE)),
    .groups="drop"
  )

contig_sample_cols <- intersect(sample_order, names(contig))
if (length(contig_sample_cols) != 56) {
  stop(
    "Contig table matched ", length(contig_sample_cols),
    " of 56 samples. Check sample names.", call.=FALSE
  )
}

missing_contig_KO <- setdiff(gene_key$KO, contig$KO)
if (length(missing_contig_KO) > 0) {
  warning(
    "These arsenic KOs are absent from the contig TPM table and will be zero: ",
    paste(missing_contig_KO, collapse=", ")
  )
}

contig_long <- contig %>%
  filter(KO %in% gene_key$KO) %>%
  select(KO, all_of(contig_sample_cols)) %>%
  pivot_longer(-KO, names_to="sample", values_to="TPM") %>%
  mutate(TPM=replace_na(as.numeric(TPM), 0)) %>%
  complete(
    KO=gene_key$KO,
    sample=sample_order,
    fill=list(TPM=0)
  ) %>%
  left_join(gene_key, by="KO") %>%
  left_join(sample_meta, by="sample") %>%
  group_by(sample) %>%
  mutate(
    total_As_TPM=sum(TPM, na.rm=TRUE),
    RA=if_else(total_As_TPM > 0, TPM / total_As_TPM * 100, 0)
  ) %>%
  ungroup() %>%
  mutate(
    logTPM=log10(TPM+1),
    sample=factor(sample, levels=rev(sample_order)),
    habitat=factor(habitat, levels=habitat_levels),
    pathway=factor(pathway, levels=pathway_levels),
    gene_label=factor(gene_label, levels=gene_levels)
  )

contig_total <- contig_long %>%
  distinct(sample, habitat, total_As_TPM) %>%
  mutate(
    sample=as.character(sample),
    contig_total_As_logTPM=log10(total_As_TPM+1)
  )

contig_tests <- gene_level_tests(contig_long, RA)

contig_bar <- contig_long %>%
  group_by(gene_label, gene, pathway, habitat) %>%
  summarise(
    mean_RA=mean(RA, na.rm=TRUE),
    se_RA=sd(RA, na.rm=TRUE)/sqrt(sum(is.finite(RA))),
    .groups="drop"
  ) %>%
  left_join(
    contig_tests$kw %>%
      select(gene_label, global_BH_FDR=BH_FDR, significance),
    by="gene_label"
  )

contig_sig <- contig_bar %>%
  group_by(gene_label, pathway, significance) %>%
  summarise(
    y=max(mean_RA+replace_na(se_RA, 0), na.rm=TRUE)*1.12+0.15,
    .groups="drop"
  ) %>%
  filter(significance != "")

# -------------------- 7. Unique host-MAG / clade mapping ----------------------
strict <- read.delim(strict_qc_tsv, check.names=FALSE)

required_strict <- c("Sequence_ID", "MAG_ID", "Final_class")
if (!all(required_strict %in% names(strict))) {
  stop(
    "Strict-QC table must contain: ",
    paste(required_strict, collapse=", "), call.=FALSE
  )
}

sequence_map <- strict %>%
  transmute(
    Sequence_ID=as.character(Sequence_ID),
    MAG_ID=as.character(MAG_ID),
    sequence_clade=map_clade(Final_class)
  ) %>%
  filter(sequence_clade %in% clade_levels)

mag_clade_map <- sequence_map %>%
  group_by(MAG_ID) %>%
  summarise(
    n_candidates=n_distinct(Sequence_ID),
    n_clades=n_distinct(sequence_clade),
    MAG_subclass=if_else(
      n_clades == 1,
      first(sequence_clade),
      "mixed K08356 subclasses"
    ),
    candidate_ids=paste(unique(Sequence_ID), collapse=";"),
    candidate_clades=paste(unique(sequence_clade), collapse=";"),
    .groups="drop"
  )

if (n_distinct(sequence_map$Sequence_ID) != 48) {
  warning(
    "Strict-QC mapping contains ",
    n_distinct(sequence_map$Sequence_ID),
    " candidate sequences rather than 48."
  )
}

if (nrow(mag_clade_map) != 47) {
  warning(
    "Unique host-MAG mapping contains ",
    nrow(mag_clade_map),
    " MAGs rather than 47."
  )
}

# ----------------------------- 8. Host-MAG TPM --------------------------------
mag_tpm <- read_excel(mag_tpm_xlsx, sheet=mag_tpm_sheet)
names(mag_tpm)[1] <- "MAG_ID"
mag_tpm <- mag_tpm %>%
  mutate(MAG_ID=as.character(MAG_ID)) %>%
  filter(!is.na(MAG_ID))

mag_sample_cols <- intersect(sample_order, names(mag_tpm))
if (length(mag_sample_cols) != 56) {
  stop(
    "MAG TPM table matched ", length(mag_sample_cols),
    " of 56 samples. Check sample names.", call.=FALSE
  )
}

mag_host_long <- mag_tpm %>%
  filter(MAG_ID %in% mag_clade_map$MAG_ID) %>%
  select(MAG_ID, all_of(mag_sample_cols)) %>%
  pivot_longer(-MAG_ID, names_to="sample", values_to="MAG_TPM") %>%
  mutate(MAG_TPM=replace_na(as.numeric(MAG_TPM), 0)) %>%
  left_join(mag_clade_map, by="MAG_ID") %>%
  left_join(sample_meta, by="sample")

# The mixed-clade MAG is kept once for the MAG functional weighting,
# but excluded from four-clade abundance sums to prevent double counting.
clade_sample <- mag_host_long %>%
  filter(MAG_subclass %in% clade_levels) %>%
  group_by(MAG_subclass, sample, habitat) %>%
  summarise(
    clade_MAG_TPM=sum(MAG_TPM, na.rm=TRUE),
    n_host_MAG=n_distinct(MAG_ID),
    .groups="drop"
  ) %>%
  complete(
    MAG_subclass=clade_levels,
    nesting(sample, habitat),
    fill=list(clade_MAG_TPM=0, n_host_MAG=0)
  ) %>%
  mutate(
    MAG_subclass=factor(MAG_subclass, levels=clade_levels),
    habitat=factor(habitat, levels=habitat_levels),
    clade_MAG_logTPM=log10(clade_MAG_TPM+1)
  )

# ----------------- 9. MAG arsenic-gene coverage / presence --------------------
if (!file.exists(mag_as_gene_file)) {
  template <- crossing(
    MAG_ID=mag_clade_map$MAG_ID,
    KO=gene_key$KO
  ) %>%
    left_join(gene_key %>% select(KO, gene, pathway), by="KO") %>%
    mutate(gene_coverage=NA_real_)

  template_file <- file.path(
    outdir, "TEMPLATE_host_MAG_arsenic_gene_coverage_long.csv"
  )
  write_csv(template, template_file)

  stop(
    "MAG arsenic-gene table is missing: ", mag_as_gene_file, "\n",
    "A template has been written to: ", template_file, "\n",
    "Fill gene_coverage and update mag_as_gene_file.",
    call.=FALSE
  )
}

mag_gene_raw <- read_table_auto(mag_as_gene_file)

mag_id_col <- intersect(
  c("MAG_ID", "mag_id", "Bin", "bin", "Genome", "genome"),
  names(mag_gene_raw)
)[1]

if (is.na(mag_id_col)) {
  stop("MAG arsenic-gene table has no recognizable MAG_ID column.", call.=FALSE)
}

names(mag_gene_raw)[names(mag_gene_raw) == mag_id_col] <- "MAG_ID"
mag_gene_raw$MAG_ID <- as.character(mag_gene_raw$MAG_ID)

if ("KO" %in% names(mag_gene_raw)) {
  coverage_col <- intersect(
    c("gene_coverage", "coverage", "copy_number", "copies", "presence", "value"),
    names(mag_gene_raw)
  )[1]

  if (is.na(coverage_col)) {
    stop(
      "Long MAG gene table must contain gene_coverage, coverage, ",
      "copy_number, copies, presence, or value.", call.=FALSE
    )
  }

  mag_gene_long <- mag_gene_raw %>%
    transmute(
      MAG_ID,
      KO=as.character(KO),
      gene_coverage=as.numeric(.data[[coverage_col]])
    )
} else {
  ko_cols <- intersect(gene_key$KO, names(mag_gene_raw))
  if (length(ko_cols) == 0) {
    stop(
      "Wide MAG gene table has no columns matching arsenic-gene KO IDs.",
      call.=FALSE
    )
  }

  mag_gene_long <- mag_gene_raw %>%
    select(MAG_ID, all_of(ko_cols)) %>%
    pivot_longer(-MAG_ID, names_to="KO", values_to="gene_coverage") %>%
    mutate(gene_coverage=as.numeric(gene_coverage))
}

mag_gene_long <- mag_gene_long %>%
  filter(
    MAG_ID %in% mag_clade_map$MAG_ID,
    KO %in% gene_key$KO
  ) %>%
  complete(
    MAG_ID=mag_clade_map$MAG_ID,
    KO=gene_key$KO,
    fill=list(gene_coverage=0)
  ) %>%
  mutate(gene_coverage=replace_na(gene_coverage, 0)) %>%
  left_join(gene_key, by="KO")

# Per sample and arsenic gene:
# weighted TPM = sum across 47 unique hosts (MAG TPM x gene coverage).
mag_weighted_long <- mag_host_long %>%
  select(MAG_ID, sample, habitat, MAG_TPM) %>%
  inner_join(
    mag_gene_long %>%
      select(MAG_ID, KO, gene, gene_label, pathway, gene_coverage),
    by="MAG_ID"
  ) %>%
  mutate(weighted_TPM=MAG_TPM*gene_coverage) %>%
  group_by(sample, habitat, KO, gene, gene_label, pathway) %>%
  summarise(
    weighted_TPM=sum(weighted_TPM, na.rm=TRUE),
    .groups="drop"
  ) %>%
  complete(
    nesting(sample, habitat),
    KO=gene_key$KO,
    fill=list(weighted_TPM=0)
  ) %>%
  left_join(
    gene_key %>% select(KO, gene_key_gene=gene,
                        gene_key_label=gene_label,
                        gene_key_pathway=pathway),
    by="KO"
  ) %>%
  mutate(
    gene=coalesce(gene, gene_key_gene),
    gene_label=coalesce(as.character(gene_label), gene_key_label),
    pathway=coalesce(as.character(pathway), as.character(gene_key_pathway))
  ) %>%
  select(-gene_key_gene, -gene_key_label, -gene_key_pathway) %>%
  group_by(sample) %>%
  mutate(
    total_weighted_As_TPM=sum(weighted_TPM, na.rm=TRUE),
    RA=if_else(
      total_weighted_As_TPM > 0,
      weighted_TPM/total_weighted_As_TPM*100,
      0
    )
  ) %>%
  ungroup() %>%
  mutate(
    logWeightedTPM=log10(weighted_TPM+1),
    sample=factor(sample, levels=rev(sample_order)),
    habitat=factor(habitat, levels=habitat_levels),
    pathway=factor(pathway, levels=pathway_levels),
    gene_label=factor(gene_label, levels=gene_levels)
  )

mag_weighted_total <- mag_weighted_long %>%
  distinct(sample, habitat, total_weighted_As_TPM) %>%
  mutate(
    sample=as.character(sample),
    MAG_weighted_total_As_logTPM=log10(total_weighted_As_TPM+1)
  )

mag_tests <- gene_level_tests(mag_weighted_long, RA)

mag_bar <- mag_weighted_long %>%
  group_by(gene_label, gene, pathway, habitat) %>%
  summarise(
    mean_RA=mean(RA, na.rm=TRUE),
    se_RA=sd(RA, na.rm=TRUE)/sqrt(sum(is.finite(RA))),
    .groups="drop"
  ) %>%
  left_join(
    mag_tests$kw %>%
      select(gene_label, global_BH_FDR=BH_FDR, significance),
    by="gene_label"
  )

mag_sig <- mag_bar %>%
  group_by(gene_label, pathway, significance) %>%
  summarise(
    y=max(mean_RA+replace_na(se_RA, 0), na.rm=TRUE)*1.12+0.15,
    .groups="drop"
  ) %>%
  filter(significance != "")

# ------------------------ 10. E/F association data ----------------------------
contig_assoc <- clade_sample %>%
  mutate(sample=as.character(sample)) %>%
  left_join(contig_total, by=c("sample", "habitat"))

mag_assoc <- clade_sample %>%
  mutate(sample=as.character(sample)) %>%
  left_join(mag_weighted_total, by=c("sample", "habitat"))

contig_assoc_stats <- association_stats(
  contig_assoc,
  contig_total_As_logTPM,
  clade_MAG_logTPM
)

mag_assoc_stats <- association_stats(
  mag_assoc,
  MAG_weighted_total_As_logTPM,
  clade_MAG_logTPM
)

# ------------------------------- 11. Panels -----------------------------------
pA <- ggplot(
  contig_bar,
  aes(gene_label, mean_RA, fill=habitat)
) +
  geom_col(
    position=position_dodge(width=0.82),
    width=0.76,
    color="grey30",
    linewidth=0.15
  ) +
  geom_errorbar(
    aes(
      ymin=pmax(0, mean_RA-se_RA),
      ymax=mean_RA+se_RA
    ),
    position=position_dodge(width=0.82),
    width=0.18,
    linewidth=0.25
  ) +
  geom_text(
    data=contig_sig,
    aes(gene_label, y, label=significance),
    inherit.aes=FALSE,
    size=3.4,
    family="Times",
    fontface="bold"
  ) +
  facet_grid(. ~ pathway, scales="free_x", space="free_x") +
  scale_fill_manual(values=habitat_cols, drop=FALSE) +
  scale_y_continuous(
    labels=function(x) paste0(x, "%"),
    expand=expansion(mult=c(0, 0.14))
  ) +
  labs(
    title="A  Contig-level arsenic-gene relative abundance",
    x=NULL,
    y="Relative abundance (%)",
    fill="Habitat"
  ) +
  base_theme +
  theme(
    axis.text.x=element_text(angle=55, hjust=1, size=6.1),
    legend.position="top",
    panel.spacing.x=unit(0.8, "mm")
  )

pC <- ggplot(
  contig_long,
  aes(gene_label, sample, fill=logTPM)
) +
  geom_tile(color="white", linewidth=0.06) +
  facet_grid(
    habitat ~ pathway,
    scales="free",
    space="free"
  ) +
  scale_fill_gradientn(
    colors=heat_cols,
    name="log10(TPM+1)"
  ) +
  labs(
    title="C  Contig-level arsenic-gene abundance across 56 samples",
    x=NULL,
    y=NULL
  ) +
  base_theme +
  theme(
    axis.text.x=element_text(angle=55, hjust=1, size=5.8),
    axis.text.y=element_text(size=4.4),
    strip.text.y=element_text(angle=0),
    panel.spacing=unit(0.7, "mm"),
    legend.position="right"
  )

pB <- ggplot(
  mag_bar,
  aes(gene_label, mean_RA, fill=habitat)
) +
  geom_col(
    position=position_dodge(width=0.82),
    width=0.76,
    color="grey30",
    linewidth=0.15
  ) +
  geom_errorbar(
    aes(
      ymin=pmax(0, mean_RA-se_RA),
      ymax=mean_RA+se_RA
    ),
    position=position_dodge(width=0.82),
    width=0.18,
    linewidth=0.25
  ) +
  geom_text(
    data=mag_sig,
    aes(gene_label, y, label=significance),
    inherit.aes=FALSE,
    size=3.4,
    family="Times",
    fontface="bold"
  ) +
  facet_grid(. ~ pathway, scales="free_x", space="free_x") +
  scale_fill_manual(values=habitat_cols, drop=FALSE) +
  scale_y_continuous(
    labels=function(x) paste0(x, "%"),
    expand=expansion(mult=c(0, 0.14))
  ) +
  labs(
    title="B  MAG-weighted arsenic-gene relative abundance",
    subtitle="Weighted TPM = host MAG TPM x gene coverage, summed over 47 unique host MAGs",
    x=NULL,
    y="Relative abundance (%)",
    fill="Habitat"
  ) +
  base_theme +
  theme(
    plot.subtitle=element_text(size=7),
    axis.text.x=element_text(angle=55, hjust=1, size=6.1),
    legend.position="top",
    panel.spacing.x=unit(0.8, "mm")
  )

pD <- ggplot(
  mag_weighted_long,
  aes(gene_label, sample, fill=logWeightedTPM)
) +
  geom_tile(color="white", linewidth=0.06) +
  facet_grid(
    habitat ~ pathway,
    scales="free",
    space="free"
  ) +
  scale_fill_gradientn(
    colors=heat_cols,
    name="log10(weighted\nTPM+1)"
  ) +
  labs(
    title="D  MAG-weighted arsenic-gene abundance across 56 samples",
    x=NULL,
    y=NULL
  ) +
  base_theme +
  theme(
    axis.text.x=element_text(angle=55, hjust=1, size=5.8),
    axis.text.y=element_text(size=4.4),
    strip.text.y=element_text(angle=0),
    panel.spacing=unit(0.7, "mm"),
    legend.position="right"
  )

pE <- ggplot(
  contig_assoc,
  aes(
    contig_total_As_logTPM,
    clade_MAG_logTPM,
    color=habitat,
    fill=habitat
  )
) +
  geom_point(size=1.25, alpha=0.75) +
  geom_smooth(
    aes(group=habitat),
    method="lm",
    formula=y ~ x,
    se=TRUE,
    linewidth=0.55,
    alpha=0.12
  ) +
  facet_wrap(
    ~MAG_subclass,
    ncol=2,
    labeller=label_wrap_gen(width=27)
  ) +
  geom_text(
    data=contig_assoc_stats$labels,
    aes(x=-Inf, y=Inf, label=label),
    inherit.aes=FALSE,
    hjust=-0.03,
    vjust=1.05,
    size=2.05,
    family="Times",
    lineheight=0.92
  ) +
  scale_color_manual(values=habitat_cols, drop=FALSE) +
  scale_fill_manual(values=habitat_cols, drop=FALSE) +
  labs(
    title="E  Contig arsenic-gene abundance versus four K08356 host-MAG clades",
    x="Total contig arsenic-gene log10(TPM + 1)",
    y="Clade-summed host MAG log10(TPM + 1)",
    color="Habitat",
    fill="Habitat"
  ) +
  base_theme +
  theme(
    legend.position="top",
    strip.text=element_text(size=6.6)
  )

pF <- ggplot(
  mag_assoc,
  aes(
    MAG_weighted_total_As_logTPM,
    clade_MAG_logTPM,
    color=habitat,
    fill=habitat
  )
) +
  geom_point(size=1.25, alpha=0.75) +
  geom_smooth(
    aes(group=habitat),
    method="lm",
    formula=y ~ x,
    se=TRUE,
    linewidth=0.55,
    alpha=0.12
  ) +
  facet_wrap(
    ~MAG_subclass,
    ncol=2,
    labeller=label_wrap_gen(width=27)
  ) +
  geom_text(
    data=mag_assoc_stats$labels,
    aes(x=-Inf, y=Inf, label=label),
    inherit.aes=FALSE,
    hjust=-0.03,
    vjust=1.05,
    size=2.05,
    family="Times",
    lineheight=0.92
  ) +
  scale_color_manual(values=habitat_cols, drop=FALSE) +
  scale_fill_manual(values=habitat_cols, drop=FALSE) +
  labs(
    title="F  MAG-weighted arsenic-gene abundance versus four K08356 host-MAG clades",
    x="MAG-weighted total arsenic-gene log10(TPM + 1)",
    y="Clade-summed host MAG log10(TPM + 1)",
    color="Habitat",
    fill="Habitat"
  ) +
  base_theme +
  theme(
    legend.position="top",
    strip.text=element_text(size=6.6)
  )

# Reference-style layout: A/C above B/D on the left; E/F on the right.
left_top <- plot_grid(
  pA, pC,
  ncol=1,
  rel_heights=c(0.72, 1.12),
  align="v",
  axis="lr"
)

left_bottom <- plot_grid(
  pB, pD,
  ncol=1,
  rel_heights=c(0.78, 1.12),
  align="v",
  axis="lr"
)

left <- plot_grid(
  left_top,
  left_bottom,
  ncol=1,
  rel_heights=c(1, 1.03)
)

right <- plot_grid(
  pE,
  pF,
  ncol=1,
  rel_heights=c(1, 1),
  align="v",
  axis="lr"
)

body <- plot_grid(
  left,
  right,
  ncol=2,
  rel_widths=c(1.78, 1),
  align="h"
)

title_panel <- ggdraw() +
  draw_label(
    "Contig- and MAG-resolved arsenic functional potential associated with four K08356 evidence subclasses",
    x=0,
    y=0.72,
    hjust=0,
    fontfamily="Times",
    fontface="bold",
    size=14
  ) +
  draw_label(
    "A/C use contig KO-level TPM; B/D use host-MAG abundance-weighted gene coverage. E/F show sample-level co-variation and do not establish enzymatic function.",
    x=0,
    y=0.22,
    hjust=0,
    fontfamily="Times",
    size=8.3
  )

figure <- plot_grid(
  title_panel,
  body,
  ncol=1,
  rel_heights=c(0.065, 1)
)

# ------------------------------ 12. Outputs -----------------------------------
prefix <- file.path(outdir, "K08356_As_six_panel_contig_MAG_latest")

ggsave(
  paste0(prefix, ".pdf"),
  figure,
  width=16.8,
  height=14.2,
  units="in",
  device=cairo_pdf
)

ggsave(
  paste0(prefix, ".svg"),
  figure,
  width=16.8,
  height=14.2,
  units="in",
  device=svg
)

ggsave(
  paste0(prefix, ".png"),
  figure,
  width=16.8,
  height=14.2,
  units="in",
  dpi=600
)

write_csv(gene_key, file.path(outdir, "arsenic_gene_key_used.csv"))
write_csv(contig_long, file.path(outdir, "A_C_contig_arsenic_gene_long.csv"))
write_csv(contig_bar, file.path(outdir, "A_contig_habitat_RA_summary.csv"))
write_csv(contig_tests$kw, file.path(outdir, "A_contig_gene_Kruskal_Wallis_BH.csv"))
write_csv(contig_tests$pairwise, file.path(outdir, "A_contig_gene_pairwise_Wilcoxon_BH.csv"))

write_csv(mag_clade_map, file.path(outdir, "K08356_unique_MAG_clade_map.csv"))
write_csv(mag_gene_long, file.path(outdir, "MAG_arsenic_gene_coverage_used.csv"))
write_csv(mag_weighted_long, file.path(outdir, "B_D_MAG_weighted_arsenic_gene_long.csv"))
write_csv(mag_bar, file.path(outdir, "B_MAG_weighted_habitat_RA_summary.csv"))
write_csv(mag_tests$kw, file.path(outdir, "B_MAG_gene_Kruskal_Wallis_BH.csv"))
write_csv(mag_tests$pairwise, file.path(outdir, "B_MAG_gene_pairwise_Wilcoxon_BH.csv"))

write_csv(clade_sample, file.path(outdir, "four_clade_host_MAG_TPM_by_sample.csv"))
write_csv(contig_assoc, file.path(outdir, "E_contig_As_vs_four_clades_data.csv"))
write_csv(contig_assoc_stats$table, file.path(outdir, "E_contig_As_vs_four_clades_Spearman_BH.csv"))
write_csv(mag_assoc, file.path(outdir, "F_MAG_weighted_As_vs_four_clades_data.csv"))
write_csv(mag_assoc_stats$table, file.path(outdir, "F_MAG_weighted_As_vs_four_clades_Spearman_BH.csv"))

qc <- c(
  paste("metagenomic_samples", nrow(sample_meta), sep="\t"),
  paste(
    "habitat_counts",
    paste(
      names(table(sample_meta$habitat)),
      as.integer(table(sample_meta$habitat)),
      sep="=",
      collapse=";"
    ),
    sep="\t"
  ),
  paste("arsenic_genes", nrow(gene_key), sep="\t"),
  paste("K08356_candidate_sequences", n_distinct(sequence_map$Sequence_ID), sep="\t"),
  paste("unique_K08356_host_MAGs", nrow(mag_clade_map), sep="\t"),
  paste(
    "mixed_clade_MAGs_excluded_from_clade_sums",
    paste(
      mag_clade_map$MAG_ID[
        mag_clade_map$MAG_subclass == "mixed K08356 subclasses"
      ],
      collapse=";"
    ),
    sep="\t"
  ),
  paste(
    "MAG_weighting_formula",
    "sum(unique host MAG TPM x arsenic gene coverage)",
    sep="\t"
  ),
  paste("figure_pdf", paste0(prefix, ".pdf"), sep="\t"),
  paste("figure_svg", paste0(prefix, ".svg"), sep="\t"),
  paste("figure_png", paste0(prefix, ".png"), sep="\t")
)

writeLines(qc, file.path(outdir, "K08356_As_six_panel_QC.txt"))
cat(paste(qc, collapse="\n"), "\n")

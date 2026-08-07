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

outdir <- "/Users/catherine/Downloads/aioA 蛋白序列建树/TPM_combo_figure_20260728/As_original_layout_corrected_latest"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

contig_xlsx <- "/Users/catherine/Downloads/Contig 热图柱状图 new/all_koid_tpm-final.xlsx"
gene_list_xlsx <- "/Users/catherine/Downloads/Contig 热图柱状图 new/精简版_CNSAs_contig热图基因清单.xlsx"
mag_tpm_xlsx <- "/Users/catherine/Downloads/MAG-bin-tpm/副本merged_tpm_matrix.xlsx"
sample_map_tsv <- "/Users/catherine/Downloads/MAG-bin-tpm/derepMAG_abundance_stats/sample_habitat_map.tsv"
strict_qc_tsv <- paste0(
  "/Users/catherine/Downloads/aioA 蛋白序列建树/MAG 补充 idra-tree/",
  "joint_reference_search_48/updated_final_evidence_table_joint_reference_strict_QC.tsv"
)
kofam_xlsx <- "/Users/catherine/Downloads/aioA 蛋白序列建树/冷泉-砷-热图/所有 MAG 携带哪些功能基因all_bin_kofamscan_filtered_split.20260121200414199.xlsx"
mag_as_gene_file <- file.path(outdir, "host_MAG_arsenic_gene_coverage_47_MAGs.csv")

habitat_levels <- c("IS", "AS", "ES", "NS")
habitat_cols <- c(IS="#66C2A5", AS="#F4A261", ES="#8DA0CB", NS="#E76F51")

clade_levels <- c(
  "canonical AioA-associated",
  "unknown/uncertain DMSOR",
  "IdrA-associated without complete DIRM-like synteny",
  "IdrA-associated with DIRM-like synteny"
)
clade_cols <- c(
  "canonical AioA-associated"="#F4A3C4",
  "unknown/uncertain DMSOR"="#4DAF4A",
  "IdrA-associated without complete DIRM-like synteny"="#377EB8",
  "IdrA-associated with DIRM-like synteny"="#8E44AD"
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

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
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

standardize_habitat <- function(x) {
  x <- as.character(x)
  x[x %in% c("BG", "Background", "background")] <- "NS"
  factor(x, levels=habitat_levels)
}

map_clade <- function(x) {
  case_when(
    x == "canonical AioA-associated" ~ "canonical AioA-associated",
    x == "synteny-supported IdrA/DIRM-like neighborhood-associated lineage" ~
      "IdrA-associated with DIRM-like synteny",
    x == "IdrA-associated phylogenetic lineage" ~
      "IdrA-associated without complete DIRM-like synteny",
    x == "unknown/uncertain DMSOR" ~ "unknown/uncertain DMSOR",
    TRUE ~ as.character(x)
  )
}

kw_bh_by_gene <- function(df, value_col) {
  value_col <- rlang::ensym(value_col)
  df %>%
    group_by(gene_label) %>%
    summarise(p=suppressWarnings(kruskal.test(!!value_col ~ habitat)$p.value),
              .groups="drop") %>%
    mutate(BH_FDR=p.adjust(p, method="BH"), significance=sig_code(BH_FDR))
}

is_pairwise_stars <- function(df, value_col) {
  value_col <- rlang::ensym(value_col)
  df %>%
    group_by(gene_label) %>%
    group_modify(~{
      out <- lapply(c("AS", "ES", "NS"), function(h) {
        x <- dplyr::pull(filter(.x, habitat == "IS"), !!value_col)
        y <- dplyr::pull(filter(.x, habitat == h), !!value_col)
        p <- if (sum(is.finite(x)) >= 2 && sum(is.finite(y)) >= 2) {
          suppressWarnings(wilcox.test(x, y, exact=FALSE)$p.value)
        } else {
          NA_real_
        }
        tibble(habitat=h, p=p)
      })
      bind_rows(out)
    }) %>%
    ungroup() %>%
    mutate(BH_FDR=p.adjust(p, method="BH"), significance=sig_code(BH_FDR)) %>%
    filter(significance != "")
}

sample_meta <- read.delim(sample_map_tsv, check.names=FALSE) %>%
  transmute(
    sample=as.character(Sample),
    habitat=standardize_habitat(Habitat),
    depth_mid=parse_depth_mid(sample)
  ) %>%
  filter(!is.na(habitat)) %>%
  distinct(sample, .keep_all=TRUE)

if (nrow(sample_meta) != 56) {
  stop("Expected 56 metagenomic samples, found ", nrow(sample_meta), call.=FALSE)
}

sample_order <- sample_meta %>% arrange(habitat, depth_mid, sample) %>% pull(sample)

gene_key_raw <- read_excel(gene_list_xlsx, sheet="精简热图基因清单")
gene_key <- gene_key_raw %>%
  mutate(
    KO=str_extract(as.character(`KO/注释标识`), "K\\d{5}"),
    gene=as.character(`推荐基因`),
    module_cn=as.character(`模块`),
    keep_flag=as.character(`主图建议`)
  ) %>%
  filter(str_detect(keep_flag, "保留"),
         str_detect(module_cn, "砷"),
         !is.na(KO),
         !is.na(gene)) %>%
  transmute(
    KO,
    gene,
    pathway=case_when(
      KO == "K08356" | str_detect(gene, regex("aioA|K08356|相关序列", TRUE)) ~
        "DMSOR-family K08356",
      str_detect(gene, regex("arrA|arsC", TRUE)) ~ "As(V) reduction",
      str_detect(gene, regex("arsM", TRUE)) ~ "Methylation",
      str_detect(gene, regex("acr3|arsB|arsA", TRUE)) ~ "As transport",
      str_detect(gene, regex("arsR|arsH", TRUE)) ~ "Regulation / detoxification",
      TRUE ~ "Other arsenic functions"
    )
  ) %>%
  distinct(KO, .keep_all=TRUE)

if (!"K08356" %in% gene_key$KO) {
  gene_key <- bind_rows(gene_key, tibble(
    KO="K08356",
    gene="K08356 (DMSOR-family AioA/IdrA-related)",
    pathway="DMSOR-family K08356"
  ))
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
  mutate(pathway=factor(pathway, levels=pathway_levels),
         gene_label=make.unique(gene)) %>%
  arrange(pathway, gene_label)
gene_levels <- gene_key$gene_label

contig <- read_excel(contig_xlsx, sheet="all_koid_tpm")
names(contig)[1] <- "KO"
contig <- contig %>%
  mutate(KO=as.character(KO)) %>%
  group_by(KO) %>%
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE)), .groups="drop")
contig_sample_cols <- intersect(sample_order, names(contig))
if (length(contig_sample_cols) != 56) {
  stop("Contig table matched ", length(contig_sample_cols), " of 56 samples.", call.=FALSE)
}

contig_long <- contig %>%
  filter(KO %in% gene_key$KO) %>%
  select(KO, all_of(contig_sample_cols)) %>%
  pivot_longer(-KO, names_to="sample", values_to="TPM") %>%
  mutate(TPM=replace_na(as.numeric(TPM), 0)) %>%
  complete(KO=gene_key$KO, sample=sample_order, fill=list(TPM=0)) %>%
  left_join(gene_key, by="KO") %>%
  left_join(sample_meta, by="sample") %>%
  group_by(sample) %>%
  mutate(total_As_TPM=sum(TPM, na.rm=TRUE),
         RA=if_else(total_As_TPM > 0, TPM/total_As_TPM*100, 0)) %>%
  ungroup() %>%
  mutate(
    logTPM=log10(TPM+1),
    sample=factor(sample, levels=rev(sample_order)),
    habitat=factor(habitat, levels=habitat_levels),
    pathway=factor(pathway, levels=pathway_levels),
    gene_label=factor(gene_label, levels=gene_levels)
  )

contig_tests <- kw_bh_by_gene(contig_long, RA)
contig_pair_stars <- is_pairwise_stars(contig_long, RA)
contig_bar <- contig_long %>%
  group_by(gene_label, gene, pathway, habitat) %>%
  summarise(mean_RA=mean(RA, na.rm=TRUE),
            se_RA=sd(RA, na.rm=TRUE)/sqrt(sum(is.finite(RA))),
            .groups="drop") %>%
  left_join(contig_tests %>% select(gene_label, global_BH_FDR=BH_FDR),
            by="gene_label")
contig_pair_annot <- contig_bar %>%
  filter(as.character(habitat) != "IS") %>%
  left_join(contig_pair_stars, by=c("gene_label", "habitat")) %>%
  filter(significance != "") %>%
  mutate(y=mean_RA + replace_na(se_RA, 0) + 2.0)

strict <- read.delim(strict_qc_tsv, check.names=FALSE)
sequence_map <- strict %>%
  transmute(Sequence_ID=as.character(Sequence_ID),
            MAG_ID=as.character(MAG_ID),
            sequence_clade=map_clade(Final_class)) %>%
  filter(sequence_clade %in% clade_levels)
mag_clade_map <- sequence_map %>%
  group_by(MAG_ID) %>%
  summarise(
    n_candidates=n_distinct(Sequence_ID),
    n_clades=n_distinct(sequence_clade),
    MAG_subclass=if_else(n_clades == 1, first(sequence_clade), "mixed K08356 subclasses"),
    candidate_ids=paste(unique(Sequence_ID), collapse=";"),
    candidate_clades=paste(unique(sequence_clade), collapse=";"),
    .groups="drop"
  )

mag_tpm <- read_excel(mag_tpm_xlsx, sheet="调整顺序后的bin的tpm")
names(mag_tpm)[1] <- "MAG_ID"
mag_tpm <- mag_tpm %>% mutate(MAG_ID=as.character(MAG_ID)) %>% filter(!is.na(MAG_ID))
mag_sample_cols <- intersect(sample_order, names(mag_tpm))
if (length(mag_sample_cols) != 56) {
  stop("MAG TPM table matched ", length(mag_sample_cols), " of 56 samples.", call.=FALSE)
}

mag_host_long <- mag_tpm %>%
  filter(MAG_ID %in% mag_clade_map$MAG_ID) %>%
  select(MAG_ID, all_of(mag_sample_cols)) %>%
  pivot_longer(-MAG_ID, names_to="sample", values_to="MAG_TPM") %>%
  mutate(MAG_TPM=replace_na(as.numeric(MAG_TPM), 0)) %>%
  left_join(mag_clade_map, by="MAG_ID") %>%
  left_join(sample_meta, by="sample")

clade_sample <- mag_host_long %>%
  filter(MAG_subclass %in% clade_levels) %>%
  group_by(MAG_subclass, sample, habitat) %>%
  summarise(clade_MAG_TPM=sum(MAG_TPM, na.rm=TRUE),
            n_host_MAG=n_distinct(MAG_ID),
            .groups="drop") %>%
  complete(MAG_subclass=clade_levels, nesting(sample, habitat),
           fill=list(clade_MAG_TPM=0, n_host_MAG=0)) %>%
  mutate(MAG_subclass=factor(MAG_subclass, levels=clade_levels),
         habitat=factor(habitat, levels=habitat_levels),
         clade_MAG_logTPM=log10(clade_MAG_TPM+1))

if (!file.exists(mag_as_gene_file)) {
  kofam <- read_excel(kofam_xlsx, sheet=1, col_names=FALSE)
  names(kofam)[1:4] <- c("MAG_ID", "protein_id", "KO", "description")
  mag_gene_long_for_write <- crossing(MAG_ID=mag_clade_map$MAG_ID, KO=gene_key$KO) %>%
    left_join(gene_key %>% select(KO, gene, pathway), by="KO") %>%
    left_join(
      kofam %>%
        mutate(MAG_ID=as.character(MAG_ID), KO=as.character(KO)) %>%
        filter(MAG_ID %in% mag_clade_map$MAG_ID, KO %in% gene_key$KO) %>%
        count(MAG_ID, KO, name="gene_coverage"),
      by=c("MAG_ID", "KO")
    ) %>%
    mutate(gene_coverage=replace_na(gene_coverage, 0))
  write_csv(mag_gene_long_for_write, mag_as_gene_file)
}

mag_gene_long <- read_csv(mag_as_gene_file, show_col_types=FALSE) %>%
  transmute(MAG_ID=as.character(MAG_ID),
            KO=as.character(KO),
            gene_coverage=as.numeric(gene_coverage)) %>%
  filter(MAG_ID %in% mag_clade_map$MAG_ID, KO %in% gene_key$KO) %>%
  complete(MAG_ID=mag_clade_map$MAG_ID, KO=gene_key$KO, fill=list(gene_coverage=0)) %>%
  left_join(gene_key, by="KO")

mag_weighted_long <- mag_host_long %>%
  select(MAG_ID, sample, habitat, MAG_TPM) %>%
  inner_join(mag_gene_long %>% select(MAG_ID, KO, gene, gene_label, pathway, gene_coverage),
             by="MAG_ID", relationship="many-to-many") %>%
  mutate(weighted_TPM=MAG_TPM*gene_coverage) %>%
  group_by(sample, habitat, KO, gene, gene_label, pathway) %>%
  summarise(weighted_TPM=sum(weighted_TPM, na.rm=TRUE), .groups="drop") %>%
  complete(nesting(sample, habitat), KO=gene_key$KO, fill=list(weighted_TPM=0)) %>%
  left_join(gene_key %>% select(KO, key_gene=gene, key_label=gene_label, key_pathway=pathway),
            by="KO") %>%
  mutate(gene=coalesce(gene, key_gene),
         gene_label=coalesce(as.character(gene_label), key_label),
         pathway=coalesce(as.character(pathway), as.character(key_pathway))) %>%
  select(-key_gene, -key_label, -key_pathway) %>%
  group_by(sample) %>%
  mutate(total_weighted_As_TPM=sum(weighted_TPM, na.rm=TRUE),
         RA=if_else(total_weighted_As_TPM > 0, weighted_TPM/total_weighted_As_TPM*100, 0)) %>%
  ungroup() %>%
  mutate(logWeightedTPM=log10(weighted_TPM+1),
         sample=factor(sample, levels=rev(sample_order)),
         habitat=factor(habitat, levels=habitat_levels),
         pathway=factor(pathway, levels=pathway_levels),
         gene_label=factor(gene_label, levels=gene_levels))

mag_tests <- kw_bh_by_gene(mag_weighted_long, RA)
mag_pair_stars <- is_pairwise_stars(mag_weighted_long, RA)
mag_bar <- mag_weighted_long %>%
  group_by(gene_label, gene, pathway, habitat) %>%
  summarise(mean_RA=mean(RA, na.rm=TRUE),
            se_RA=sd(RA, na.rm=TRUE)/sqrt(sum(is.finite(RA))),
            .groups="drop") %>%
  left_join(mag_tests %>% select(gene_label, global_BH_FDR=BH_FDR),
            by="gene_label")
mag_pair_annot <- mag_bar %>%
  filter(as.character(habitat) != "IS") %>%
  left_join(mag_pair_stars, by=c("gene_label", "habitat")) %>%
  filter(significance != "") %>%
  mutate(y=mean_RA + replace_na(se_RA, 0) + 2.0)

non_k08356_gene <- gene_key %>% filter(KO != "K08356") %>% pull(KO)
clade_gene_weighted <- mag_host_long %>%
  filter(MAG_subclass %in% clade_levels) %>%
  select(MAG_ID, MAG_subclass, sample, habitat, MAG_TPM) %>%
  inner_join(mag_gene_long %>% filter(KO %in% non_k08356_gene) %>%
               select(MAG_ID, KO, gene_coverage),
             by="MAG_ID", relationship="many-to-many") %>%
  mutate(weighted_TPM=MAG_TPM*gene_coverage) %>%
  group_by(MAG_subclass, sample, habitat) %>%
  summarise(weighted_As_no_K08356_TPM=sum(weighted_TPM, na.rm=TRUE), .groups="drop") %>%
  complete(MAG_subclass=clade_levels, nesting(sample, habitat),
           fill=list(weighted_As_no_K08356_TPM=0)) %>%
  mutate(MAG_subclass=factor(MAG_subclass, levels=clade_levels),
         habitat=factor(habitat, levels=habitat_levels),
         weighted_As_no_K08356_logTPM=log10(weighted_As_no_K08356_TPM+1))

pA <- ggplot(contig_bar, aes(gene_label, mean_RA, fill=habitat)) +
  geom_col(position=position_dodge(width=0.82), width=0.76, color="grey30", linewidth=0.15) +
  geom_errorbar(aes(ymin=pmax(0, mean_RA-se_RA), ymax=mean_RA+se_RA),
                position=position_dodge(width=0.82), width=0.18, linewidth=0.25) +
  geom_text(data=contig_pair_annot, aes(gene_label, y, label=significance, group=habitat),
            position=position_dodge(width=0.82), inherit.aes=FALSE,
            family="Times", fontface="bold", size=3.2) +
  facet_grid(. ~ pathway, scales="free_x", space="free_x") +
  scale_fill_manual(values=habitat_cols, drop=FALSE) +
  scale_y_continuous(labels=function(x) paste0(x, "%"), expand=expansion(mult=c(0,0.14))) +
  labs(title="A  Contig-level arsenic-gene relative abundance",
       x=NULL, y="Relative abundance (%)", fill="Habitat") +
  base_theme +
  theme(axis.text.x=element_text(angle=55, hjust=1, size=6.1),
        legend.position="top", panel.spacing.x=unit(0.8, "mm"))

pB <- ggplot(contig_long, aes(gene_label, sample, fill=logTPM)) +
  geom_tile(color="white", linewidth=0.06) +
  facet_grid(habitat ~ pathway, scales="free", space="free") +
  scale_fill_gradientn(colors=heat_cols, name="log10(TPM+1)") +
  labs(title="B  Contig-level arsenic-gene abundance across 56 samples",
       x=NULL, y=NULL) +
  base_theme +
  theme(axis.text.x=element_text(angle=55, hjust=1, size=5.8),
        axis.text.y=element_text(size=4.4), strip.text.y=element_text(angle=0),
        panel.spacing=unit(0.7, "mm"), legend.position="right")

pC <- ggplot(clade_sample, aes(habitat, clade_MAG_logTPM, color=MAG_subclass)) +
  geom_jitter(width=0.14, height=0, size=1.25, alpha=0.72) +
  facet_wrap(~MAG_subclass, nrow=1, scales="free_y", labeller=label_wrap_gen(width=24)) +
  scale_color_manual(values=clade_cols, drop=FALSE, guide="none") +
  labs(title="C  Four K08356 host-MAG clade abundance across habitats",
       x=NULL, y="Clade-summed host MAG log10(TPM + 1)") +
  base_theme +
  theme(axis.text.x=element_text(angle=0, size=6.5), strip.text=element_text(size=6.5))

pD <- ggplot(mag_bar, aes(gene_label, mean_RA, fill=habitat)) +
  geom_col(position=position_dodge(width=0.82), width=0.76, color="grey30", linewidth=0.15) +
  geom_errorbar(aes(ymin=pmax(0, mean_RA-se_RA), ymax=mean_RA+se_RA),
                position=position_dodge(width=0.82), width=0.18, linewidth=0.25) +
  geom_text(data=mag_pair_annot, aes(gene_label, y, label=significance, group=habitat),
            position=position_dodge(width=0.82), inherit.aes=FALSE,
            family="Times", fontface="bold", size=3.2) +
  facet_grid(. ~ pathway, scales="free_x", space="free_x") +
  scale_fill_manual(values=habitat_cols, drop=FALSE) +
  scale_y_continuous(labels=function(x) paste0(x, "%"), expand=expansion(mult=c(0,0.14))) +
  labs(title="D  MAG-weighted arsenic-gene relative abundance",
       subtitle="Weighted TPM = host MAG TPM x gene copy number, summed over 47 unique host MAGs",
       x=NULL, y="Relative abundance (%)", fill="Habitat") +
  base_theme +
  theme(plot.subtitle=element_text(size=7),
        axis.text.x=element_text(angle=55, hjust=1, size=6.1),
        legend.position="top", panel.spacing.x=unit(0.8, "mm"))

pE <- ggplot(mag_weighted_long, aes(gene_label, sample, fill=logWeightedTPM)) +
  geom_tile(color="white", linewidth=0.06) +
  facet_grid(habitat ~ pathway, scales="free", space="free") +
  scale_fill_gradientn(colors=heat_cols, name="log10(weighted\nTPM+1)") +
  labs(title="E  MAG-weighted arsenic-gene abundance across 56 samples",
       x=NULL, y=NULL) +
  base_theme +
  theme(axis.text.x=element_text(angle=55, hjust=1, size=5.8),
        axis.text.y=element_text(size=4.4), strip.text.y=element_text(angle=0),
        panel.spacing=unit(0.7, "mm"), legend.position="right")

pF <- ggplot(clade_gene_weighted,
             aes(habitat, weighted_As_no_K08356_logTPM, color=MAG_subclass)) +
  geom_jitter(width=0.14, height=0, size=1.25, alpha=0.72) +
  facet_wrap(~MAG_subclass, nrow=1, scales="free_y", labeller=label_wrap_gen(width=24)) +
  scale_color_manual(values=clade_cols, drop=FALSE, guide="none") +
  labs(title="F  Arsenic-gene weighted abundance carried by four K08356 clades",
       subtitle="K08356 itself is excluded from weighted arsenic-gene abundance",
       x=NULL, y="Clade-carried As genes log10(weighted TPM + 1)") +
  base_theme +
  theme(plot.subtitle=element_text(size=7),
        axis.text.x=element_text(angle=0, size=6.5), strip.text=element_text(size=6.5))

top_left <- plot_grid(pA, pB, pC, ncol=1, rel_heights=c(0.72, 1.04, 0.58), align="v", axis="lr")
bottom_left <- plot_grid(pD, pE, pF, ncol=1, rel_heights=c(0.78, 1.04, 0.58), align="v", axis="lr")
body <- plot_grid(top_left, bottom_left, ncol=1, rel_heights=c(1, 1.02), align="v")
title_panel <- ggdraw() +
  draw_label("Contig- and MAG-resolved arsenic functional potential associated with K08356 evidence subclasses",
             x=0, y=0.72, hjust=0, fontfamily="Times", fontface="bold", size=14) +
  draw_label("A/B/C use contig-level evidence and host-MAG clade abundance; D/E/F use 47 unique K08356 host MAGs. F excludes K08356 itself.",
             x=0, y=0.22, hjust=0, fontfamily="Times", size=8.3)
figure <- plot_grid(title_panel, body, ncol=1, rel_heights=c(0.055, 1))

prefix <- file.path(outdir, "K08356_original_layout_corrected_latest")
ggsave(paste0(prefix, ".pdf"), figure, width=16.8, height=14.2, units="in", device=cairo_pdf)
ggsave(paste0(prefix, ".svg"), figure, width=16.8, height=14.2, units="in", device=svg)
ggsave(paste0(prefix, ".png"), figure, width=16.8, height=14.2, units="in", dpi=600)

write_csv(gene_key, file.path(outdir, "arsenic_gene_key_used.csv"))
write_csv(contig_long, file.path(outdir, "A_B_contig_arsenic_gene_long.csv"))
write_csv(contig_bar, file.path(outdir, "A_contig_habitat_RA_summary.csv"))
write_csv(contig_tests, file.path(outdir, "A_contig_gene_Kruskal_Wallis_BH.csv"))
write_csv(contig_pair_stars, file.path(outdir, "A_contig_IS_vs_other_Wilcoxon_BH.csv"))
write_csv(mag_clade_map, file.path(outdir, "K08356_unique_MAG_clade_map.csv"))
write_csv(mag_gene_long, file.path(outdir, "MAG_arsenic_gene_coverage_used.csv"))
write_csv(mag_weighted_long, file.path(outdir, "D_E_MAG_weighted_arsenic_gene_long.csv"))
write_csv(mag_bar, file.path(outdir, "D_MAG_weighted_habitat_RA_summary.csv"))
write_csv(mag_tests, file.path(outdir, "D_MAG_gene_Kruskal_Wallis_BH.csv"))
write_csv(mag_pair_stars, file.path(outdir, "D_MAG_IS_vs_other_Wilcoxon_BH.csv"))
write_csv(clade_sample, file.path(outdir, "C_four_clade_host_MAG_TPM_by_sample.csv"))
write_csv(clade_gene_weighted, file.path(outdir, "F_four_clade_nonK08356_As_weighted_TPM_by_sample.csv"))

qc <- c(
  paste("metagenomic_samples", nrow(sample_meta), sep="\t"),
  paste("habitat_counts", paste(names(table(sample_meta$habitat)), as.integer(table(sample_meta$habitat)), sep="=", collapse=";"), sep="\t"),
  paste("arsenic_genes", nrow(gene_key), sep="\t"),
  paste("K08356_candidate_sequences", n_distinct(sequence_map$Sequence_ID), sep="\t"),
  paste("unique_K08356_host_MAGs", nrow(mag_clade_map), sep="\t"),
  paste("mixed_clade_MAGs_excluded_from_clade_sums", paste(mag_clade_map$MAG_ID[mag_clade_map$MAG_subclass == "mixed K08356 subclasses"], collapse=";"), sep="\t"),
  paste("F_excludes_K08356", "yes", sep="\t"),
  paste("figure_pdf", paste0(prefix, ".pdf"), sep="\t"),
  paste("figure_svg", paste0(prefix, ".svg"), sep="\t"),
  paste("figure_png", paste0(prefix, ".png"), sep="\t")
)
writeLines(qc, file.path(outdir, "K08356_original_layout_corrected_QC.txt"))
cat(paste(qc, collapse="\n"), "\n")

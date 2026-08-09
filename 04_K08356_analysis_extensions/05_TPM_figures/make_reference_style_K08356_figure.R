suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(stringr)
  library(scales)
  library(grid)
})

outdir <- '/Users/catherine/Downloads/aioA 蛋白序列建树/TPM_combo_figure_20260728/reference_style_K08356'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
contig_xlsx <- '/Users/catherine/Downloads/Contig 热图柱状图 new/all_koid_tpm-final.xlsx'
mag_xlsx <- '/Users/catherine/Downloads/MAG-bin-tpm/副本merged_tpm_matrix.xlsx'
sample_map_tsv <- '/Users/catherine/Downloads/MAG-bin-tpm/derepMAG_abundance_stats/sample_habitat_map.tsv'
strict_qc_tsv <- '/Users/catherine/Downloads/aioA 蛋白序列建树/MAG 补充 idra-tree/joint_reference_search_48/updated_final_evidence_table_joint_reference_strict_QC.tsv'

hab_cols <- c(IS='#66c2a5', AS='#f4a261', ES='#8da0cb', NS='#e76f51')
class_cols <- c(
  'canonical AioA-associated'='#f4a3c4',
  'IdrA-associated with DIRM-like synteny'='#8e44ad',
  'IdrA-associated without complete DIRM-like synteny'='#377eb8',
  'unknown/uncertain DMSOR'='#4daf4a',
  'mixed K08356 subclasses'='#7f7f7f'
)
class_levels <- names(class_cols)

parse_depth <- function(x) {
  m <- str_match(x, '[-_](\\d+)-(\\d+)$')
  (as.numeric(m[,2]) + as.numeric(m[,3])) / 2
}
map_class <- function(x) case_when(
  x == 'canonical AioA-associated' ~ 'canonical AioA-associated',
  x == 'synteny-supported IdrA/DIRM-like neighborhood-associated lineage' ~ 'IdrA-associated with DIRM-like synteny',
  x == 'IdrA-associated phylogenetic lineage' ~ 'IdrA-associated without complete DIRM-like synteny',
  x == 'unknown/uncertain DMSOR' ~ 'unknown/uncertain DMSOR',
  TRUE ~ x
)
fmt_p <- function(p) ifelse(is.na(p), 'NA', ifelse(p < 0.001, '<0.001', sprintf('%.3f', p)))
wrap_subclass <- function(x) {
  stringr::str_wrap(as.character(x), width=26)
}

sample_meta <- read.delim(sample_map_tsv, check.names=FALSE) %>%
  rename(sample=Sample, habitat=Habitat) %>%
  mutate(depth_mid=parse_depth(sample), habitat=factor(habitat, levels=c('IS','AS','ES','NS')))
contig <- read_excel(contig_xlsx, sheet='all_koid_tpm') %>% rename(feature=KO_ID)
mag <- read_excel(mag_xlsx, sheet='调整顺序后的bin的tpm')
names(mag)[1] <- 'MAG_ID'
mag <- mag %>% filter(!is.na(MAG_ID))
strict <- read.delim(strict_qc_tsv, check.names=FALSE)

sample_cols <- intersect(intersect(names(contig), names(mag)), sample_meta$sample)
sample_order <- sample_meta %>% filter(sample %in% sample_cols) %>% arrange(habitat, depth_mid, sample) %>% pull(sample)

# Contig-level total K08356.
contig_k08356 <- contig %>% filter(feature == 'K08356') %>% select(feature, all_of(sample_cols)) %>%
  pivot_longer(-feature, names_to='sample', values_to='TPM') %>%
  left_join(sample_meta, by='sample') %>%
  mutate(logTPM=log10(TPM+1), sample=factor(sample, levels=sample_order))
kw_A <- suppressWarnings(kruskal.test(logTPM ~ habitat, data=contig_k08356))

# Unique MAG classification. Mixed subclass MAG retained once in heatmap, excluded from subclass boxplot.
seq_map <- strict %>% transmute(Sequence_ID, MAG_ID, sequence_subclass=map_class(Final_class), Final_class, neighborhood_class, strict_synteny_QC_pass)
mag_class_map <- seq_map %>% group_by(MAG_ID) %>% summarise(
  n_candidates=n_distinct(Sequence_ID),
  n_subclasses=n_distinct(sequence_subclass),
  MAG_subclass=ifelse(n_subclasses == 1, first(sequence_subclass), 'mixed K08356 subclasses'),
  candidate_ids=paste(unique(Sequence_ID), collapse=';'),
  candidate_subclasses=paste(unique(sequence_subclass), collapse=';'),
  .groups='drop'
) %>% mutate(MAG_subclass=factor(MAG_subclass, levels=class_levels))

host_mag_ids <- unique(mag_class_map$MAG_ID)
mag_host_long <- mag %>% filter(MAG_ID %in% host_mag_ids) %>% select(MAG_ID, all_of(sample_cols)) %>%
  pivot_longer(-MAG_ID, names_to='sample', values_to='TPM') %>%
  left_join(mag_class_map, by='MAG_ID') %>% left_join(sample_meta, by='sample') %>%
  mutate(logTPM=log10(TPM+1), sample=factor(sample, levels=sample_order))

class_sample <- mag_host_long %>% filter(MAG_subclass != 'mixed K08356 subclasses') %>%
  group_by(MAG_subclass, sample, habitat) %>% summarise(TPM=sum(TPM, na.rm=TRUE), logTPM=log10(TPM+1), n_host_MAG=n_distinct(MAG_ID), .groups='drop') %>%
  mutate(MAG_subclass=factor(MAG_subclass, levels=class_levels[1:4]))

class_sample_complete <- class_sample %>%
  complete(MAG_subclass=factor(class_levels[1:4], levels=class_levels[1:4]),
           sample=factor(sample_order, levels=sample_order),
           fill=list(TPM=0, logTPM=0, n_host_MAG=0)) %>%
  left_join(sample_meta %>% select(sample, habitat), by='sample') %>%
  mutate(habitat=coalesce(habitat.x, habitat.y), habitat=factor(habitat, levels=c('IS','AS','ES','NS'))) %>%
  select(MAG_subclass, sample, habitat, TPM, logTPM, n_host_MAG)

relative_abundance <- class_sample_complete %>%
  group_by(sample) %>%
  mutate(sample_total_TPM=sum(TPM, na.rm=TRUE),
         relative_abundance=if_else(sample_total_TPM > 0, TPM/sample_total_TPM, 0)) %>%
  ungroup()

kw_B <- class_sample %>%
  group_by(MAG_subclass) %>%
  summarise(kruskal_p=suppressWarnings(kruskal.test(logTPM ~ habitat)$p.value),
            .groups='drop') %>%
  mutate(BH_FDR=p.adjust(kruskal_p, method='BH'),
         p_label=paste0('BH-FDR = ', fmt_p(BH_FDR)),
         MAG_subclass=factor(MAG_subclass, levels=class_levels[1:4]))

kw_B_pos <- class_sample %>%
  group_by(MAG_subclass) %>%
  summarise(x=1.05, y=max(logTPM, na.rm=TRUE)*0.96, .groups='drop') %>%
  left_join(kw_B, by='MAG_subclass')

# Heatmap with 47 unique host MAG rows.
heat_df <- mag_host_long %>% mutate(MAG_subclass=factor(MAG_subclass, levels=class_levels), sample=factor(sample, levels=sample_order))
row_order <- heat_df %>% group_by(MAG_ID, MAG_subclass) %>% summarise(mean_log=mean(logTPM, na.rm=TRUE), .groups='drop') %>% arrange(MAG_subclass, desc(mean_log), MAG_ID) %>% pull(MAG_ID)
heat_df <- heat_df %>% mutate(MAG_ID=factor(MAG_ID, levels=rev(unique(row_order))))

# Compact heatmap for contig K08356 as a one-row tile strip.
contig_heat <- contig_k08356 %>% mutate(row='K08356', sample=factor(sample, levels=sample_order))

base_theme <- theme_classic(base_family='Times', base_size=9) +
  theme(plot.title=element_text(face='bold', size=11), axis.text.x=element_text(size=7), axis.text.y=element_text(size=7), legend.title=element_text(size=8), legend.text=element_text(size=7), plot.margin=margin(3,4,3,4))

pA <- ggplot(contig_k08356, aes(habitat, logTPM, fill=habitat)) +
  geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.72, color='grey25', linewidth=0.28) +
  geom_jitter(aes(color=habitat), width=0.12, size=1.75, alpha=0.85, show.legend=FALSE) +
  scale_fill_manual(values=hab_cols, drop=FALSE) +
  scale_color_manual(values=hab_cols, drop=FALSE) +
  annotate('text', x=Inf, y=Inf, label=paste0('Kruskal-Wallis p = ', fmt_p(kw_A$p.value)),
           hjust=1.05, vjust=1.5, size=3.0, family='Times') +
  labs(title='A  Contig-level K08356 abundance', x=NULL, y='log10(TPM + 1)') +
  base_theme +
  theme(legend.position='none')

pC <- ggplot(contig_heat, aes(sample, row, fill=logTPM)) +
  geom_tile(color='white', linewidth=0.18) +
  scale_fill_gradientn(colors=c('#2166ac','#f7f7f7','#f4a582','#b2182b'), name='log10(TPM+1)') +
  labs(title='B  Contig-level K08356 TPM', x=NULL, y=NULL) + base_theme +
  theme(axis.text.x=element_text(angle=65, hjust=1, size=5.5), axis.text.y=element_text(size=8), legend.position='right', plot.title=element_text(size=10))

pB <- ggplot(class_sample, aes(habitat, logTPM, fill=habitat)) +
  geom_boxplot(width=0.56, outlier.shape=NA, alpha=0.72, color='grey25', linewidth=0.25) +
  geom_jitter(aes(color=habitat), width=0.12, height=0, size=1.25, alpha=0.68, show.legend=FALSE) +
  geom_text(data=kw_B_pos, aes(x=x, y=y, label=p_label), inherit.aes=FALSE,
            hjust=0, vjust=1, family='Times', size=2.6) +
  facet_wrap(~MAG_subclass, ncol=2, scales='free_y', labeller=label_wrap_gen(width=24)) +
  scale_fill_manual(values=hab_cols, drop=FALSE) + scale_color_manual(values=hab_cols, drop=FALSE) +
  labs(title='D  Host-MAG abundance by habitat and K08356 subclass',
       x=NULL, y='Subclass-summed host MAG log10(TPM + 1)') + base_theme +
  theme(axis.text.x=element_text(angle=0, hjust=0.5, size=6.4),
        strip.background=element_rect(fill='grey94', color='grey80'),
        strip.text=element_text(face='bold', size=7.0),
        legend.position='top')

pD <- ggplot(heat_df, aes(sample, MAG_ID, fill=logTPM)) +
  geom_tile(color='white', linewidth=0.06) + facet_grid(MAG_subclass ~ ., scales='free_y', space='free_y', labeller=label_wrap_gen(width=32)) +
  scale_fill_gradientn(colors=c('#2166ac','#f7f7f7','#f4a582','#b2182b'), name='log10(TPM+1)') +
  labs(title='C  47 unique K08356-bearing host MAGs', x=NULL, y=NULL) + base_theme +
  theme(strip.background=element_rect(fill='grey94', color='grey80'), strip.text.y=element_text(angle=0, face='bold', size=6.7), axis.text.x=element_text(angle=65, hjust=1, size=5.2), axis.text.y=element_text(size=4.35), legend.position='right')

left <- plot_grid(pA, pC, pD, ncol=1, rel_heights=c(0.58,0.30,1.48), align='v')
body <- plot_grid(left, pB, ncol=2, rel_widths=c(1.72,1.0), align='h')
title <- ggdraw() +
  draw_label('K08356 abundance at contig and host-MAG levels', x=0, y=0.72, hjust=0, fontfamily='Times', fontface='bold', size=15) +
  draw_label('Reference-style layout with corrected statistical units: relative abundance and tests use four K08356 subclasses after excluding mixed-subclass MAGs; heatmap keeps each host MAG once.', x=0, y=0.24, hjust=0, fontfamily='Times', size=9)
fig <- plot_grid(title, body, ncol=1, rel_heights=c(0.09,1))

pdf_file <- file.path(outdir, 'reference_style_K08356_abundance.pdf')
png_file <- file.path(outdir, 'reference_style_K08356_abundance.png')
svg_file <- file.path(outdir, 'reference_style_K08356_abundance.svg')
ggsave(pdf_file, fig, width=15.2, height=12.2, units='in', device=cairo_pdf)
ggsave(png_file, fig, width=15.2, height=12.2, units='in', dpi=300)
ggsave(svg_file, fig, width=15.2, height=12.2, units='in', device=svg)

write.csv(contig_k08356, file.path(outdir,'contig_K08356_TPM_by_sample.csv'), row.names=FALSE)
write.csv(mag_class_map, file.path(outdir,'unique_MAG_subclass_map_47.csv'), row.names=FALSE)
write.csv(class_sample, file.path(outdir,'subclass_unique_MAG_TPM_by_sample.csv'), row.names=FALSE)
write.csv(relative_abundance, file.path(outdir,'subclass_relative_abundance_by_sample.csv'), row.names=FALSE)
write.csv(kw_B, file.path(outdir,'subclass_host_MAG_Kruskal_BH.csv'), row.names=FALSE)
write.csv(heat_df %>% select(MAG_ID, MAG_subclass, sample, habitat, TPM, logTPM, n_candidates, candidate_ids, candidate_subclasses), file.path(outdir,'unique_MAG_heatmap_long.csv'), row.names=FALSE)
qc <- c(
  paste('shared_samples', length(sample_cols), sep='\t'),
  paste('contig_K08356_samples', nrow(contig_k08356), sep='\t'),
  paste('K08356_candidate_sequences', nrow(seq_map), sep='\t'),
  paste('unique_host_MAGs', nrow(mag_class_map), sep='\t'),
  paste('mixed_subclass_MAG_retained_once', paste(mag_class_map$MAG_ID[mag_class_map$MAG_subclass=='mixed K08356 subclasses'], collapse=';'), sep='\t'),
  paste('mixed_subclass_MAG_excluded_from_relative_abundance_and_tests', paste(mag_class_map$MAG_ID[mag_class_map$MAG_subclass=='mixed K08356 subclasses'], collapse=';'), sep='\t'),
  paste('figure_pdf', pdf_file, sep='\t')
)
writeLines(qc, file.path(outdir,'reference_style_K08356_QC.txt'))
cat(paste(qc, collapse='\n'), '\n')

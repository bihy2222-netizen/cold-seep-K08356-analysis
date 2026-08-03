# Clade 4 IS association: manuscript-ready text

## English results

Across 56 environmental samples, host MAGs carrying Clade 4 K08356 sequences
were detected in 20 of 21 IS samples (95.2%), but in only 1 of 35 non-IS samples
(2.9%; Fisher's exact test, P = 5.46 x 10^-13; conditional odds ratio = 415.5,
exact 95% CI = 32.7-28,834.5). Clade 4 host-MAG TPM differed significantly among
habitats (Kruskal-Wallis, P = 4.78 x 10^-10), with significantly higher abundance
in IS than in AS, ES, and NS after Benjamini-Hochberg correction. These results
identify Clade 4 as a strongly IS-associated genomic lineage rather than an
IS-exclusive lineage.

Because all 20 Clade 4 host genomes were initially assembled from IS, this
comparison specifically demonstrates enrichment in IS of populations matching
these IS-derived genomes. It does not exclude more divergent Clade 4 populations
in other habitats and does not provide evidence of metabolic activity or
elemental-cycle flux.

## Chinese results

在 56 个环境样本中，携带 Clade 4 K08356 序列的宿主 MAG 在 21 个 IS 样本中的
20 个被检出（95.2%），而在 35 个非 IS 样本中仅有 1 个被检出（2.9%；Fisher
精确检验，P = 5.46 x 10^-13；条件优势比 = 415.5，精确 95% CI =
32.7-28,834.5）。Clade 4 宿主 MAG 的 TPM 在四种生境间存在显著差异
（Kruskal-Wallis，P = 4.78 x 10^-10）；经 Benjamini-Hochberg 校正后，IS 中的
丰度均显著高于 AS、ES 和 NS。该结果支持将 Clade 4 定义为与 IS 强烈相关的
基因组谱系，而非 IS 特异谱系。

由于全部 20 个 Clade 4 宿主基因组最初均由 IS 样本组装获得，上述比较严格
证明的是：与这些 IS 来源基因组相匹配的种群在 IS 中富集。该结果不能排除
AS、ES 或 NS 中存在因序列差异而未被当前参考捕获的远缘 Clade 4 成员，也不
代表代谢活性或元素循环通量。

## Methods and figure-legend wording

For each sample `s` and metabolic feature `m`, the abundance-weighted genomic
potential was calculated as `W_sm = sum(TPM_si * Score_im) / sum(TPM_si)` across
detected host MAGs. Habitat-level colors are equal-weight means of `W_sm` across
samples in which at least one host was detected. Samples without host detection
are treated as missing (`NA`), not as a metabolic score of zero. Bubble size is
the sample detection prevalence of host MAGs with >=50% partial reconstruction.
Only habitat-feature combinations with prevalence >0 are plotted. Blank cells
indicate that no qualifying host MAG was detected; they do not demonstrate
functional absence from the habitat.

Host detection required TPM >0 after CoverM filtering with >=10% covered
fraction, >=95% read identity, >=75% aligned-read fraction, and 0.1/0.9 read-end
trimming. The current nonparametric tests do not adjust for station/core or
depth; a mixed-effects sensitivity analysis is still required if these samples
contain repeated depth layers from the same station or core.

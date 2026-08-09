# Overall dRepMAG figure notes

This folder contains the finalized overall dRepMAG figures for K08356-bearing host MAG abundance.

Main statistical unit: sample.
Input abundance: `merged_tpm_matrix.txt` MAG TPM. Each sample column is already TPM-normalized, with column sums close to 1,000,000.
Biological meaning: MAG TPM describes the abundance of host derepMAGs carrying K08356-related genes; it is not K08356 gene/contig TPM.
Branch panels: if one derepMAG contains multiple K08356 copies assigned to different branches, that MAG TPM is split equally among its K08356 copy rows before branch-level summation.
Significance: Kruskal-Wallis tests are run on sample-level log10(TPM + 1). Pairwise Wilcoxon tests use BH adjustment and are retained in the CSV/Excel tables.

Key results for wording:
- Total K08356-bearing host MAG TPM differs significantly among habitats: Kruskal-Wallis p = 1.32e-6.
- aioA-like-associated host MAG TPM differs significantly among habitats: Kruskal-Wallis p = 2.20e-7; IS is significantly higher than AS, ES, and NS in pairwise BH-adjusted tests.
- unknown AioA-like / uncertain DMSOR host MAG TPM differs significantly among habitats: Kruskal-Wallis p = 4.01e-4; IS is higher than AS and ES, but IS vs NS is not significant after BH adjustment.
- IdrA-associated host MAG TPM is not significant among habitats after copy-number adjustment: Kruskal-Wallis p = 0.121.
- canonical aioA-associated host MAG TPM is not significant among habitats: Kruskal-Wallis p = 0.877.

Recommended placement:
- `overall_derepMAG_total_boxplot_square` and `overall_derepMAG_branch_boxplot_square` are the cleanest main Figure A/B files.
- `overall_derepMAG_main_figure` is a one-file combined overview, but it has one blank facet slot because it contains five panels.
- `overall_derepMAG_branch_sample_order_stacked_bar` is better used as a supplementary/sample-order distribution figure; its x-axis order matches `/Users/catherine/Downloads/metaphlan 作图/class34.csv 作图/最终版图片/class_level_grouped_final_v2.svg` and the previous grouped-analysis order.
- Grouped-derep results should be added later as an appendix/sensitivity analysis, not as a replacement for the overall dRepMAG abundance figure.

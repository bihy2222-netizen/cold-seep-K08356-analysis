# Group-derep MAG48 abundance-weighted metabolism

## Validated inputs

- 49 K08356 sequence memberships
- 48 unique host MAGs
- clade counts 3/7/19/20
- 56 metagenomic samples
- 2,688 unique sample x target-MAG TPM rows
- 48 x 40 locked metabolic score matrix

The merged source contains four separate group references: IS 315 MAGs, AS 128,
ES 252, and NS 156. Each reference was profiled against all 56 samples. Every
target MAG occurs in one source reference and has one TPM value in every sample.

## Main figure semantics

`group_derep49_MAG48_abundance_weighted_metabolism` uses:

- fill: the mean of within-sample TPM-weighted metabolic scores. For sample
  `s` and feature `m`, `W_sm = sum(TPM_si * Score_im) / sum(TPM_si)` across
  detected host MAGs, followed by an equal-weight mean across detected samples
  in each habitat;
- size: sample detection prevalence of host MAGs with >=50% partial
  reconstruction;
- blank cell: no sample detected a host MAG reaching the >=50% threshold for
  that feature. Its plotting prevalence and weighted score are both `NA`, and
  no point is drawn. Blank cells do not demonstrate functional absence from
  the habitat.

Only the 116 of 640 clade x habitat x feature combinations with positive
qualifying-host prevalence are sent to `geom_point()`. `scale_size_area()` maps
area from prevalence and omits the zero break, so a true low-frequency signal
(for example, 1/14 = 7.14%) remains visible while zero detection remains blank.

The exact bubble-size legend is: `Sample detection prevalence of host MAGs with
>=50% partial reconstruction (%)`.

Host detection is defined as TPM >0 after CoverM filtering with >=10% covered
fraction, >=95% read identity, >=75% aligned-read fraction, and 0.1/0.9 read-end
trimming. Thus TPM >0 is not an unfiltered single-read criterion.

The unweighted locked 50% figure is copied beside it with an explicit
`unweighted` filename.

`S1_9-12_bin1` contributes to Clades 1 and 3. The main figure keeps both
memberships to reproduce 3/7/19/20, and
`dual_clade_MAG_exclusion_abundance_sensitivity.tsv` reports the exclusion
sensitivity.

## Reference-consistency boundary

- Clades 1-3 combine MAGs quantified in three different background references.
  Their abundance-weighted panels are descriptive and cannot establish habitat
  enrichment.
- All 20 Clade 4 MAGs occur in the IS reference. Their comparison across the 56
  samples therefore uses one consistent reference.

Clade 4 host MAGs were detected in 20/21 IS samples, 1/14 AS samples, 0/13 ES
samples, and 0/8 NS samples. The sample-level Kruskal-Wallis test on
`log10(summed TPM + 1)` gave p=4.778e-10. Pairwise Wilcoxon tests with BH
correction supported IS versus AS, ES, and NS, while AS/ES/NS contrasts were not
significant.

For detection versus nondetection, Clade 4 occurred in 20/21 IS samples (95.2%)
and 1/35 non-IS samples (2.9%). Fisher's exact test gave p=5.465e-13, with a
conditional odds ratio of 415.5 and an exact 95% confidence interval of
32.7-28,834.5.

This identifies Clade 4 as strongly IS-associated, not IS-exclusive. Because all
20 host genomes were initially assembled from IS, the strict interpretation is
that populations matching these IS-derived Clade 4 genomes are enriched in IS.
The analysis cannot exclude more divergent Clade 4 populations in AS, ES, or NS.
It also does not establish habitat causality, spatial heterogeneity independent
of station/core/depth, metabolic activity, or elemental-cycle flux.

## Outputs

- `results/group_derep49_MAG48_abundance_weighted_metabolism.{png,pdf,svg}`
- `results/Clade4_same_reference_host_MAG_TPM_by_habitat.{png,pdf,svg}`
- `results/abundance_weighted_clade_habitat_feature_summary.tsv`
- `results/Clade4_same_reference_Kruskal_Wallis.tsv`
- `results/Clade4_same_reference_pairwise_Wilcoxon_BH.tsv`
- `results/Clade4_detection_prevalence_by_habitat.tsv`
- `results/Clade4_IS_vs_nonIS_Fisher_exact.tsv`
- `results/figure_plotting_cells_positive_prevalence.tsv`
- `results/clade_reference_consistency_audit.tsv`
- `results/ABUNDANCE_WEIGHTED_METABOLISM_QC.txt`
- `results/MANUSCRIPT_RESULT_CLADE4_IS_ASSOCIATION.md`

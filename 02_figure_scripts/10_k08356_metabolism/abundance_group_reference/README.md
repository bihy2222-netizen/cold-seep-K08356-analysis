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

## Run

```bash
export GROUP_TPM_LONG=/path/to/IS_AS_ES_NS_MAG_TPM_long.tsv
export GROUP_TPM_SUMMARY=/path/to/IS_AS_ES_NS_MAG_TPM_summary.tsv
export SAMPLE_METADATA_CSV=/path/to/sample_metadata_used.csv
Rscript plot_abundance_weighted_metabolism_group_derep49.R
```

`METABOLIC_SCORE_TSV`, `K08356_MAPPING_TSV`, and
`ABUNDANCE_METABOLISM_OUT_DIR` can also override their repository-relative
defaults. The long TPM input is not committed.

## Main figure semantics

`group_derep49_MAG48_abundance_weighted_metabolism` uses:

- fill: TPM-weighted mean metabolic score;
- size: percentage of samples in a habitat with TPM >0 for at least one host
  MAG passing the legacy >=50% partial-reconstruction threshold;
- hollow point: no sample-level detection of a >=50% functional host.

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

This supports a provisional, reference-consistent concentration of Clade 4 host
abundance in IS samples. It does not yet establish habitat causality, spatial
heterogeneity independent of station/depth, activity, or elemental-cycle flux.

## Outputs

- `results/group_derep49_MAG48_abundance_weighted_metabolism.{png,pdf,svg}`
- `results/Clade4_same_reference_host_MAG_TPM_by_habitat.{png,pdf,svg}`
- `results/abundance_weighted_clade_habitat_feature_summary.tsv`
- `results/Clade4_same_reference_Kruskal_Wallis.tsv`
- `results/Clade4_same_reference_pairwise_Wilcoxon_BH.tsv`
- `results/clade_reference_consistency_audit.tsv`
- `results/ABUNDANCE_WEIGHTED_METABOLISM_QC.txt`

# GTDB family adjustment interpretation

## Exact taxonomy set

GTDB-Tk release r226 taxonomy is available for all 48 current K08356-bearing
MAGs. Forty-seven IDs come from the existing all-derep classification. The
group-derep-only `R2111_N500_0-10_bin13` was rerun with GTDB-Tk 2.4.1 and r226,
then merged by exact MAG ID. It is classified as
Actinomycetota/Acidimicrobiia/UBA5794/ZC4RG35/JAWWBL01.

## Confounding structure

- 21 GTDB families are assigned across 48 MAGs.
- 15 families contain only one MAG.
- Three families occur in more than one habitat.
- Two families occur in more than one clade; one is the dual-clade host
  `S1_9-12_bin1` by itself.
- Rhodobacteraceae contains 22 MAGs, all from IS: five Clade 3 and 17 Clade 4.
- After excluding `S1_9-12_bin1`, the completeness + habitat + clade + family
  model matrix has rank 25 of 27 and 22 residual degrees of freedom.
- Habitat retains three independently estimable degrees of freedom, whereas
  clade retains only one after the other covariates are included.

These properties mean that family-aware models can diagnose confounding but
cannot fully separate host family, clade, and habitat effects.

## Locked-feature sensitivity

Within IS Rhodobacteraceae only:

- `arsC`: Clade 3 60.0%, Clade 4 82.35%; Fisher p=0.548.
- Independent `arrA`: absent from both groups.
- Independent `arxA/aioA`: absent from both groups.
- Cobinamide-to-cobalamin reconstruction: mean 72.5% in Clade 3 and 78.68% in
  Clade 4; Wilcoxon p=0.533.
- The completeness-adjusted clade coefficients are also not significant for
  `arsC` or B12 after BH correction.

This does not prove equivalence. The Clade 3 group contains only five MAGs, and
the analysis addresses genomic potential rather than activity.

## Permitted interpretation

The exact GTDB family table and family-aware sensitivity are complete. The
locked Clade 4-IS arsenic/B12 description remains valid as a descriptive host
metabolic-potential pattern, but it is not a demonstrated clade effect independent
of host family. Habitat enrichment and spatial heterogeneity remain deferred
until the unified 56-sample abundance matrix is available.

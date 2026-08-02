# Locked manuscript wording

## Analysis status

- Metabolic-potential figure: verified.
- Abundance/TPM figure: not verified.
- Completeness sensitivity: verified; broad pathway scores remain sensitive.
- Exact GTDB family table and family-aware sensitivity: verified; residual
  family-habitat-clade confounding remains.

The current figure is interpreted as an unadjusted association in genomic
metabolic potential. It does not demonstrate activity, elemental-cycle coupling,
or a habitat-driven effect.

## Canonical AioA evidence

Two of the three phylogenetically assigned canonical AioA sequences exceeded the
METABOLIC `aioA.hmm` threshold of 800. The third obtained a METABOLIC score of
762.0 and was retained based on independent sequence identity, phylogenetic
placement, and an independent joint-screen `Combined_score` of 969.5 (threshold
>=640), rather than METABOLIC HMM detection.

## Legacy reconstruction thresholds

For continuity with the pre-specified legacy algorithm, the main visualization
used a >=50% partial gene-set-reconstruction threshold. A >=75% threshold was
evaluated as a sensitivity analysis. Neither threshold was interpreted as proof
of a complete pathway, and no mechanistic inference relied solely on the 50%
cutoff.

## Clade 4-IS result

Among the 20 Clade 4-IS MAGs, `arsC` was detected in 17/20 (85%), whereas `arrA`
was detected in 2/20 (10%). No independent `arxA`/`aioA` marker remained after
excluding the focal K08356 sequences. All 20 MAGs exceeded the legacy 50%
partial-reconstruction threshold for cobinamide-to-cobalamin biosynthesis, with
a mean gene-set coverage of 79.375%.

Clade 4 MAGs showed widespread arsenic detoxification potential, whereas
respiratory arsenate reduction was uncommon and no independent
arsenite-oxidation marker remained after exclusion of the focal K08356
sequences.

## Sensitivity boundaries

Removing the 19 KOs unavailable in the current KO result catalog from compatible
denominators caused 41 MAG-feature threshold flips at 50% and 10 at 75%. The
flips were confined to TCA, second-stage TCA, pyruvate oxidation, reductive TCA,
denitrification, and urea-cycle rows. Threshold-dependent interpretation of these
rows is therefore not robust to annotation-catalog compatibility.

Excluding the dual-copy MAG `S1_9-12_bin1` made the n=1 Clade 1-AS and Clade 3-AS
cells empty. Clade 4-IS and all other populated cells were unchanged.

Completeness scaling under a random-gene-loss assumption caused 83
MAG-feature threshold flips at 50% and 35 at 75%. The locked Clade 4-IS
arsenic-marker prevalences and B12 prevalence at 50% were unchanged, but this
scaling is a diagnostic sensitivity rather than an unbiased correction.

GTDB-Tk r226 family taxonomy matched all 48 MAG IDs. Fifteen of 21 families were
represented by one MAG, only three families spanned habitats, and only two
spanned clades. Within the 22 IS Rhodobacteraceae MAGs, no locked arsenic/B12
feature differed significantly between Clade 3 (n=5) and Clade 4 (n=17), before
or after completeness adjustment. These analyses do not establish that host
family or completeness confounding has been eliminated.

Habitat enrichment, IS hotspot, and spatial-heterogeneity wording remains
deferred until all 56 samples have been mapped to one common 48-MAG derep
reference and the resulting abundance matrix has passed its ID and method audit.

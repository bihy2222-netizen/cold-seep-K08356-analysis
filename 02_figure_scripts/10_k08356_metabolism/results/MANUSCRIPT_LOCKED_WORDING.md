# Locked manuscript wording

## Analysis status

- Metabolic-potential figure: verified.
- Abundance/TPM figure: not verified.
- Completeness and GTDB-family adjustment: not completed.

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

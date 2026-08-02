# Completeness sensitivity interpretation

## Scope

All 48 MAGs pass the predefined completeness >=50% and contamination <=10%
filter, but completeness remains imbalanced among habitats and clades. This
supplement asks whether observed gene-set scores are sensitive to that remaining
variation. It does not alter the locked main plot.

For multi-gene sets, the sensitivity value is:

`min(100, observed gene-set coverage / MAG completeness fraction)`

This assumes that missing genes are lost at random as MAG completeness declines.
That assumption may be false for any individual pathway, so adjusted scores are
diagnostic rather than corrected ground truth. Single-gene arsenic-marker
carriage remains binary; completeness cannot justify imputing an absent marker.

## Results

- The exact matrix remains 48 MAG x 40 features (1,920 unique rows).
- Completeness scaling changes 83 MAG-feature classifications at the 50%
  threshold and 35 at the 75% threshold.
- Descriptive OLS models use raw score as outcome and include completeness,
  habitat, and clade after excluding the dual-clade host `S1_9-12_bin1`.
- Ten of 40 features retain an association with completeness after BH correction.
  These include denitrification, dissimilatory sulfate reduction, several carbon
  gene sets, nitrite plus ammonia to nitrogen, pyruvate oxidation, and urea
  cycle. These rows require caution in between-group interpretation.
- The three Clade 4-IS arsenic-marker prevalences do not change.
- Clade 4-IS cobinamide-to-cobalamin B12 reconstruction remains 20/20 at the
  legacy 50% threshold; its mean changes from 79.375% raw to 86.556% under the
  sensitivity scaling. At 75%, prevalence changes from 90% to 95%.

## Permitted interpretation

The locked Clade 4-IS arsenic/B12 result is robust at the prespecified 50%
threshold. The main score matrix remains an unadjusted genomic-potential
association, and broad pathway differences should not be described as free of
completeness effects. Family-level taxonomy is evaluated separately.

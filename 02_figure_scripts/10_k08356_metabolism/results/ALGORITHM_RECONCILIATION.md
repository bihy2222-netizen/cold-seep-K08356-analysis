# Previous-algorithm reconciliation

## Formula retained

For every unique MAG and predefined gene set:

`gene-set coverage (%) = 100 * number of distinct genes with count > 0 /
number of distinct genes in the set`

For every habitat-clade cell, fill is the arithmetic mean of the MAG-level
coverage values. Bubble size is the percentage of unique MAGs at or above the
selected coverage threshold. The main threshold is 50%; 75% is a sensitivity
analysis. Gene copy number does not increase coverage.

## Version changes

The PDF example used 48 K08356 sequences from 47 MAGs. The current canonical
set has 49 sequences from 48 MAGs because
`R2111_N500_0-10_bin13-k141_4404301_3` was added to Clade 2. The formula is
unchanged, but all values are recalculated from the current 48-MAG KO calls.

## Regression examples

- Clade 4-IS, `Cobalamin biosynthesis, cobinamide => cobalamin`: the PDF
  reported 100% MAG prevalence and 78.1% mean coverage. The current run gives
  100% and 79.375%, respectively.
- Clade 4-IS, arsenate reduction: the PDF example reported approximately 80%
  prevalence and 40% mean coverage. Current HMM-component scoring gives 85%
  and 47.5%.
- Arsenite oxidation is scored only after excluding all focal K08356 IDs. Two
  focal `aioA` hits are excluded; no non-focal arxA/aioA hit remains, so the
  current arsenite-oxidation score is zero.

The small numerical differences are annotation/version effects, not a change
to the MAG-level averaging or prevalence formulas.

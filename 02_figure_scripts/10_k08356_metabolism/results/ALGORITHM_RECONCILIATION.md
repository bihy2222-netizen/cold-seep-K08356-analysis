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
- Clade 4-IS, the revised binary marker rows give 10% `arrA` carriage for
  respiratory arsenate reduction and 85% `arsC` carriage for arsenate
  detoxification/resistance. These biological meanings are no longer combined.
- Arsenite oxidation is scored only after excluding all focal K08356 IDs. Two
  focal `aioA` hits are excluded; no non-focal arxA/aioA hit remains, so the
  current arsenite-oxidation score is zero.
- The third canonical sequence is present in the scoring audit and returns a
  raw METABOLIC `aioA.hmm` score of 762.0. It was omitted from the original
  `tblout` because the run used `-T 800`; this is not an identifier mismatch.

The B12 difference is an annotation/version effect, not a change to the
MAG-level averaging or prevalence formulas. The arsenic difference also
reflects the explicit switch from a multi-gene color interpretation to binary
independent-marker carriage.

## Definition-catalog audit

The old table contains 508 KO definition records and two non-KO EC records.
`EC:1.20.4.1` denotes ArsC-type detoxification/resistance and `EC:1.20.99.1`
denotes donor-dependent ArrA-type respiratory arsenate reduction. The EC records
are confined to the legacy arsenic set, which is replaced in the
figure by explicit HMM-marker rows. Thirty KO definition records (19 distinct
KOs) do not occur in the current per-MAG KO result catalog. The main figure keeps
them in the denominator for exact old-algorithm reproduction; affected coverage
values are conservative and the impact is enumerated in the output audit table.

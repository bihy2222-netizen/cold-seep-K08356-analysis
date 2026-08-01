# Publication-risk revision record

## Resolved in the current figure

1. The canonical count is explicitly fixed at 49 K08356 sequences from 48
   unique MAGs, with clade counts 3/7/19/20. Older mappings are audited and
   marked superseded.
2. Focal K08356 IDs are excluded from binary per-MAG arsenic HMM-marker
   carriage. This removes two circular `aioA` hits; no non-focal arxA/aioA
   marker remains. The third canonical sequence scores 762.0 against the same
   HMM and was below the original `-T 800` reporting threshold, not mismatched.
3. Arsenic markers are split by interpretation: respiratory arsenate reduction
   (`arrA`), arsenate detoxification/resistance (`arsC`), and arsenite oxidation
   (`arxA`/`aioA`). They are not labelled as multi-gene coverage.
4. Bubble prevalence is calculated after deduplicating host MAGs within each
   habitat-clade-feature cell.
5. Each panel shows sequence and unique-MAG sample sizes. One-MAG cells are
   marked as descriptive only.
6. Gray hatching and `NA: no members` denote no clade member; hollow zero points
   denote members with no reconstructed feature.
7. S13 and S15 are assigned to ES using both the original sample metadata and
   the group-derep input manifest. Two stale classification labels are logged.
8. A 75% gene-set-coverage sensitivity figure accompanies the 50% main figure.
9. The calculation has been restored to the previous distinct-gene set
   algorithm documented in `MAG和Contig结果核验.pdf`; it does not use METABOLIC
   module-step coverage.
10. Flagellar Assembly is labelled as partial 36-KO reconstruction and is not
   interpreted as complete motility.
11. Output names distinguish 49 focal sequences from 48 host MAGs.
12. Captions limit interpretation to genomic potential and state that gene sets
    can overlap and are not independent.

## Remaining limitation

The legacy 510-record denominator is reproduced exactly. It contains 508 KO
records and two EC identifiers. Thirty records (19 distinct KOs) are absent from
the current KO result catalog and therefore cannot score above zero in this run;
affected gene-set coverage is conservatively depressed. The exact IDs and
feature-level denominator impact are reported in the definition audit tables.

Completeness and contamination are available and audited for all 48 MAGs from
the source metaWRAP tables. Completeness is uneven by habitat and clade, so the
50/10 filter alone does not eliminate this confounding. A validated family-level GTDB taxonomy table was not
found in the inspected group-dRep or METABOLIC paths. Therefore, a family- and
completeness-adjusted clade model is not reported. It should only be added after
the exact 48 MAG IDs are matched to a GTDB family table; CheckM lineage is too
coarse to substitute for that analysis.

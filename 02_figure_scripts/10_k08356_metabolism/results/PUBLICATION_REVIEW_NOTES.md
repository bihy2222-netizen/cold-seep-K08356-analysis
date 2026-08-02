# Publication-risk revision record

## Resolved in the current figure

1. The canonical count is explicitly fixed at 49 K08356 sequences from 48
   unique MAGs, with clade counts 3/7/19/20. Older mappings are audited and
   marked superseded.
2. Focal K08356 IDs are excluded from binary per-MAG arsenic HMM-marker
   carriage. This removes two circular `aioA` hits; no non-focal arxA/aioA
   marker remains. The third canonical sequence scores 762.0 against the same
   HMM and was below the original `-T 800` reporting threshold, not mismatched.
   Its independent joint-screen `Combined_score` is 969.5 (>=640), and its
   canonical assignment is supported by sequence identity and phylogeny rather
   than METABOLIC HMM detection.
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
   Both are partial-reconstruction thresholds and do not establish a complete
   pathway; no mechanistic claim relies solely on the 50% threshold.
9. The calculation has been restored to the previous distinct-gene set
   algorithm documented in `MAG和Contig结果核验.pdf`; it does not use METABOLIC
   module-step coverage.
10. Flagellar Assembly is labelled as partial 36-KO reconstruction and is not
   interpreted as complete motility.
11. Output names distinguish 49 focal sequences from 48 host MAGs.
12. Captions limit interpretation to genomic potential and state that gene sets
    can overlap and are not independent.
13. Excluding dual-copy MAG `S1_9-12_bin1` removes only the descriptive n=1
    Clade 1-AS and Clade 3-AS cells. Clade 4-IS and every other populated cell
    remain unchanged.
14. The locked Clade 4-IS result is `arsC` 17/20, `arrA` 2/20, independent
    `arxA`/`aioA` 0/20, and cobinamide-to-cobalamin B12 partial reconstruction
    >=50% in 20/20 with 79.375% mean coverage.
15. Completeness sensitivity is complete. Random-gene-loss scaling produces 83
    threshold flips at 50% and 35 at 75%; the locked Clade 4-IS 50% result is
    unchanged.
16. GTDB-Tk r226 taxonomy is matched to the exact 48 MAG IDs. The added
    `R2111_N500_0-10_bin13` is classified as family ZC4RG35 using GTDB-Tk 2.4.1
    and the same r226 database.
17. Family-aware design-rank and within-Rhodobacteraceae sensitivities are
    reported. They diagnose strong host-family confounding but do not remove it.

## Remaining limitation

The legacy 510-record denominator is reproduced exactly. It contains 508 KO
records and two EC identifiers. Thirty records (19 distinct KOs) are absent from
the current KO result catalog and therefore cannot score above zero in this run;
affected gene-set coverage is conservatively depressed. The exact IDs and
feature-level denominator impact are reported in the definition audit tables.
The bias is uneven: a catalog-compatible denominator causes 41 MAG-feature
classification flips at 50% and 10 at 75%. These occur in TCA, second-stage TCA,
pyruvate oxidation, reductive TCA, denitrification, and urea-cycle rows. The
arsenic-marker rows and locked Clade 4-IS B12 result are unaffected, but the
threshold-sensitive carbon/nitrogen rows cannot support standalone mechanistic
claims.

Completeness and contamination are available and audited for all 48 MAGs from
the source metaWRAP tables. Completeness is uneven by habitat and clade, so the
50/10 filter alone does not eliminate this confounding. Ten of 40 raw feature
scores remain associated with completeness after BH correction in descriptive
models containing habitat and clade.

The exact GTDB table contains 21 assigned families, but 15 are singletons, only
three cross habitats, and only two cross clades. The full completeness + habitat
+ clade + family design has rank 25/27; only one clade degree of freedom is
independent of the other covariates. Rhodobacteraceae is especially entangled:
all 22 representatives are from IS, with five in Clade 3 and 17 in Clade 4.
Within that family and habitat, the locked features do not differ significantly,
but the five-MAG Clade 3 group limits power. These results permit a family-aware
sensitivity statement, not a claim that host background has been controlled
away.

The remaining pre-submission hard task is recovery and audit of the 56-sample
TPM result against one common 48-MAG derep reference. Habitat enrichment, IS
hotspot, spatial heterogeneity, activity, and element-cycle coupling remain
unresolved until that abundance analysis is complete.

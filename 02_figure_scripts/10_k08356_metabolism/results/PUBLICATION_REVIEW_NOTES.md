# Publication-risk revision record

## Resolved in the current figure

1. The canonical count is explicitly fixed at 49 K08356 sequences from 48
   unique MAGs, with clade counts 3/7/19/20. Older mappings are audited and
   marked superseded.
2. Focal K08356 IDs are excluded from arsenite-oxidation marker scoring. This
   removes two circular `aioA` hits; no non-focal arxA/aioA marker remains.
3. Bubble prevalence is calculated after deduplicating host MAGs within each
   habitat-clade-feature cell.
4. Each panel shows sequence and unique-MAG sample sizes. One-MAG cells are
   marked as descriptive only.
5. Gray hatching denotes no clade member; hollow zero points denote members with
   no reconstructed feature.
6. S13 and S15 are assigned to ES using both the original sample metadata and
   the group-derep input manifest. Two stale classification labels are logged.
7. A 75% gene-set-coverage sensitivity figure accompanies the 50% main figure.
8. The calculation has been restored to the previous distinct-gene set
   algorithm documented in `MAG和Contig结果核验.pdf`; it does not use METABOLIC
   module-step coverage.
9. Flagellar Assembly is included as the previous 36-KO gene set.
10. Output names distinguish 49 focal sequences from 48 host MAGs.
11. Captions limit interpretation to genomic potential and state that gene sets
    can overlap and are not independent.

## Remaining limitation

Completeness and contamination are available and audited for all 48 MAGs from
the source metaWRAP tables. Completeness is uneven by habitat and clade, so the
50/10 filter alone does not eliminate this confounding. A validated family-level GTDB taxonomy table was not
found in the inspected group-dRep or METABOLIC paths. Therefore, a family- and
completeness-adjusted clade model is not reported. It should only be added after
the exact 48 MAG IDs are matched to a GTDB family table; CheckM lineage is too
coarse to substitute for that analysis.

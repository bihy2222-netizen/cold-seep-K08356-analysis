# GTDB family taxonomy audit

The current K08356 set contains 48 unique MAGs. The pre-existing all-derep
GTDB-Tk classification contains 47 of those exact IDs. The group-derep-only MAG
`R2111_N500_0-10_bin13` is classified separately and then merged by exact MAG ID.

Both sources use:

- GTDB-Tk 2.4.1
- GTDB release r226
- `classify_wf --skip_ani_screen`

The raw source tables are retained in `raw_gtdb/`. Run
`../audit_GTDB_family_adjustment.R` only after both the overall and added-MAG
summary tables are present. The script asserts an exact 48/48 match before
writing taxonomy, confounding, model-rank, and within-Rhodobacteraceae outputs.

Family is taken from the `f__` rank in the GTDB classification string. Blank
family assignments remain `Unclassified_family`; CheckM lineage is not used as
a substitute.

Family, habitat, and clade are strongly entangled in this data set. Any adjusted
test is descriptive and must pass the design-rank audit. The script also reports
a focused IS-only, Rhodobacteraceae-only Clade 3 versus Clade 4 sensitivity for
the locked arsenic/B12 features.

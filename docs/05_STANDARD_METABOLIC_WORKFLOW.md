# Standard METABOLIC-G workflow

The reproducible K08356 MAG workflow is documented in:

`01_metagenomic_workflow/04_metabolism/METABOLIC_group_derep49/README.md`

The key distinction is analytical scope:

- `all_50_10_dRep/.../bin_METABOLIC` is a completed 808-MAG background run.
- The original subset has 48 K08356 proteins from 47 unique MAGs.
- Habitat-wise dereplication yields 49 K08356 proteins from 48 unique MAGs and
  adds `R2111_N500_0-10_bin13`, so a separate run is required.

Use the K08356 hit table to select unique MAG-level protein inputs. Do not pass
the 49 individual K08356 proteins to METABOLIC-G; the tool requires each MAG's
complete protein set (`*.cds.faa`) to reconstruct metabolic modules and element
cycles. `S1_9-12_bin1` has two K08356 proteins and is included once, giving 48
MAG inputs.

The 48-MAG METABOLIC-G v4.0 run completed on 2026-08-01 in 35 minutes 28
seconds. All cardinality and completion checks passed; see the workflow
manifest directory for the exact inputs, parameters, and validation summary.

## Publication figure scoring

The four-habitat figure workflow is in
`02_figure_scripts/10_k08356_metabolism/`. It restores the previous
distinct-gene-set coverage formula, calculates prevalence from unique MAGs, and
generates 50% and 75% threshold figures. The displayed matrix contains 48 MAGs
by 40 features (1920 unique rows): 37 non-arsenic legacy gene sets plus separate
binary `arrA`, `arsC`, and `arxA/aioA` marker-carriage rows.

Run `audit_submission_closeout.R` after the plot and quality scripts. It asserts
the canonical 49-sequence/48-MAG set, 3/7/19/20 clade counts, exact ID agreement
among the current tree, neighborhood, classification, and metabolism artifacts,
and the 48 x 40 score-matrix uniqueness contract. It deliberately reports the
abundance artifact as unverified until a current 49-sequence TPM table or figure
is supplied.

The definition audit preserves all 510 old-algorithm records for reproducibility
but reports that 30 records (19 distinct KOs) are unavailable in the current KO
result catalog. Affected coverage values are therefore conservative. MAGs all
pass the 50/10 quality threshold, but no completeness- or GTDB-family-adjusted
model is claimed.

The catalog-compatible denominator sensitivity analysis must accompany the
legacy figure: 41 MAG-feature calls flip at the 50% threshold and 10 at 75%,
restricted to TCA, second-stage TCA, pyruvate oxidation, reductive TCA,
denitrification, and urea-cycle rows. The locked Clade 4-IS arsenic-marker and
cobinamide-to-cobalamin B12 results are unaffected.

`S1_9-12_bin1` is the only MAG represented in two clades. Its exclusion makes
the descriptive n=1 Clade 1-AS and Clade 3-AS cells empty; Clade 4-IS and all
other populated cells remain unchanged. See `MANUSCRIPT_LOCKED_WORDING.md` and
the supplementary sensitivity tables in the figure workflow's `results/`
directory.

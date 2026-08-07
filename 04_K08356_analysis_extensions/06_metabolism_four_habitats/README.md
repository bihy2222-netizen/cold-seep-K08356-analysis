# Four-habitat metabolism plot for K08356seq49_MAG48

This workflow reconstructs the four-habitat, four-clade metabolic bubble plot
from METABOLIC-G v4.0 results for 48 unique group-dereplicated MAGs. The set
contains 49 focal K08356 proteins because `S1_9-12_bin1` contributes one
Canonical AioA sequence and one IdrA phylogenetic sequence.

## Canonical version

- K08356 sequences: 49
- Unique host MAGs: 48
- Clade sequence counts: 3 / 7 / 19 / 20
- This version supersedes the older 48-sequence/47-MAG mappings. The exact
  history is recorded in `results/version_count_audit.tsv`.

## Locked status

- Metabolic-potential figure: verified.
- Abundance/TPM figure: not verified; the current 49-sequence folder contains no
  recoverable 56-sample TPM artifact.
- Completeness and GTDB-family adjustment: not completed.

Accordingly, this figure supports unadjusted genomic metabolic-potential
associations only, not activity, element-cycle coupling, or habitat-driven effects.

## Inputs

- `METABOLIC_group_derep49_result.xlsx`: completed METABOLIC-G result workbook.
- `MAG49_input_manifest.tsv`: authoritative MAG input list and checksums.
- `../four_clade_classification_49.tsv`: sequence-level clade assignment.
- `sample_group_corrected.csv`: original sample-to-habitat metadata.
- `KEGG_identifier_result/*.cds.result.txt`: per-MAG KO calls used for the
  previous distinct-gene coverage algorithm.
- `gene_set_definitions_previous_algorithm.tsv`: the 510 distinct
  source-sheet/module/gene rows underlying the previous 38 gene sets.
- `quality/raw_stats/*.stats`: source metaWRAP quality tables for the 48 MAGs.

The current clade table contains stale AS labels for `S13_0-2_bin3` and
`S15_8-10_bin20`. Both the original sample metadata and the group-derep input
manifest place them in ES. The plot uses that validated assignment and records
the two corrections in `results/habitat_assignment_QC.tsv`.

## Previous algorithm and safeguards

- For each MAG and gene set, coverage is `100 * distinct detected genes /
  distinct predefined genes`. Extra copies of a detected gene do not increase
  coverage. This is the same algorithm documented in the July 2026
  `MAG和Contig结果核验.pdf`.
- Bubble fill is mean gene-set coverage across unique host MAGs in a
  habitat-clade cell. Gene sets overlap and are not independent biological
  observations.
- Bubble size is true MAG prevalence, not sequence-membership prevalence.
- For continuity with the pre-specified legacy algorithm, the main plot uses a
  >=50% partial gene-set-reconstruction threshold. A >=75% threshold is generated
  as a sensitivity analysis. Neither threshold alone establishes complete pathway
  presence, and no mechanistic inference relies solely on the 50% cutoff.
- Arsenic is split into three binary per-MAG METABOLIC HMM-marker rows:
  respiratory arsenate reduction (`arrA`), arsenate detoxification/resistance
  (`arsC`), and arsenite oxidation (`arxA`/`aioA`). These are not multi-gene
  coverage scores. All 49 focal K08356 sequence IDs are excluded before scoring
  to avoid circular evidence.
- The three canonical sequences score 1192.6, 1177.1, and 762.0 against the same
  METABOLIC `aioA.hmm`. The original run used `-T 800`, so only the first two
  appeared in its `tblout`. The third is retained by independent sequence/tree
  evidence and the joint-screen `Combined_score` of 969.5 (>=640), not by
  METABOLIC `aioA.hmm` detection. Phylogenetic placement and METABOLIC detection
  are separate evidence systems.
- Gray hatching means that a habitat has no member of that clade. A small hollow
  point means clade members are present but the feature score is zero.
- Every habitat-clade panel reports both sequence and unique MAG counts. Cells
  with one MAG are marked as descriptive only.
- Flagellar Assembly is a partial 36-KO gene-set reconstruction and does not by
  itself establish complete motility.
- The figure describes genomic potential only. It does not establish activity,
  transcript abundance, metabolic flux, or ecological coupling.

## Run

```bash
export METABOLIC_RESULT_XLSX=/path/to/METABOLIC_group_derep49_result.xlsx
export KEGG_IDENTIFIER_DIR=/path/to/KEGG_identifier_result
Rscript plot_K08356seq49_MAG48_metabolism_four_habitats.R
Rscript audit_MAG48_quality.R
Rscript audit_submission_closeout.R
```

The plotting script validates 49 sequence memberships, 48 unique MAGs, clade
counts of 3/7/19/20, 48 KO result files, 510 gene-definition rows, 38 previous
gene sets, 48 x 40 = 1920 unique MAG-feature rows, and exact agreement among the
MAG sets before writing output. The 38-set catalog contains one legacy arsenate
set; the figure replaces it with three biologically distinct HMM-marker rows, so
40 features are displayed.

The 510 definition records comprise 508 KO records and two legacy EC identifiers:
`EC:1.20.4.1` is ArsC-type detoxification/resistance, whereas `EC:1.20.99.1`
corresponds to donor-dependent respiratory ArrA activity. Thirty definition
records representing 19 distinct KOs are absent from the current per-MAG KO
result catalog. They remain
in the denominator to reproduce the previous algorithm exactly, making affected
gene-set coverage values conservative but not uniformly so. A catalog-compatible
denominator sensitivity analysis produces 41 MAG-feature threshold flips at 50%
and 10 at 75%, concentrated in TCA, second-stage TCA, pyruvate oxidation,
reductive TCA, denitrification, and urea-cycle rows. These rows must not be used
for threshold-dependent mechanistic claims without the sensitivity table. The
Clade 4-IS arsenic-marker rows and target cobinamide-to-cobalamin B12 row are
unchanged.

`S1_9-12_bin1` contributes one sequence to Clade 1 and one to Clade 3. Excluding
this dual-copy MAG makes the two original n=1 Clade 1-AS and Clade 3-AS cells
empty; Clade 4-IS and all other populated cells are unchanged.

## Locked biological result

Among 20 Clade 4-IS MAGs, `arsC` occurs in 17/20 (85%), `arrA` in 2/20
(10%), and no independent `arxA`/`aioA` marker remains after focal K08356
exclusion. All 20/20 exceed the legacy 50% partial-reconstruction threshold for
`Cobalamin biosynthesis, cobinamide => cobalamin`, with 79.375% mean coverage.
This supports a widespread arsenic detoxification/resistance background, not
widespread arsenic respiration or independent arsenite oxidation.

## Main outputs

- `results/K08356seq49_MAG48_previous_gene_set_coverage50.svg`
- `results/K08356seq49_MAG48_previous_gene_set_coverage75_sensitivity.svg`
- `results/K08356seq49_MAG48_MAG_clade_feature_scores.tsv`
- `results/K08356seq49_MAG48_quality_audit.tsv`
- `results/arsenic_marker_hit_focal_exclusion_QC.tsv`
- `results/canonical_AioA_METABOLIC_HMM_score_audit.tsv`
- `results/gene_set_definition_ID_audit_summary.tsv`
- `results/gene_set_non_KO_and_unmatched_definitions.tsv`
- `results/gene_set_unavailable_KO_threshold_sensitivity_by_feature.tsv`
- `results/gene_set_unavailable_KO_MAG_threshold_flips.tsv`
- `results/dual_copy_MAG_exclusion_sensitivity.tsv`
- `results/dual_copy_MAG_exclusion_sensitivity_summary.tsv`
- `results/Clade4_IS_locked_arsenic_B12_results.tsv`
- `results/FINAL_LOCK_STATUS.tsv`
- `results/MANUSCRIPT_LOCKED_WORDING.md`
- `results/artifact_49_sequence_sync_audit.tsv`
- `results/confounding_adjustment_status.tsv`
- `results/PUBLICATION_REVIEW_NOTES.md`

The metaWRAP audit supplies completeness and contamination for all 48 MAGs.
Completeness is imbalanced across habitats (median AS 60.79%, ES 71.22%, IS
88.25%, NS 79.02%) and clades, so passing the 50/10 filter does not remove
quality confounding.
No validated GTDB family-level host table was found in the inspected group-dRep
or METABOLIC paths, so a model such as `module ~ clade + family + habitat +
completeness` is not reported. CheckM lineage is retained as coarse QC metadata,
not substituted for family-level taxonomy.

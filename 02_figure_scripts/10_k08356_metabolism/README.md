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

## Inputs

- Completed METABOLIC-G result workbook supplied through
  `METABOLIC_RESULT_XLSX`.
- `01_metagenomic_workflow/04_metabolism/METABOLIC_group_derep49/manifests/MAG49_input_manifest.tsv`:
  authoritative MAG input list and checksums.
- `input/four_clade_classification_49.tsv`: sequence-level clade assignment.
- `sample_group_corrected.csv`: original sample-to-habitat metadata.
- Per-MAG `KEGG_identifier_result/*.cds.result.txt` files supplied through
  `KEGG_IDENTIFIER_DIR` for the previous distinct-gene coverage algorithm.
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
- The main plot uses a 50% gene-set-coverage threshold; a 75% sensitivity
  plot is generated in parallel.
- Arsenate reduction and arsenite oxidation use METABOLIC arsenic HMM hits.
  All 49 focal K08356 sequence IDs are excluded before scoring to avoid circular
  evidence. Two focal `aioA` hits are excluded and no non-focal arxA/aioA hit
  remains in this run.
- Gray hatching means that a habitat has no member of that clade. A small hollow
  point means clade members are present but the feature score is zero.
- Every habitat-clade panel reports both sequence and unique MAG counts. Cells
  with one MAG are marked as descriptive only.
- Flagellar Assembly is one of the previous gene sets and contains 36 KOs.
- The figure describes genomic potential only. It does not establish activity,
  transcript abundance, metabolic flux, or ecological coupling.

## Run

```bash
export METABOLIC_RESULT_XLSX=/path/to/METABOLIC_result.xlsx
export KEGG_IDENTIFIER_DIR=/path/to/KEGG_identifier_result
export K08356_CLADE_TSV="$PWD/02_figure_scripts/10_k08356_metabolism/input/four_clade_classification_49.tsv"
export MAG49_MANIFEST_TSV="$PWD/01_metagenomic_workflow/04_metabolism/METABOLIC_group_derep49/manifests/MAG49_input_manifest.tsv"
export METABOLIC_PLOT_OUT_DIR="$PWD/02_figure_scripts/10_k08356_metabolism/results"

Rscript 02_figure_scripts/10_k08356_metabolism/plot_K08356seq49_MAG48_metabolism_four_habitats.R
Rscript 02_figure_scripts/10_k08356_metabolism/audit_MAG48_quality.R
```

The plotting script validates 49 sequence memberships, 48 unique MAGs, clade
counts of 3/7/19/20, 48 KO result files, 510 gene-definition rows, 38 previous
gene sets, and exact agreement among the MAG sets before writing output. The
38-set catalog contains one legacy arsenate set; the figure replaces it with
separate arsenate-reduction and arsenite-oxidation HMM-component scores, so 39
features are displayed.

The METABOLIC workbook and per-MAG KO result files are not committed; they are
server-side/generated inputs. The fixed gene-set definitions, metaWRAP quality
source tables, editable SVGs, and plot-level source tables are committed.

## Main outputs

- `results/K08356seq49_MAG48_previous_gene_set_coverage50.svg`
- `results/K08356seq49_MAG48_previous_gene_set_coverage75_sensitivity.svg`
- `results/K08356seq49_MAG48_MAG_clade_feature_scores.tsv`
- `results/K08356seq49_MAG48_quality_audit.tsv`
- `results/arsenic_marker_hit_focal_exclusion_QC.tsv`
- `results/PUBLICATION_REVIEW_NOTES.md`

The metaWRAP audit supplies completeness and contamination for all 48 MAGs.
Completeness is imbalanced across habitats (median AS 60.79%, ES 71.22%, IS
88.25%, NS 79.02%) and clades, so passing the 50/10 filter does not remove
quality confounding.
No validated GTDB family-level host table was found in the inspected group-dRep
or METABOLIC paths, so a model such as `module ~ clade + family + habitat +
completeness` is not reported. CheckM lineage is retained as coarse QC metadata,
not substituted for family-level taxonomy.

# K08356 MAG metabolism across four habitats and four clades

This workflow plots the completed METABOLIC-G v4.0 analysis of the habitat-wise
group-dereplicated K08356 MAG set. It distinguishes 49 K08356 sequence
memberships from 48 unique host MAGs and validates clade counts of 3/7/19/20.

## Inputs

1. The completed 48-MAG `METABOLIC_result.xlsx` workbook.
2. `input/four_clade_classification_49.tsv` for sequence-level clades.
3. The authoritative MAG input manifest at
   `01_metagenomic_workflow/04_metabolism/METABOLIC_group_derep49/manifests/MAG49_input_manifest.tsv`.

The classification table contains two stale habitat labels: `S13_0-2_bin3` and
`S15_8-10_bin20` are labelled AS there, but their group-derep source is ES. The
plot uses the input manifest as habitat authority and records both corrections
in `results/habitat_assignment_QC.tsv`.

## Run

```bash
export METABOLIC_RESULT_XLSX=/path/to/METABOLIC_result.xlsx
export K08356_CLADE_TSV="$PWD/02_figure_scripts/10_k08356_metabolism/input/four_clade_classification_49.tsv"
export MAG49_MANIFEST_TSV="$PWD/01_metagenomic_workflow/04_metabolism/METABOLIC_group_derep49/manifests/MAG49_input_manifest.tsv"
export METABOLIC_PLOT_OUT_DIR="$PWD/02_figure_scripts/10_k08356_metabolism/results"
Rscript 02_figure_scripts/10_k08356_metabolism/plot_MAG49_metabolism_four_habitats.R
```

KEGG modules are scored by detected `KEGGModuleStepHit` steps. Arsenate
reduction and arsenite oxidation are scored from METABOLIC `FunctionHit`
presence. Bubble fill is the mean reconstruction score; bubble size is the
percentage of K08356 sequence memberships whose host MAG score is at least
50%.

The editable SVG and the exact plotted source tables are committed under
`results/`. The METABOLIC workbook, PNG, and PDF stay outside Git because they
are generated or server-side artifacts.

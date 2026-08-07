# K08356 analysis extensions

This directory records the K08356-specific analyses developed after the core
metagenome pipeline was assembled. Scripts are grouped by analytical role while
retaining their original subdirectory names where those names encode figure or
revision provenance.

## Directory map

- `01_MAG_abundance_figures_ml/`: gene-level and MAG-level abundance reports,
  tree audits, gene-island plots, nutrient association, machine learning,
  Sankey plots, metabolic-module analyses, and Clade 4 ecological summaries.
- `02_hmm_tree_workflow/`: 48/49-sequence HMM comparison and standardized small
  AioA/IdrA/DMSOR tree workflow.
- `03_gene_neighborhood_workflow/`: dereplicated-MAG neighborhood extraction and
  synteny summaries.
- `04_MAG49_integration/`: integration of HMM scores and MAG neighborhood tables.
- `05_TPM_figures/`: contig- and MAG-level K08356 TPM figure revisions.
- `06_metabolism_four_habitats/`: four-habitat metabolism and quality audits.
- `07_tree_visualization/`: four-clade tree and tree-aligned gene-island helpers.

## Reproducibility status

Scripts in the numbered core workflow are configuration-driven and recommended
for reuse. Scripts in this extension directory are analysis snapshots retained
for full provenance. Some snapshots contain historical absolute input paths;
set those paths to the current project layout before running. Their statistical
logic and figure settings were intentionally preserved.

Large data, generated figures, R package libraries, HMM models, and intermediate
results are excluded from Git. No credentials or access tokens are included.

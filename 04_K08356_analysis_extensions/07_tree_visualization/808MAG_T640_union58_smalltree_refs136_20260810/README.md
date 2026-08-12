# 808 MAG three-model T640 union58 phylogeny and gene islands

This directory contains the accepted small-tree analysis for the union of 58 high-scoring candidates from three HMM models across 808 dereplicated MAG protein catalogs. The candidates come from 56 MAGs and were analyzed together with 136 references (194 sequences total).

The 58 proteins are a candidate union, not 58 functionally confirmed IdrA proteins. Tree placement, HMM scores, and gene-neighborhood evidence are integrated into four conservative working classes:

- Canonical AioA clade: 3
- Unknown DMSOR clade: 8
- IdrA-associated clade: 27
- IdrA clade: 20

The 48 candidates overlapping the prior strict-QC evidence set retain their prior classification. The 10 new T640 candidates remain provisional: eight are assigned to the IdrA-associated working class and two to Unknown DMSOR based on tree placement and available neighborhood evidence.

## Key outputs

- `02_tree/MAG_T640_union58_plus_refs136.treefile`: final ML tree with SH-aLRT/UFBoot labels
- `figures/MAG_T640_union58_four_clades_gene_islands.svg`: editable tree plus aligned gene-island panel
- `figures/MAG_T640_union58_four_clades_gene_islands.png`: high-resolution preview
- `04_tables/candidate58_four_clade_integrated_evidence.tsv`: per-candidate classification evidence
- `04_neighborhood/union58_gene_island_manifest.tsv`: gene-island plotting manifest

IQ-TREE selected `LG+R6` by BIC. The final run used 1,000 SH-like aLRT replicates and 1,000 ultrafast bootstrap replicates.

## Reproduction

The tree command is recorded in `scripts/05_run_MAG_T640_union58_smalltree.sh`. Neighborhood extraction, evidence integration, and SVG assembly are implemented in scripts `08` through `10`. Scripts `09` and `10` resolve inputs relative to this directory.

Functional assignments remain subject to further gene-neighborhood, domain, and manual biological validation.

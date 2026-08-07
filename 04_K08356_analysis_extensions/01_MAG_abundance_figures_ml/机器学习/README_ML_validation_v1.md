# K08356 ML-assisted validation v1

This is a small local test focused on MAG-level K08356 gene-neighborhood features.

## Source
- ORF coordinates and annotations were derived from `/Users/catherine/Downloads/MAG-bin-tpm/K08356_MAG_gene_island_arrows_20260708/K08356_four_branch_gene_neighborhood_coordinates.csv`.
- The underlying sample-level ORF/GFF source path provided by the user is `/home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/CDS/cat_result/cat_others`.
- No new files were written to the supercomputer source directory.

## Feature design
- Each row in the feature matrix is one target K08356 neighborhood.
- Predictors use only neighboring ORFs; the target K08356 ORF itself is excluded from model features to avoid label leakage.
- Both count and binary features were generated for functional classes used in the gene-island figure.

## First-pass modeling
- PCA and heatmap are descriptive.
- PERMANOVA tests whether gene-neighborhood feature composition differs among branches.
- Leave-one-out nearest-centroid classification is used as a lightweight validation suitable for small n.
- Decision-tree feature importance is exploratory only.

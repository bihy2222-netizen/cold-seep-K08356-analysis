# Protein Sequence-Space Validation v1

This run validates K08356-related DMSOR branch separation using a reproducible local sequence-feature embedding.

## Inputs
- FASTA: /Users/catherine/Downloads/aioA 蛋白序列建树/tree11_no_uncertain_DryadAioAIdrAfinal/all_proteins.faa
- tree11 contig clade table: /Users/catherine/Downloads/aioA 蛋白序列建树/tree11_no_uncertain_DryadAioAIdrAfinal/tree11_clade_detail.csv
- MAG branch label table for MAG-derived flags: /Users/catherine/Downloads/MAG-bin-tpm/机器学习/00_input/K08356_branch_label_table.csv

## Method
- Protein sequences were represented by amino-acid composition plus dipeptide frequency features.
- Features were scaled before PCA, PCoA/MDS, centroid distance, silhouette, and PERMANOVA analyses.
- The main quantitative summaries use the 230 ColdSeep contig K08356 targets with tree11 clade labels.
- Reference sequences are retained in PCA/PCoA/UMAP plots as contextual anchors.
- This is a local sequence-space validation, not an ESM-2 protein language model run.
- If the R `umap` package is available to Rscript, UMAP is generated automatically. Otherwise the UMAP figure states that PCoA/MDS is the fallback.

## Interpretation boundary
- This analysis can support whether branches occupy distinguishable protein sequence-feature space.
- It cannot directly prove substrate specificity, catalytic activity, or enzyme function.
- The current run uses the tree11 contig-level K08356 target set and clade labels.
- The 48 MAG-derived K08356 labels are used only to flag which contig targets are also MAG-derived.
- ESM-2/ProtT5 can be added later as an upgraded protein language model embedding, but was not run in this local v1.

## Key outputs
- protein_sequence_feature_embedding_matrix.csv
- protein_sequence_PCA_scores.csv
- protein_sequence_UMAP_scores.csv if available
- branch_centroid_distance_matrix.csv
- silhouette_summary_by_branch.csv
- PERMANOVA_embedding_distance_by_branch.csv

UMAP available to Rscript in this run: TRUE

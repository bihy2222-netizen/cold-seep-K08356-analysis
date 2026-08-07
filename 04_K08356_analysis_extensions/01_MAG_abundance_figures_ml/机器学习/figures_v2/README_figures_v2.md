# figures_v2 notes

Generated revised gene-neighborhood ML figures:
- A_feature_prevalence_heatmap.pdf/svg
- B_Jaccard_PCoA_PERMANOVA.pdf/svg
- C_RF_feature_importance_top10.pdf/svg
- D_RF_branchwise_LOOCV_accuracy.pdf/svg
- E_RF_label_permutation_test.pdf/svg
- Figure_gene_neighborhood_ML_four_panel.pdf/svg

PERMANOVA: R2 = 0.238, F = 4.570, p = 0.001.
RF label permutation test: observed balanced accuracy = 0.599; random mean = 0.236; p = 0.005; permutations = 200.

Interpretation: use Random Forest feature importance as exploratory support, not causal proof.

Heatmap stars indicate branch-associated binary feature differences based on Fisher's exact tests with Benjamini-Hochberg FDR correction (* FDR < 0.05; ** FDR < 0.01; *** FDR < 0.001). Stars are shown on the feature row because the test is feature-level across branches.

The RF label-permutation test used 200 label permutations and should be treated as an exploratory small-sample validation. It supports that the observed balanced accuracy is higher than expected under random branch labels, but it does not imply causal enzyme-function prediction.

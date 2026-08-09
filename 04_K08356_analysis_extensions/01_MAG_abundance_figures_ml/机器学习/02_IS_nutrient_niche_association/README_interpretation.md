# IS internal nutrient-DMSOR niche association analysis

This analysis links nutrient/organic matter gradients with K08356-related DMSOR gene- and MAG-level abundance within IS samples only.

## Boundary
- This is an IS-internal microenvironment association analysis.
- It must not be interpreted as a four-habitat environmental driver analysis because AS/ES/NS nutrient data are not available here.

## Data
- Nutrient table: `/Users/catherine/Downloads/营养盐热图绘制/varechem.csv`.
- MAG branch TPM: `K08356_branch_host_MAG_TPM_by_sample.csv`.
- Gene-level K08356 TPM: `K08356_gene_total_sample_long.csv`.
- Sample matching: 21 nutrient samples; 21 matched to gene TPM; 21 matched to MAG TPM.

## Derived nutrient variables
- DIN = NH4_mgkg + NO3_mgkg.
- NH4_NO3_ratio = NH4_mgkg / NO3_mgkg.
- TOC_TOP_ratio = TOC_mgkg / TOP_mgkg.

## Abundance transformation
- All K08356 abundance metrics were transformed as log10(TPM + 1) for correlation plots.

## Figure interpretation
- `nutrient_layer_boxplots`: IS nutrient and organic matter variation among sediment layers.
- `nutrient_K08356_spearman_heatmap`: Spearman rho between nutrient variables and log-transformed K08356 abundance metrics.
- `selected_scatterplots`: automatically selected nutrient-DMSOR pairs with relatively high |rho| and lower p/FDR, prioritizing non-canonical and total K08356 metrics.
- `exploratory_RDA_nutrient_branch_abundance`: exploratory constrained ordination using NH4, NO3, TOP and TOC only.

## Recommended wording
These results should be described as nutrient-DMSOR niche association within IS samples, supporting micro-scale environmental heterogeneity potentially related to non-canonical DMSOR-bearing MAG distributions. They should not be described as causal evidence or as drivers of four-habitat differences.

# K08356 host taxonomy Sankey

Input: /Users/catherine/Downloads/营养盐热图绘制/results_revision/02_K08356_host_origin/K08356_sequence_MAG_taxonomy_master.tsv

Taxonomic ranks shown: Phylum -> Class -> Order -> Family -> Genus -> Species.
Weights are K08356 sequence counts from the 48-row K08356 MAG host taxonomy master table.
Species was parsed from the GTDB-style Classification field. Empty species annotations were written as s__unclassified.
For readability, species observed in fewer than 2 sequences were collapsed to s__rare species (<2 seq) in the plotted figure.
The full uncollapsed path count table is retained in results/K08356_host_taxonomy_path_counts_full_species.csv.

Statistical note:
- Clade x Rhodobacteraceae-vs-other Fisher exact test: P = 9.054e-06; Cramer's V = 0.691.
- Clade 4 contained 17 of 20 K08356 loci in Rhodobacteraceae-assigned MAGs; standardized residual r = 4.6.
- These statistics describe enrichment within the 48 candidate-locus dataset and do not represent taxon abundance in the complete 808-MAG background.

Outputs:
- figures/K08356_host_taxonomy_sankey_phylum_to_species.pdf
- figures/K08356_host_taxonomy_sankey_phylum_to_species.svg
- figures/K08356_host_taxonomy_sankey_phylum_to_species.png
- figures/K08356_host_taxonomy_sankey_phylum_to_species.eps
- figures/K08356_host_taxonomy_sankey_clade_to_species.pdf
- figures/K08356_host_taxonomy_sankey_clade_to_species.svg
- figures/K08356_host_taxonomy_sankey_clade_to_species.png
- figures/K08356_host_taxonomy_sankey_clade_to_species.eps
- figures/K08356_host_taxonomy_sankey_clade_to_species_with_n.pdf
- figures/K08356_host_taxonomy_sankey_clade_to_species_with_n.svg
- figures/K08356_host_taxonomy_sankey_clade_to_species_with_n.png
- figures/K08356_host_taxonomy_sankey_clade_to_species_with_n.eps

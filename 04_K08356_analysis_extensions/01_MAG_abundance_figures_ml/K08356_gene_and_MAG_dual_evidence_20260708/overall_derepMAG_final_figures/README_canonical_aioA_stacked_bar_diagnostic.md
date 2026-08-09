# canonical aioA-associated stacked-bar diagnostic

The branch stacked bar is not an origin-only occurrence plot. It shows read-recruited MAG TPM across all samples.

For canonical aioA-associated K08356, the branch map contains three derepMAG entries:
- S1_9-12_bin1 | source habitat: AS | original gene id: S1_9-12_bin1-k141_637215_2
- SY366YW-4-8_bin7 | source habitat: IS | original gene id: SY366YW-4-8_bin7-k141_14765_13
- R2111_S300_0-10_bin19 | source habitat: NS | original gene id: R2111_S300_0-10_bin19-k141_825759_6

In IS samples, nonzero canonical-associated TPM is contributed mainly by:
- SY366YW-4-8_bin7 | source habitat: IS | nonzero IS samples: 6 | sum IS TPM: 85318.902 | max IS sample: SY366YW-4-8 | max TPM: 81186.07
- R2111_S300_0-10_bin19 | source habitat: NS | nonzero IS samples: 14 | sum IS TPM: 26434.036 | max IS sample: SY366YB-8-12 | max TPM: 4875.858

Interpretation: the IS-origin canonical MAG is present mainly at SY366YW-4-8, but the NS-origin R2111_S300 canonical derepMAG also recruits reads in multiple IS samples. Therefore orange segments in several IS bars are abundance/recruitment signals from the MAG TPM matrix, not additional IS-origin canonical bins.

Recommended wording: `canonical aioA-associated host MAG abundance was not significantly different among habitats (Kruskal-Wallis p = 0.877); although only one canonical-associated MAG was recovered from IS, low-level read recruitment to canonical-associated derepMAG representatives was observed across multiple IS samples.`

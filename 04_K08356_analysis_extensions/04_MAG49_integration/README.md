# K08356 MAG49 neighborhood integration

This directory preserves the original 48-MAG tables under `source_48/` and adds
the group-derep-specific MAG candidate:

`R2111_N500_0-10_bin13-k141_4404301_3`

## Evidence added for the new MAG

- Combined HMM score: 780.4 (published 640 comparison passed)
- IdrA HMM score: 874.2
- AioA HMM score: 487.3
- Small-tree sister: `IS|SY368YW-8-12_bin16-k141_314641_4`
- Sister-node support: 100/100
- Nearest non-MAG reference: `ISMEJFig5|Unknown_clade|MBT97654.1|Dehalococcoidia_bacterium`
- Nearest-reference tree distance: 1.138459
- Reliable neighborhood: K08355/AioB followed by K08356/AioA-like on the same strand
- Joint reference search: accepted canonical AioB hit; no accepted P-like protein
- Strict synteny QC: failed because two distinct P-like genes are absent
- Final conservative class: `unknown/uncertain DMSOR` (Medium confidence)

## Main outputs

- `outputs/K08356_MAG49_integrated_evidence.xlsx`
- `outputs/updated_final_evidence_table_49_strict_QC.tsv`
- `outputs/gene_neighborhood_49_summary.tsv`
- `outputs/gene_neighborhood_49_long.tsv`
- `outputs/strict_synteny_qc_49_reclassified.tsv`
- `outputs/HMM_scores_49.tsv`
- `outputs/joint_reference_accepted_hits_long_49.tsv`
- `outputs/new_MAG_joint_reference_accepted_hits.tsv`
- `outputs/validation_summary.txt`

## Rebuild

Run `build_mag49_tables.py` first, then `build_mag49_workbook.mjs`.

The original 48-MAG files are not overwritten.

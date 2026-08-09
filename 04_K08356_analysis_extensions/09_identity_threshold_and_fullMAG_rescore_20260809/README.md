# Four-clade identity and full-MAG HMM discovery

This workflow separates two questions that must not be conflated:

1. Alignment-derived pairwise amino-acid identity among the 49 curated proteins
   from 48 unique MAGs.
2. Profile-HMM-driven discovery across all proteins from 808 dereplicated MAGs.

K08356 is not used as the discovery gate. It is retained as an independent
conventional-annotation comparator for later recall and mixed-function audits.

## Curated-core analysis

The fixed local groups are named `canonical AioA-associated` (n=3),
`synteny-supported strict DIRM-like IdrA` (n=20), `partial IdrA-associated`
(n=19), and `AioA-like or unresolved DMSOR` (n=7). The last two are analytical
categories rather than experimentally validated enzyme assignments. The primary
identity is EMBOSS Needle full-global
identity including one-sided gaps. Untrimmed MAFFT and MMseqs2 local identity
with bidirectional coverage provide independent sensitivity analyses.

The old 48-protein set and latest 49-protein set differ by one group-derep
candidate. The latest unit of analysis is 49 proteins from 48 MAGs; one MAG
contains two phylogenetically distinct candidates.

## Full-MAG discovery

The two-stage HMM design avoids an unfiltered `--max` scan of 1.89 million
proteins:

1. Normal accelerated HMMER discovery at `-T 100 --domT 20` for IdrA, AioA,
   and combined profiles.
2. Union extraction followed by exact three-profile rescoring with
   `--max -T 0 --domT 0`.

Existing T640 outputs are preserved. Candidate function remains unresolved
until phylogenetic and gene-neighborhood validation is complete.

## Scripts

- `00_run_core49_external_aligners.sh`: MAFFT, trimAl, MMseqs2, and EMBOSS.
- `01_compute_core49_pairwise_identity.py`: QC and global/local pairwise tables.
- `02_integrate_mafft_mmseqs_identity.py`: MAFFT and MMseqs2 integration.
- `03_summarize_integrated_identity_metrics.py`: clade summaries and outliers.
- `04_parse_emboss_needleall.py`: canonical Needle identity tables and an
  implementation-level comparison against the Biopython global alignment.
- `04_run_full808_two_stage_HMM_rescore.sh`: T100 discovery and exact rescoring.
- `01_summarize_union58_rescore.py`: urgent T640-union score/coverage matrix.
- `02_summarize_T100_union_rescore.py`: 720-candidate union and K08356 overlap.

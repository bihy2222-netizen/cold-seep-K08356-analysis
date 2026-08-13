# Metatranscriptome helper scripts

## `fetch_tree_dmsor_cds.py`

Recover nucleotide CDS for DMSOR-family protein references used in the protein tree.
This is needed because RNA read mapping must use nucleotide references, not protein
HMMs or protein FASTA files.

Example:

```bash
ROOT=/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra/maggie_expression_validation_20260812
BROAD=/home/ps/ps1/data/bihongyu/cold_seep/ASproteintree/aioA-like_tree_results_20260629_final/aioA-like-contigstree-5/broad_idra_dmsor_refs.large_subunit.filtered.faa
ARR=/home/ps/ps1/data/bihongyu/cold_seep/ASproteintree/aioA-like_tree_results_20260629_final/aioA-like-contigstree-5/ColdSeep_K08356_plus_Fig5_plus_broad_DMSOR_refs_refDedup_shortLabels_allNonstandardToX.faa

python3 scripts/metatranscriptome/fetch_tree_dmsor_cds.py   --broad-faa "$BROAD"   --arr-faa "$ARR"   --outdir "$ROOT"   --per-class 6
```

Outputs include a fetch audit table, a manifest of accepted CDS references, and
FASTA files for CDS/proteins. Keep failed WP/SwissProt records in the audit file;
they document accession classes that cannot be directly converted to CDS via
NCBI `fasta_cds_na`.

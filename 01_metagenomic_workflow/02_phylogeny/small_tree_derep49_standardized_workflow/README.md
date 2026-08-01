# Small Tree Derep49 K08356 Workflow

This folder standardizes the small AioA/IdrA/Unknown derep49 workflow for the tree:

`tree_small_AioA_IdrA_Unknown_group_derep49_refs136`

## What Is Included

- `scripts/01_run_hmm_supplement.sh`: original remote HMM `-T 640` workflow.
- `scripts/02_merge_hmm_scores.py`: original HMM table merge script.
- `scripts/extract_gene_neighborhood_48.py`: original neighborhood extraction helper.
- `scripts/make_itol_order_gene_island_outputs.py`: standardized small-tree plotting and summary-table script.
- `input/`: iTOL SVG, small-tree Newick, tip metadata, strict-QC table, and gene-arrow manifest used to recover top-to-bottom order and regenerate local outputs.
- `output/remote_T640_20260801/`: archived remote rerun of the HMM threshold workflow.
- `output/figures/`: beautified tree SVG and iTOL-order gene-island arrow SVG.
- `output/tables/`: iTOL row order, T640 check table, and IdrA-PP/gene-island prediction table.

## Remote 640 Rerun

The 640 threshold workflow was rerun on `2026-08-01` on:

`/home/ps/ps1/data/bihongyu/cold_seep`

Command pattern:

```bash
cd /home/ps/ps1/data/bihongyu/cold_seep/HMM_supplement
export PATH=/home/ps/anaconda3/envs/hmmer/bin:$PATH
OUT=result_small_tree_derep49_T640_20260801
CAND=../illu/binning/K08356_four_group_annotation/strict_kofam_best/all_groups.K08356.faa
ALLMAG=../illu/binning/all_50_10_dRep/dereplicated_genomes/bin_DRAM/annotation/genes.faa
HMMDIR=../Neighborhood_Analyses/hmm
bash 01_run_hmm_supplement.sh "$CAND" "$ALLMAG" "$HMMDIR" "$OUT" 8
```

Remote run summary:

- Small-tree candidate FASTA: 49 sequences.
- all MAG protein database: 129,785 sequences.
- combined candidate HMM hits: 49/49.
- idrA candidate HMM hits: 49/49.
- aioA candidate HMM hits: 49/49.
- all MAG `combined_iriA_aioA.hmm -T 640` hits: 48.
- new T640 hits outside the 49 candidate set by normalized ID: 0.

Important note: the extra NS sequence `NS|R2111_N500_0-10_bin13-k141_4404301_3` passes the candidate HMM search with combined score `780.4`, idrA score `874.2`, and aioA score `487.3`. It is absent from the all MAG T640 table because the all MAG `genes.faa` table does not contain that sequence under the normalized ID used by the merge script.

## Gene-Island Prediction Rule

The standardized prediction table uses the iTOL SVG row order and the existing strict-QC neighborhood/gene-arrow evidence:

- `IdrA-PP gene island: strong`: focal A plus B-related and at least two P-like CDS arrows.
- `IdrA-PP gene island: partial`: focal A plus B-related and one P-like CDS arrow.
- `B-linked, P-like missing`: focal A plus B-related but no P-like CDS arrow.
- `canonical Aio-like`: focal A plus canonical AioB.
- `unresolved/truncated`: insufficient local gene-island evidence.

The P-like labels are conservative positional labels from joint-reference search, not definitive IdrP1/IdrP2 names.

## Regenerate Local Outputs

From this repository root:

```bash
bash 01_metagenomic_workflow/02_phylogeny/small_tree_derep49_standardized_workflow/run_workflow.sh
```

Main outputs:

- `output/figures/c0J-IDrhpB2odQ8jRmUo-Q_beautified.svg`
- `output/figures/small_tree_gene_island_arrows_iTOL_svg_order.svg`
- `output/tables/small_tree_T640_threshold_check.tsv`
- `output/tables/small_tree_idrA_PP_gene_island_predictions.tsv`

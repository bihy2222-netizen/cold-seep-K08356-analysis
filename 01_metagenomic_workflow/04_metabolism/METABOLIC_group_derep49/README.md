# METABOLIC-G workflow for K08356 MAG sets

This workflow separates two analyses that were previously easy to conflate:

1. The existing overall-derep METABOLIC-G run used 808 MAG protein files. Its
   K08356 subset contains 48 K08356 proteins from 47 unique MAGs.
2. Habitat-wise dereplication recovered 49 K08356 proteins from 48 unique MAGs.
   This set includes `R2111_N500_0-10_bin13`, which is absent from the
   overall-derep input by design and therefore needs a separate run.

`S1_9-12_bin1` contains two K08356 proteins. METABOLIC-G is a MAG-level tool, so
its complete protein set is supplied once. The standardized run therefore has
48 input `*.cds.faa` files, not 49 duplicated MAG inputs.

The scripts do not copy genome proteins or METABOLIC databases into Git. They
create a checked symlink input directory, preserve a manifest, run METABOLIC-G,
and validate output cardinalities.

## Remote paths used for the 2026-08-01 run

```text
K08356 hit table:
/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/K08356_four_group_annotation/strict_kofam_best/all_groups.K08356_hits.tsv

Group-derep roots:
/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/{IS,AS,ES,NS}/{GROUP}_50_10_dRep/dereplicated_genomes

METABOLIC-G v4.0:
/home/ps/Research/data/working_dir/wangpudi/METABOLIC/METABOLIC-G.pl

Conda environment:
METABOLIC_v4.0
```

## Run

```bash
WORK=/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/K08356_four_group_annotation/METABOLIC_group_derep49_20260801
BINNING=/home/ps/ps1/data/bihongyu/cold_seep/illu/binning
HITS="$BINNING/K08356_four_group_annotation/strict_kofam_best/all_groups.K08356_hits.tsv"

python3 scripts/prepare_group_derep_input.py \
  --hits "$HITS" \
  --binning-root "$BINNING" \
  --input-dir "$WORK/input_faa" \
  --manifest "$WORK/manifests/MAG49_input_manifest.tsv" \
  --expected-mag-count 48

bash scripts/run_metabolic_g.sh \
  --input-dir "$WORK/input_faa" \
  --output-dir "$WORK/output" \
  --log-dir "$WORK/logs" \
  --threads 100 \
  --module-cutoff 0.75 \
  --kofam-db full \
  --conda-env METABOLIC_v4.0 \
  --metabolic-script /home/ps/Research/data/working_dir/wangpudi/METABOLIC/METABOLIC-G.pl

python3 scripts/validate_metabolic_output.py \
  --input-dir "$WORK/input_faa" \
  --output-dir "$WORK/output" \
  --expected-mag-count 48 \
  --summary "$WORK/manifests/MAG49_output_validation.tsv" \
  --strict
```

## Fixed parameters

- METABOLIC-G version: 4.0
- Input mode: translated MAG proteins (`-in`)
- Threads: 100
- KEGG module cutoff: 0.75
- KOfam database: full
- Prodigal: not rerun because group-derep `*.cds.faa` files are supplied

## Completed run

The standardized run started at `2026-08-01T20:45:57+08:00` and completed at
`2026-08-01T21:21:27+08:00` in 35 minutes 28 seconds. Its validated output is:

```text
/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/K08356_four_group_annotation/METABOLIC_group_derep49_20260801/output
```

All strict checks passed: 48 MAG inputs, 96 KEGG identifier files, 192
nutrient-cycling PDFs, six per-category spreadsheets, a non-empty final
workbook, and the METABOLIC-G completion marker. Exact run parameters and
validation results are preserved in `manifests/MAG49_run_config.tsv` and
`manifests/MAG49_output_validation.tsv`.

## Expected output cardinalities

The completed 808-MAG background run produced two KEGG identifier files and
four nutrient-cycling PDFs per input MAG. The validator checks the same
invariants for the 49-protein/48-MAG run, together with the final workbook and
completion marker.

Only scripts, manifests, audit records, and concise validation summaries should
be committed. Full logs, METABOLIC intermediate files, databases, PDFs, and
large workbooks stay on the analysis server.

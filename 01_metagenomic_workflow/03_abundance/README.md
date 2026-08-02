# Unified MAG48 abundance workflow

All 56 metagenomic samples must be mapped to the same exact set of 48
K08356-bearing MAGs. Habitat-specific derep genome directories are not valid
references for cross-habitat abundance comparisons because their genome
composition differs.

## Prepare and validate

```bash
RUN_COVERM=0 bash 02_run_unified_MAG48_coverm_tpm.sh
```

This creates one symlink-only `reference_MAG48/`, a 56-sample read/habitat
manifest, an exact 48-MAG checksum manifest, and a CoverM settings manifest.
Preparation stops unless all 56 sample pairs and all 48 MAG FASTAs exist.
The validated 2026-08-02 preparation manifests are retained in `manifests/`.
Eight legacy sample-ID typos corrected against the actual CoverM commands and
FASTQ paths are recorded in `sample_ID_correction_audit.tsv`; no habitat label
changes.

## Run

```bash
RUN_COVERM=1 THREADS=100 bash 02_run_unified_MAG48_coverm_tpm.sh
```

Every sample uses the same settings: TPM, minimum covered fraction 0.1, read
identity 95%, aligned-read fraction 75%, and 0.1/0.9 end trimming. Existing
non-empty sample outputs are skipped, so the command can be resumed.

After all samples finish, `03_audit_unified_MAG48_tpm.R` asserts an exact
56 x 48 matrix (2,688 unique sample-MAG rows) and writes long, matrix, clade, and
dual-clade sensitivity tables.

`S1_9-12_bin1` is handled three ways:

- non-exclusive Clade 1 and Clade 3 summaries, which must not be summed;
- one exclusive `Multi-clade host` category;
- a sensitivity summary excluding the MAG.

The old `01_run_group_derep_coverm_tpm.sh` name is retained only as a compatible
entry point and forwards to this unified workflow. Existing habitat-specific
TPM outputs are process records only and must not enter cross-habitat statistics.

Statistical tests and ecological plots are deferred until the exact matrix is
complete. The 56 samples, not the 48 MAGs, are the environmental replicates.

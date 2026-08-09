# Script manifest and reproducibility status

The repository currently contains 151 shell, Python, and R scripts. The
machine-generated SHA-256 inventory is stored in
`docs/05_SCRIPT_MANIFEST_SHA256.tsv`.

## Status classes

- `01_metagenomic_workflow/`: current reusable workflow. New primary scripts
  use command-line arguments or environment variables and should be preferred.
- `02_figure_scripts/`: original community-level figure scripts, retained in
  their established manuscript order.
- `04_K08356_analysis_extensions/`: 100 K08356-specific scripts added from the
  analysis workspace. These include both current scripts and exact historical
  analysis snapshots.

Fifty extension scripts retain one or more historical `/Users/catherine` or
`/home/ps` paths. They are intentionally labelled as provenance snapshots in
the extension README. Before reuse, replace those paths with the current data
layout. Statistical methods and plotting choices were not mechanically changed
during repository organization.

Credential audit patterns covered login passwords, GitHub token prefixes,
explicit password/token assignments, and the historical HPC host address. No
credential-like values were detected in the committed script set.

Generated figures, intermediate tables, HMM models, protein FASTA files, raw
reads, R libraries, and cache directories remain excluded by `.gitignore`.

## Syntax validation

- All 22 shell scripts passed `bash -n`.
- All 35 Python scripts compiled successfully with Python 3.
- 89 of 94 R scripts parsed successfully.
- Five pre-existing legacy/archive R files did not parse: three contain pasted
  shell-session commands and two contain embedded NUL characters. They remain
  only for provenance and are not recommended execution targets.

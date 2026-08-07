# 808-MAG AioA/IdrA HMM screen

This workflow screens the complete protein sets of dereplicated MAGs before any
K08356-based candidate filtering. It therefore tests whether the KOfam-first
workflow missed AioA/IdrA-related proteins.

## Models and interpretation

Obtain the three HMM files from the first-author repository:

<https://github.com/vmreyes12/Neighborhood_Analyses/tree/main/hmm>

Required files:

- `combined_iriA_aioA.hmm`: primary discovery model. The author repository
  explicitly associates this model with a score threshold of 640.
- `iriA_new.hmm`: 17-sequence IdrA-specific model. A score threshold of 640 is
  retained here only as a sensitivity analysis because the repository does not
  document a separate calibrated cutoff for this model.
- `aioA.hmm`: 19-sequence AioA comparator model; its 640-score screen is also
  exploratory.

An HMM hit is a candidate, not a final functional assignment. Candidate naming
must be confirmed with an AioA/IdrA/DMSOR reference tree and gene-neighborhood
evidence.

## Usage

```bash
bash 01_run_808MAG_IdrA_HMM_screen.sh \
  /path/to/dereplicated_MAG_faa \
  /path/to/Neighborhood_Analyses/hmm \
  /path/to/output_directory \
  32 \
  808 \
  0
```

Arguments are FAA directory, model directory, output directory, CPU count,
expected MAG count, and expected protein count. Set the final argument to `0`
to disable the protein-count assertion. Input files must end in `.cds.faa`.

The script prefixes each protein identifier with its MAG identifier, checks for
empty files and duplicate target IDs, copies and hashes the models, runs three
independent `-T 640` searches, and writes tidy hit and summary tables.

Compare the full-MAG hits with an existing candidate FASTA and extract newly
detected proteins:

```bash
python 02_compare_extract_new_candidates.py \
  --hits /path/to/output/03_tables/all_T640_hits_long.tsv \
  --known-fasta /path/to/existing_candidates.faa \
  --combined-fasta /path/to/output/00_input/all_808MAG_proteins_prefixed.faa \
  --output-dir /path/to/output/03_tables/candidate_comparison
```

New hits are separated into primary combined-model hits and exploratory
IdrA-specific-only hits. Both groups still require phylogenetic and neighborhood
validation.

## Evidence boundary

The primary screen reproduces the public author-repository setting
`combined_iriA_aioA.hmm + T640`. The IdrA-specific and AioA-specific screens are
reported separately so their uncalibrated exploratory thresholds cannot be
mistaken for the publication's primary cutoff.

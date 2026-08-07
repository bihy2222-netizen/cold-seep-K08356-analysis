# Standardized group-derep DMSOR evidence workflow

Run on the analysis server:

```bash
cd /path/to/standardized_derep_DMSOR_pipeline
bash run_standardized_DMSOR_pipeline.sh
```

The workflow keeps these evidence layers separate:

- `01_hmm/`: author combined HMM at `-T 640`, plus descriptive IdrA/AioA scores.
- `02_neighborhood/`: target-centered ±10 ORFs and edge-censoring status.
- `03_reference_search/`: matches to literature-validated IdrB/AioB-related and P-like references.
- `04_tables/synteny_evidence.tsv`: independent synteny score and `synteny_status`.

Do not use `synteny_status` to define `phylogenetic_class`. Generate
`final_evidence_grade` only after both independent layers have been reviewed.

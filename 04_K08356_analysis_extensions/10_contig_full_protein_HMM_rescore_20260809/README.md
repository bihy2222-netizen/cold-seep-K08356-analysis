# Full-contig protein HMM discovery

This workflow scans all predicted proteins from the sample-level metagenomic
contigs without using K08356 as a discovery gate. It is separate from the
completed 808-derepMAG scan.

The source directory contains 57 FAA files, but
`SY457BB8-12_AA.faa` is an exact truncated prefix of
`SY457BB-8-12_AA.faa` and ends in a partial FASTA header. It is retained in the
audit manifest as `excluded_truncated_duplicate`; 56 complete sample FAA files
are included in the actual scan.

## Workflow

1. Audit every FAA, detect within-file duplicate IDs, add the sample prefix,
   and record input and merged-library SHA-256 values.
2. Search the complete prefixed protein library with the IdrA, AioA, and
   combined profiles at `-T 100 --domT 20` using normal HMMER acceleration.
3. Extract the three-model candidate union.
4. Re-score only the candidate union with each profile using
   `--max -T 0 --domT 0`.
5. Report sample, contig, ORF, full score, best-domain score, HMM coverage,
   protein coverage, and IdrA-minus-AioA delta score.

Broad T100 hits are DMSOR-related discovery candidates, not IdrA assignments.
Functional classification requires a reference phylogeny and gene-neighborhood
validation.

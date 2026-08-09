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

## NCBI pause point

`04_prepare_ncbiblast_and_tree_inputs.py` combines the completed MAG and contig
score tables only for audited candidate preparation. Priority NCBI files use
the conservative criterion that at least one exact IdrA, AioA, or combined
full-sequence score is at least 640, followed by 100% amino-acid sequence
dereplication. Every original source remains in the membership tables.

The script also copies the exact 136-reference FASTA used by the previous small
tree and records its SHA-256. It does not run alignment, phylogeny, functional
naming, or gene-neighborhood analysis. Those steps remain paused until the MAG
and contig NCBI Protein BLAST HitTable files are returned.

# Maggie DIAMOND blastx competitive pilot

This analysis tests whether nucleotide divergence caused Bowtie2 to miss distant
IdrA homologs. Clean metatranscriptomic cDNA reads are searched against the 63
protein competitive DMSOR reference library with DIAMOND blastx 2.0.14.

## Search parameters

```text
--very-sensitive --evalue 1e-5 --query-cover 60 --top 5 --threads 16
```

R1 and R2 are searched separately. Read IDs are normalized by removing a terminal
`/1`, `/2`, or equivalent whitespace mate suffix before read-pair aggregation.
`alignments`, unique read ends, unique pair IDs, and pairs whose two ends receive
the same best-family assignment are reported separately.

Each read end is assigned to its highest-bitscore family. If the best scores from
two different families differ by less than 10 bits, the read is classified as
`ambiguous_DMSOR`. Strict and partial IdrA are always reported separately.

Protein reference intervals and merged coverage breadth are reported. The
`conserved_region_only` field is currently `not_assessed_no_domain_coordinates`
because a validated Mo-bisPGD/Pfam coordinate manifest was not available in the
runtime environment. This field must not be interpreted as PASS or FAIL.

Short or localized hits alone are not considered expression evidence. Robust
strict IdrA support requires strict IdrA to remain the competitive best family,
clear separation from the second-best family, and broad multi-region coverage.

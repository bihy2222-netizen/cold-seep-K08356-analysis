# DMSOR competitive nucleotide references

These tables document external tree-reference DMSOR protein accessions whose
nucleotide CDS were recovered through NCBI `fasta_cds_na` and validated by exact
translation back to the tree protein sequence.

Files:

- `tree_reference_cds_manifest.tsv`: accepted NapA/ArrA/DmsA CDS references.
- `cds_fetch_audit.tsv`: success and failure records from accession-to-CDS fetch.
- `translation_validation.tsv`: translation status for accepted references.

Do not use protein HMMs or protein FASTA as Bowtie2 targets. Build the final
Maggie competitive mapping reference from nucleotide CDS only.

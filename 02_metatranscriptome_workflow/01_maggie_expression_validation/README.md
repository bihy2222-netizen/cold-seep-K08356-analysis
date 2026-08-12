# Maggie metatranscriptome expression validation

This folder contains the reference-building step for RNA read mapping to confirmed
IdrA-associated nucleotide CDS references.

## Method boundary

Protein HMMs such as `iriA_new.hmm`, `aioA.hmm`, or
`combined_iriA_aioA.hmm` are not nucleotide mapping targets. Use them to define
protein candidates first, recover the corresponding CDS (`.fna`), then map RNA
reads to the CDS with Bowtie2 or another nucleotide aligner.

Accepted expression validation requires:

1. high-confidence IdrA-associated CDS references;
2. exact CDS-to-protein translation checks;
3. a competitive DMSOR nucleotide reference containing non-target families such
   as AioA, NapA, ArrA, and DmsA;
4. RNA read mapping metrics such as mapped pairs, coverage breadth, mean depth,
   and TPM/RPKM-like normalized abundance.

HMM bitscores alone are not expression evidence.

## Build references

Example server paths used in the 2026-08-12 audit:

```bash
ROOT=/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra/maggie_expression_validation_20260812
T=/home/ps/ps1/data/bihongyu/cold_seep/HMM_supplement/808MAG_T640_union58_smalltree_refs136_20260810/04_tables/candidate58_four_clade_integrated_evidence.tsv
N=/home/ps/ps1/data/bihongyu/cold_seep/HMM_supplement/808MAG_T640_union58_smalltree_refs136_20260810/04_neighborhood/union58_neighborhood_long.tsv
PROD=/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/all_50_10_dRep/dereplicated_genomes/bin_prodigal
FAA=/home/ps/ps1/data/bihongyu/cold_seep/illu/binning/all_50_10_dRep/dereplicated_genomes/bin_faa

python3 scripts/build_maggie_refs.py \
  --evidence "$T" \
  --neighborhood "$N" \
  --prodigal-dir "$PROD" \
  --faa-dir "$FAA" \
  --outdir "$ROOT"
```

## Expected outputs

- `high_confidence_IdrA_CDS.fna`
- `high_confidence_IdrA_proteins.faa`
- `competitive_DMSOR_CDS.fna`
- `competitive_DMSOR_proteins.faa`
- `02_high_confidence_IdrA_locus_audit.tsv`
- `02_high_confidence_IdrA_translation_check.tsv`
- `02_reference_manifest.tsv`
- `02_reference_duplicate_audit.tsv`
- `BLOCKED_REFERENCE_LIBRARY.md`

If reliable nucleotide CDS for required non-target DMSOR families are missing,
the workflow should stop before accepted competitive Bowtie2 validation.

# Maggie metatranscriptome expression validation

This folder contains the reproducible reference-building, nucleotide mapping,
and protein-level competitive validation steps for metatranscriptomic evidence
of IdrA-associated genes.

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

RNA-seq FASTQ reads are cDNA nucleotide reads. They can be mapped directly to
validated, strand-normalized IdrA CDS sequences with Bowtie2; translation of
the reads is not required for this primary analysis. DIAMOND blastx is an
independent supplement for detecting distant homologs that nucleotide mapping
may miss.

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


## Add external tree-reference DMSOR CDS

The broad DMSOR tree may contain NapA, ArrA, and DmsA protein references. For
competitive RNA read mapping, recover their nucleotide CDS through NCBI
`fasta_cds_na` and keep only records that translate exactly back to the tree
protein reference.

```bash
BROAD=/home/ps/ps1/data/bihongyu/cold_seep/ASproteintree/aioA-like_tree_results_20260629_final/aioA-like-contigstree-5/broad_idra_dmsor_refs.large_subunit.filtered.faa
ARR=/home/ps/ps1/data/bihongyu/cold_seep/ASproteintree/aioA-like_tree_results_20260629_final/aioA-like-contigstree-5/ColdSeep_K08356_plus_Fig5_plus_broad_DMSOR_refs_refDedup_shortLabels_allNonstandardToX.faa

python3 scripts/fetch_tree_dmsor_cds.py   --broad-faa "$BROAD"   --arr-faa "$ARR"   --outdir "$ROOT"   --per-class 6
```

In the 2026-08-12 run, this recovered exact-translation nucleotide references
for NapA (6), ArrA (4), and DmsA (6). These were merged with the MAG-derived
DMSOR CDS to create `competitive_DMSOR_plus_treeRefs_CDS.fna` for the accepted
Maggie competitive mapping pilot.

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

## Maggie Bowtie2 expression validation

`run_maggie_mapping_pilot.py` performs two separate paired-end mappings:

1. primary mapping against 20 strict, synteny-supported IdrA-associated CDS;
2. competitive mapping against 63 DMSOR-family CDS, including strict and
   partial IdrA, AioA, unknown DMSOR, NapA, ArrA, and DmsA.

Example:

```bash
python3 scripts/run_maggie_mapping_pilot.py \
  --root "$MAGGIE_ROOT" \
  --sample SRR19238834_JL_0.1 \
  --r1 "$CLEAN_R1" \
  --r2 "$CLEAN_R2" \
  --primary-ref "$MAGGIE_ROOT/high_confidence_IdrA_CDS.fna" \
  --competitive-ref "$MAGGIE_ROOT/competitive_DMSOR_plus_treeRefs_CDS.fna" \
  --threads 16
```

The workflow writes sorted BAM/BAI files, flagstat, idxstats, depth tables,
per-CDS summaries, primary-versus-competitive comparisons, and IGV inputs.
Short localized mappings or generic K08356/DMSOR matches are not sufficient to
claim IdrA expression.

## DIAMOND blastx competitive validation

`run_maggie_diamond_blastx.sh` runs the original two QDN samples against the
63-protein competitive library. `summarize_maggie_diamond_blastx.py` reports
alignments, unique read ends, normalized pair IDs, best and second-best family,
bitscore difference, protein intervals, and merged reference coverage breadth.

R1 and R2 are searched separately. A read-end assignment is marked ambiguous
when the best scores from different DMSOR families differ by less than 10 bits.
Strict and partial IdrA are never merged. See
`DIAMOND_BLASTX_METHODS.md` for the exact parameters and current limitation of
the conserved-domain overlap field.

## JL_0.1 paired DNA-RNA validation

The Jiaolong cold seep sample `JL_0.1` represents 8-10 cm sediment and has a
paired DNA and RNA record under BioSample `SAMN27753712`:

| Molecule | Correct run | Library source | Role |
| --- | --- | --- | --- |
| DNA | `SRR19020591` | METAGENOMIC | Screen K08356 and establish IdrA genomic presence |
| RNA | `SRR19238834` | METATRANSCRIPTOMIC | Maggie expression validation |

Do not label `SRR19020591` as RNA. Both records use the non-standard strategy
label `OTHER`; their molecular types are supported by the NCBI experiment title
and the matching ENA/NCBI `LibrarySource` fields.

`download_and_run_SRR19238834_JL_0.1_RNA.sh` performs fresh sequential R1/R2
downloads with at most five attempts per mate. Each mate must pass official byte
count, MD5, and `gzip -t`; analysis starts only after both mates pass. It then
runs fastp, primary/competitive Bowtie2, and the sample-specific DIAMOND blastx
workflow. Existing invalid final FASTQ files are not overwritten.

The corrected crosswalk, RNA FASTQ manifest, and metadata audit are retained in
`metadata_examples/`. Raw NCBI/ENA responses remain in the analysis audit
directory on the server and are not committed to this repository.

## Repository contents

```text
01_maggie_expression_validation/
├── README.md
├── DIAMOND_BLASTX_METHODS.md
├── metadata_examples/
│   ├── JL_0.1_DNA_RNA_crosswalk.tsv
│   ├── JL_0.1_METADATA_AUDIT.md
│   └── JL_0.1_RNA_fastq_manifest.tsv
└── scripts/
    ├── build_maggie_refs.py
    ├── download_and_run_SRR19238834_JL_0.1_RNA.sh
    ├── fetch_tree_dmsor_cds.py
    ├── run_maggie_diamond_blastx.sh
    ├── run_maggie_mapping_pilot.py
    ├── run_SRR19238834_JL_0.1_diamond_blastx.sh
    └── summarize_maggie_diamond_blastx.py
```

# South China Sea metatranscriptome IdrA validation

This folder contains a resumable workflow for screening South China Sea cold-seep sediment metatranscriptomes for IdrA/DIRM-like transcriptional evidence.

Primary target projects:

- `PRJNA739005`: Qiongdongnan/QDN metatranscriptomes, first priority.
- `PRJNA831433`: Jiaolong and Haima metatranscriptomes, second priority.
- `PRJNA738468`: S11/Haima-related metatranscriptomes, supporting comparison.

The workflow intentionally keeps two evidence routes separate:

1. Direct read mapping to curated strict/extended/competitive IdrA CDS references.
2. De novo transcript assembly followed by Prodigal and `idra.hmm` search.

Important: `idra.hmm` is a protein profile-HMM. It must not be used as a Bowtie2 mapping target for RNA reads. Bowtie2 maps RNA reads to nucleotide CDS references such as `idra_strict_20_CDS.fna`. The HMM route is only:

```text
RNA reads -> transcript assembly -> Prodigal proteins -> hmmsearch idra.hmm
```

Do not interpret all K08356/DMSOR hits as IdrA. Final labels must be assigned after phylogeny, neighborhood/context checks, and comparison against AioA/other DMSOR references.

## Server Root

```bash
/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra
```

The scripts must not write into the earlier public metagenome validation folder:

```bash
/home/ps/ps1/data/bihongyu/cold_seep/public_validation/Dong_NewbornSeep_PRJNA1313825_20260803
```

## Suggested Order

```bash
cd /home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra

# 1. Check resources, software, HMM candidates, and reference candidates.
bash scripts/00_check_environment_and_inputs.sh

# 2. Fetch RunInfo and ENA FASTQ metadata, then build QDN/JL/Haima target tables.
python3 scripts/01_fetch_filter_scs_metatranscriptomes.py

# 3. Prepare/validate nucleotide CDS reference panels before Bowtie2.
python3 scripts/02_prepare_idra_cds_references.py \
  05_idra_reference/idra_reference_metadata.tsv

# 4. Review 00_metadata/PRJNA739005_QDN_RNA_fastq_manifest.tsv.

# 5. Download QDN RNA only after confirmation.
CONFIRM_QDN_RNA_DOWNLOAD=YES \
bash scripts/03_download_fastq_from_manifest.sh \
  00_metadata/PRJNA739005_QDN_RNA_fastq_manifest.tsv \
  02_raw_fastq/PRJNA739005_QDN_RNA

# 6. Run fastp, optional rRNA removal, Bowtie2 expression, assembly, Prodigal, and HMM.
bash scripts/04_run_qc_mapping_assembly_hmm.sh \
  00_metadata/PRJNA739005_QDN_RNA_fastq_manifest.tsv

# 7. Summarize expression and HMM evidence.
python3 scripts/05_summarize_idra_expression.py
python3 scripts/06_summarize_hmm_candidates.py
python3 scripts/07_make_final_evidence_matrix.py
```

## Core Outputs

```text
00_metadata/all_runs.tsv
00_metadata/all_RNA_runs.tsv
00_metadata/SCS_sediment_metatranscriptome_runs.tsv
00_metadata/QDN_priority_runs.tsv
00_metadata/download_manifest.tsv
00_metadata/excluded_runs.tsv
11_summary/fastp_summary.tsv
11_summary/rRNA_removal_summary.tsv
11_summary/idra_bowtie_expression.tsv
11_summary/idra_HMM_transcript_candidates.tsv
11_summary/idra_final_evidence_matrix.tsv
11_summary/SCS_metatranscriptome_sample_manifest.tsv
11_summary/README_results_CN.md
```

## Safety Notes

- Download scripts require explicit confirmation variables.
- Existing verified FASTQ files are reused.
- Existing MEGAHIT output directories are not overwritten.
- A missing or ambiguous `idra.hmm` should stop HMM interpretation but not metadata auditing.
- The historical combined `iriA/aioA` threshold `-T 640` is not automatically valid for a standalone `idra.hmm`.
- Bowtie2 panel mapping reports `IdrA_FPM`, `IdrA_RPKM`, read/fragment counts, mean depth, and breadth. It is not whole-transcriptome TPM unless reads are mapped to a complete gene catalog.

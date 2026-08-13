# 22_N30_16 MAG-context and per-base validation

This staged workflow implements the DNA side of Maggie's proposed evidence
chain for `22_N30_16 / SRR35213524`. It reuses the existing `pilot4_v2`
assembly, Prodigal, HMM, and binning outputs without overwriting them.

## Evidence boundary

`SRR35213524` is a DNA metagenome. Its read coverage is a genomic-presence and
technical positive control, not transcriptional evidence. Supplementary Table
1 reports metatranscriptomic sequencing for the physical sample, but no reliable
public RNA run has been linked. The RNA stage remains blocked until genuine RNA
FASTQ files are identified.

A combined-model score of at least 640 defines a high-scoring K08356/DMSOR
candidate only. A strict DIRM-like IdrA call additionally requires a strict-core
tree position, an IdrB-related neighbor, at least two P-like neighbors, compact
same-contig organization, and an observable non-edge-censored neighborhood.

## Output root

```text
pilot4_v2/22_N30_16/10_maggie_mag_context_validation/
├── 00_audit
├── 01_candidate_reconciliation
├── 02_neighborhood
├── 03_contig_bin_assignment
├── 04_personalized_references
├── 05_dna_presence_control_mapping
├── 06_per_base_depth
├── 07_figures
├── 08_igv
└── 09_final_tables
```

## Run order

1. Run `01_audit_22_N30_16_inputs.py` on the server. It inventories likely
   inputs, hashes them, and records ambiguous/missing file classes.
2. Supply the exact audited paths to `02_reconcile_22_N30_16_candidates.py`.
   This creates an exact-ID HMM matrix, candidate FASTA files, and +/-10 ORF
   neighborhoods. Tree fields remain `pending`.
3. Place candidates in the established broad DMSOR tree and annotate neighbors
   against the validated IdrB/AioB and P-like references. Record those results
   in the review table without changing candidate IDs. Use
   `plot_candidate_neighborhoods.py` to render one SVG/PDF/PNG arrow diagram per
   reviewed candidate.
4. Run `03_finalize_22_N30_16_candidates.py`. It enforces the strict evidence
   rules, maps candidate contigs exactly to bins, and builds strand-normalized
   personalized CDS/protein/contig references.
5. Run `04_run_22_N30_16_DNA_positive_control.sh` to competitively map clean DNA
   reads. Whole-assembly contig mapping is the primary localization evidence;
   all-bin mapping supplies host-MAG context; four individual CDS are an
   auxiliary count/sensitivity view; intact native clusters and candidate
   contigs support IGV and boundary inspection. It emits
   `samtools depth -aa` tables, per-reference summaries, BAM/BAI, and IGV files.

`run_22_N30_16_staged_workflow.sh` exposes these stages as `audit`,
`reconcile`, `annotate-neighborhood`, `finalize`, and `dna-map`. Each stage
requires explicit audited input paths through environment variables; no fuzzy
path selection is performed after the initial inventory.

The driver never downloads RNA and never labels DNA coverage as expression.

## Required manual review columns

The candidate review table must use these phylogenetic classes:

- `phylogenetically strict-core IdrA-associated`
- `phylogenetically partial IdrA-associated`
- `canonical AioA`
- `other DMSOR-family`
- `unresolved/truncated`

The finalizer derives the final integrated class; it does not accept a manually
typed final strict label without the required synteny fields.

For every accepted strict cluster, four separate CDS records (`IdrA`, `IdrB`,
`P-like-1`, and `P-like-2`) are retained for gene-level counts and per-base depth.
The native genomic interval and full contig are separate references for IGV and
boundary-spanning read/pair inspection. Neighbor coverage is supporting evidence
only and cannot rescue an IdrA CDS lacking reliable coverage.

The four-CDS mini-reference is a sensitivity/auxiliary view. Primary locus-level
coverage is extracted by genomic coordinates from the whole-assembly competitive
BAM, with a separate row for each of IdrA, IdrB, P-like-1, and P-like-2. This
allows unbinned and near-homologous assembly contigs to participate in the
competition and prevents a small target-only reference from inflating support.
The all-bin competitive BAM is retained separately for host-MAG context.

Both raw and MAPQ >= 20 depth tables are retained. Each strict-cluster CDS also
gets its own compressed per-base table. Per-reference alignment metrics separate
alignments, read names, proper-pair template names, MAPQ thresholds, and NM-tag
mismatch rate; these terms must not be treated as interchangeable counts.

Reads are aligned competitively against the complete assembly. To control disk
use, `samtools depth -aa` is then extracted only for reviewed candidate contigs
from that whole-assembly BAM; the alignment competition itself is not reduced.

Every DNA table, figure, and report is labeled `DNA presence control` or
`DNA-derived coverage; not transcriptional evidence`. Junction-spanning reads
and pairs are kept separate. A single read must physically span the complete
intergenic interval; a pair is counted only when its two mates overlap the two
adjacent genes. Paired support is auxiliary and is never interpreted as
co-transcription in this DNA stage.

# JL_0.1 paired DNA-RNA metadata audit

Audit date: 2026-08-13

JL_0.1 is a Jiaolong cold seep sediment sample from 8-10 cm below the seafloor.
ENA and NCBI independently associate both runs with BioSample SAMN27753712
(SRA sample SRS12837947) under PRJNA831433.

## Corrected run assignment

- SRR19020591 / SRX15092418: DNA metagenome. NCBI title is `metagenome of
  Jiaolong seep sediment`; LibrarySource is METAGENOMIC; layout is PAIRED;
  instrument is Illumina HiSeq 2500.
- SRR19238834 / SRX15299616: RNA metatranscriptome. NCBI title is
  `metatranscriptome of Jiaolong seep sediment`; LibrarySource is
  METATRANSCRIPTOMIC; layout is PAIRED; instrument is Illumina MiSeq.

Both records use the non-standard LibraryStrategy value OTHER. The molecular
assignment therefore relies on the concordant NCBI experiment title and
ENA/NCBI LibrarySource, not on rewriting OTHER as WGS or RNA-Seq.

No pre-existing SRR19020591, SRR19238834, or JL_0.1 FASTQ/analysis file was
found under the cold-seep project root before this audit.

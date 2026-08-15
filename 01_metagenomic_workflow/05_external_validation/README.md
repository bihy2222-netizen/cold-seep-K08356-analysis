# External metagenome validation

This workflow closes the evidence loop for the 82 public `T640` candidates
from `22_N10_3`, `22_N30_16`, and `23_N10_9`. Here, `T640` means a combined
K08356/DMSOR HMM score of at least 640. It is not an IdrA functional label.

`classify_external82_candidates.py` integrates:

- the existing public-candidate/ref IQ-TREE tree;
- the previously reviewed 48 core anchors from the four-clade evidence table;
- combined, `iriA_new`, and `aioA` HMM scores and domain coverage;
- exact Prodigal/GFF coordinates and contig-edge status;
- validated IdrB-related and P-like protein searches;
- exact contig-to-bin assignments, CheckM2 quality, CoverM abundance, and GTDB taxonomy.

The strict neighborhood rule is a compact, same-strand, non-overlapping
`candidate - IdrB-related - P-like-1 - P-like-2` arrangement at consecutive
ORF positions on either genomic side of the candidate. A truncated outer
neighborhood window is reported separately and does not invalidate four genes
that are themselves present and complete. P-like proteins are not renamed
IdrP1/IdrP2.

Outputs retain tree-anchor distances and the classes represented below the
candidate/nearest-anchor MRCA. Low-support or mixed-anchor placements remain
reviewable and are not silently promoted to strict IdrA.

```bash
bash run_external82_terminal_classification.sh
bash run_22_N30_16_dna_presence_control.sh
bash run_background255_gtdbtk_r226.sh
```

The second command reuses the existing indexed whole-assembly BAM. It reports
DNA-derived per-base coverage for each candidate and each member of complete
four-gene neighborhoods. These outputs are labelled as presence controls and
must not be interpreted as expression or co-transcription evidence.

The wrapper contains server paths for the current cold-seep project. Override
them through environment variables when reusing the workflow elsewhere.

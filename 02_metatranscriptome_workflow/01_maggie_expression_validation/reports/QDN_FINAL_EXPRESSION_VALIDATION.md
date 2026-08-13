# QDN competitive expression validation

## Final result

Neither QDN metatranscriptome provides strict IdrA-associated transcriptional
support after nucleotide- and protein-level competitive validation.

| Sample | Bowtie2 | DIAMOND blastx | Final classification |
| --- | --- | --- | --- |
| `SRR16610255 / QDN-W01B-1949` | No strict20 reads; localized DmsA-like mapping | Four isolated strict-reference hits without paired-end support; the competitive signal was dominated by DmsA/NapA | `reassigned to another DMSOR family` |
| `SRR16610253 / QDN-W04B-4900` | No strict20 reads; localized DmsA/NapA-like mapping | Two isolated strict-reference hits without paired-end support; the competitive signal was dominated by DmsA/NapA | `ambiguous DMSOR-family hit` |

The isolated strict-reference DIAMOND hits had bitscores of approximately
34-36, identities of approximately 34-43%, and covered only 3.4-5.3% of their
reference proteins. No read pair had both mates assigned to strict IdrA. These
short, localized hits therefore do not meet the criteria for tentative or
robust IdrA-associated expression.

## Manuscript-ready wording

> Protein-level translated searches identified only two to four isolated reads
> with weak similarity to strict IdrA references. These hits covered 3.4-5.3%
> of the reference proteins, lacked paired-end support, and were predominantly
> reassigned to DmsA- or NapA-like proteins in the competitive database. Thus,
> translated searches did not reveal remotely homologous IdrA transcripts
> missed by nucleotide mapping.

## Scope

This is a competitive negative result for these two QDN samples, not evidence
that IdrA-associated genes are absent from the wider QDN cold-seep system.


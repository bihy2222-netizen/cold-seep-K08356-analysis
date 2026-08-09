#!/usr/bin/env python3
"""Integrate untrimmed-MAFFT and MMseqs2 identities with global pairs."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--workdir", type=Path, required=True)
    return p.parse_args()


def norm_id(value: str) -> str:
    return value.split("|")[-1]


def main():
    args = parse_args()
    work = args.workdir.resolve()
    metadata = pd.read_csv(work / "00_input/sequence_group_metadata.tsv", sep="\t")
    group = dict(zip(metadata.sequence_id, metadata.final_function_group))
    order = metadata.sequence_id.tolist()

    aln_records = list(SeqIO.parse(work / "02_mafft/K08356_49_sequences.mafft_linsi.faa", "fasta"))
    aligned = {norm_id(r.id): str(r.seq).upper() for r in aln_records}
    if set(aligned) != set(order):
        raise SystemExit("MAFFT IDs do not match the curated 49-sequence metadata")
    aln_length = len(next(iter(aligned.values())))
    if any(len(s) != aln_length for s in aligned.values()):
        raise SystemExit("MAFFT alignment has inconsistent sequence lengths")
    complete_cols = [i for i in range(aln_length) if all(aligned[sid][i] != "-" for sid in order)]

    mafft_rows = []
    for a, b in itertools.combinations(order, 2):
        s1, s2 = aligned[a], aligned[b]
        non_double = [(x, y) for x, y in zip(s1, s2) if not (x == "-" and y == "-")]
        no_gap = [(x, y) for x, y in non_double if x != "-" and y != "-"]
        complete = [(s1[i], s2[i]) for i in complete_cols]
        id_non_double = sum(x == y for x, y in non_double)
        id_no_gap = sum(x == y for x, y in no_gap)
        id_complete = sum(x == y for x, y in complete)
        mafft_rows.append({
            "seq1": a, "seq2": b,
            "group1": group[a], "group2": group[b],
            "MAFFT_alignment_length": aln_length,
            "MAFFT_non_double_gap_columns": len(non_double),
            "MAFFT_no_gap_columns": len(no_gap),
            "MAFFT_complete_deletion_columns": len(complete),
            "MAFFT_identity_excluding_double_gaps": 100 * id_non_double / len(non_double),
            "MAFFT_identity_excluding_all_gap_columns": 100 * id_no_gap / len(no_gap),
            "MAFFT_pairwise_deletion_p_distance": 1 - id_no_gap / len(no_gap),
            "MAFFT_complete_deletion_identity": 100 * id_complete / len(complete),
            "MAFFT_complete_deletion_p_distance": 1 - id_complete / len(complete),
        })
    mafft = pd.DataFrame(mafft_rows)
    mafft.to_csv(work / "02_mafft/MAFFT_pairwise_identity_long.tsv", sep="\t", index=False)
    matrix = pd.DataFrame(np.eye(len(order)) * 100, index=order, columns=order)
    for row in mafft.itertuples():
        matrix.loc[row.seq1, row.seq2] = row.MAFFT_identity_excluding_double_gaps
        matrix.loc[row.seq2, row.seq1] = row.MAFFT_identity_excluding_double_gaps
    matrix.to_csv(work / "02_mafft/MAFFT_identity_matrix.tsv", sep="\t", index_label="sequence_id")

    cols = [
        "query", "target", "pident", "alnlen", "qstart", "qend", "qlen",
        "tstart", "tend", "tlen", "evalue", "bits"
    ]
    mm = pd.read_csv(work / "04_local_identity/mmseqs_all_vs_all.raw.tsv", sep="\t", names=cols)
    mm["query"] = mm["query"].map(norm_id)
    mm["target"] = mm["target"].map(norm_id)
    mm = mm[mm["query"] != mm["target"]].copy()
    mm["pair_key"] = mm.apply(lambda r: "\t".join(sorted((r["query"], r["target"]))), axis=1)
    mm = mm.sort_values(["pair_key", "bits"], ascending=[True, False]).drop_duplicates("pair_key")
    mm[["seq1", "seq2"]] = mm.pair_key.str.split("\t", expand=True)
    mm["group1"] = mm.seq1.map(group)
    mm["group2"] = mm.seq2.map(group)
    mm["qcov"] = 100 * mm.alnlen / mm.qlen
    mm["tcov"] = 100 * mm.alnlen / mm.tlen
    mm["min_cov"] = mm[["qcov", "tcov"]].min(axis=1)
    mm["max_cov"] = mm[["qcov", "tcov"]].max(axis=1)
    mm_out = mm[[
        "seq1", "seq2", "group1", "group2", "pident", "alnlen", "qlen", "tlen",
        "qstart", "qend", "tstart", "tend", "qcov", "tcov", "min_cov", "max_cov",
        "evalue", "bits"
    ]].rename(columns={"pident": "MMseqs2_local_identity", "bits": "MMseqs2_bitscore"})
    if len(mm_out) != 1176:
        raise SystemExit(f"Expected 1176 non-redundant MMseqs pairs, observed {len(mm_out)}")
    mm_out.to_csv(work / "04_local_identity/local_identity_all_pairs.tsv", sep="\t", index=False)
    mm_out[mm_out.min_cov >= 80].to_csv(work / "04_local_identity/local_identity_cov80.tsv", sep="\t", index=False)
    mm_out[mm_out.min_cov >= 90].to_csv(work / "04_local_identity/local_identity_cov90.tsv", sep="\t", index=False)
    local_matrix = pd.DataFrame(np.eye(len(order)) * 100, index=order, columns=order)
    for row in mm_out.itertuples():
        local_matrix.loc[row.seq1, row.seq2] = row.MMseqs2_local_identity
        local_matrix.loc[row.seq2, row.seq1] = row.MMseqs2_local_identity
    local_matrix.to_csv(work / "04_local_identity/MMseqs2_local_identity_matrix.tsv", sep="\t", index_label="sequence_id")

    global_df = pd.read_csv(work / "03_global_identity/global_pairwise_identity_long.tsv", sep="\t")
    global_df["pair_key"] = global_df.apply(lambda r: "\t".join(sorted((r.seq1, r.seq2))), axis=1)
    mafft["pair_key"] = mafft.apply(lambda r: "\t".join(sorted((r.seq1, r.seq2))), axis=1)
    mm_out["pair_key"] = mm_out.apply(lambda r: "\t".join(sorted((r.seq1, r.seq2))), axis=1)
    mm_merge = mm_out.drop(columns=["seq1", "seq2", "group1", "group2"]).rename(columns={
        "qcov": "MMseqs2_qcov", "tcov": "MMseqs2_tcov",
        "min_cov": "MMseqs2_min_cov", "max_cov": "MMseqs2_max_cov",
    })
    merged = (
        global_df.merge(mafft.drop(columns=["seq1", "seq2", "group1", "group2"]), on="pair_key")
        .merge(mm_merge, on="pair_key")
    )
    merged.to_csv(work / "06_statistics/global_MAFFT_MMseqs_pairwise_comparison.tsv", sep="\t", index=False)
    metrics = [
        "global_identity_with_gaps", "global_identity_no_gaps",
        "MAFFT_identity_excluding_double_gaps", "MAFFT_identity_excluding_all_gap_columns",
        "MMseqs2_local_identity", "MMseqs2_min_cov"
    ]
    merged[metrics].corr(method="pearson").to_csv(
        work / "06_statistics/identity_metric_Pearson_correlation.tsv", sep="\t"
    )
    merged[metrics].corr(method="spearman").to_csv(
        work / "06_statistics/identity_metric_Spearman_correlation.tsv", sep="\t"
    )
    print(f"MAFFT alignment: {len(order)} sequences, {aln_length} columns")
    print(f"Complete-deletion columns: {len(complete_cols)}")
    print(f"MMseqs2 non-redundant pairs: {len(mm_out)}")
    print(f"MMseqs2 min_cov >=80%: {(mm_out.min_cov >= 80).sum()}")
    print(f"MMseqs2 min_cov >=90%: {(mm_out.min_cov >= 90).sum()}")


if __name__ == "__main__":
    main()

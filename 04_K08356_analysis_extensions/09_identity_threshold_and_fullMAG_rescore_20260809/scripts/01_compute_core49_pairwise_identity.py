#!/usr/bin/env python3
"""Compute reproducible pairwise identities for the curated 49-protein set."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices


GROUP_MAP = {
    "Canonical AioA": "canonical_AioA",
    "DIRM-synteny IdrA": "strict_DIRM_like_IdrA",
    "IdrA phylogenetic": "partial_IdrA_associated",
    "Uncertain DMSOR": "Aio_like_or_unresolved",
}
GROUP_ORDER = list(GROUP_MAP.values())
GROUP_COLORS = {
    "canonical_AioA": "#D55E00",
    "strict_DIRM_like_IdrA": "#0072B2",
    "partial_IdrA_associated": "#009E73",
    "Aio_like_or_unresolved": "#CC79A7",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--workdir", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20260809)
    return p.parse_args()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def alignment_counts(aln, seq1: str, seq2: str, matrix) -> dict[str, float]:
    coords = np.asarray(aln.coordinates)
    identical = similar = nongap = gaps = aln_len = 0
    q_aligned = t_aligned = 0
    for k in range(coords.shape[1] - 1):
        a0, a1 = int(coords[0, k]), int(coords[0, k + 1])
        b0, b1 = int(coords[1, k]), int(coords[1, k + 1])
        da, db = a1 - a0, b1 - b0
        q_aligned += da
        t_aligned += db
        if da and db:
            if da != db:
                raise RuntimeError("Unexpected unequal aligned block")
            aa, bb = seq1[a0:a1], seq2[b0:b1]
            identical += sum(x == y for x, y in zip(aa, bb))
            similar += sum(matrix[x, y] > 0 for x, y in zip(aa, bb))
            nongap += da
            aln_len += da
        else:
            gap = max(da, db)
            gaps += gap
            aln_len += gap
    return {
        "alignment_length": aln_len,
        "identical_count": identical,
        "similar_count": similar,
        "nongap_aligned_columns": nongap,
        "gap_count": gaps,
        "q_aligned_residues": q_aligned,
        "t_aligned_residues": t_aligned,
    }


def make_aligner(mode: str):
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = Align.PairwiseAligner(mode=mode)
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner, matrix


def describe(values: np.ndarray, bootstrap: int, rng) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    medians = np.median(
        rng.choice(values, size=(bootstrap, len(values)), replace=True), axis=1
    )
    return {
        "n_pairs": len(values),
        "minimum": np.min(values),
        "P05": np.quantile(values, 0.05),
        "Q1": np.quantile(values, 0.25),
        "median": np.median(values),
        "mean": np.mean(values),
        "Q3": np.quantile(values, 0.75),
        "P95": np.quantile(values, 0.95),
        "maximum": np.max(values),
        "SD": np.std(values, ddof=1) if len(values) > 1 else 0.0,
        "IQR": np.quantile(values, 0.75) - np.quantile(values, 0.25),
        "median_bootstrap_CI_low": np.quantile(medians, 0.025),
        "median_bootstrap_CI_high": np.quantile(medians, 0.975),
    }


def save_matrix(long_df: pd.DataFrame, ids: list[str], value: str, path: Path):
    matrix = pd.DataFrame(np.eye(len(ids)) * 100.0, index=ids, columns=ids)
    for row in long_df.itertuples():
        matrix.loc[row.seq1, row.seq2] = getattr(row, value)
        matrix.loc[row.seq2, row.seq1] = getattr(row, value)
    matrix.to_csv(path, sep="\t", index_label="sequence_id")
    return matrix


def save_fig(fig, stem: Path):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=350, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    work = args.workdir.resolve()
    for d in (
        "01_qc", "02_mafft", "03_global_identity", "04_local_identity",
        "05_hmm_scores", "06_statistics", "07_figures", "08_threshold",
        "09_report", "logs"
    ):
        (work / d).mkdir(parents=True, exist_ok=True)

    faa49 = work / "00_input/K08356_49_sequences.faa"
    faa48 = work / "00_input/K08356_48_core_sequences.faa"
    meta_path = work / "00_input/four_clade_classification_49.tsv"
    seq_records = list(SeqIO.parse(faa49, "fasta"))
    normalized_ids = {r.id: r.id.split("|")[-1] for r in seq_records}
    seqs = {
        normalized_ids[r.id]: str(r.seq).rstrip("*").upper()
        for r in seq_records
    }
    raw_seqs = {
        normalized_ids[r.id]: str(r.seq).upper()
        for r in seq_records
    }
    original_headers = {normalized_ids[r.id]: r.id for r in seq_records}
    if len(seqs) != 49 or len(seqs) != len(seq_records):
        raise SystemExit(f"Expected 49 unique proteins, observed {len(seqs)}")

    meta = pd.read_csv(meta_path, sep="\t", dtype=str).fillna("")
    meta["final_function_group"] = meta["clade_display"].map(GROUP_MAP)
    if meta["final_function_group"].isna().any():
        raise SystemExit("Unmapped clade labels in metadata")
    meta_ids = set(meta.Sequence_ID)
    if set(seqs) != meta_ids:
        raise SystemExit(
            f"FAA/metadata mismatch: FAA-only={sorted(set(seqs)-meta_ids)}, "
            f"metadata-only={sorted(meta_ids-set(seqs))}"
        )
    meta = meta.set_index("Sequence_ID", drop=False)

    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
    qc_rows = []
    hash_to_ids = defaultdict(list)
    for sid in meta.Sequence_ID:
        seq = seqs[sid]
        raw = raw_seqs[sid]
        digest = hashlib.sha256(seq.encode()).hexdigest()
        hash_to_ids[digest].append(sid)
        invalid = sorted(set(seq) - allowed)
        qc_rows.append({
            "sequence_id": sid,
            "original_fasta_id": original_headers[sid],
            "MAG_id": meta.loc[sid, "MAG_ID"],
            "protein_length_aa": len(seq),
            "terminal_stop_removed": "yes" if raw.endswith("*") else "no",
            "internal_stop_count": raw[:-1].count("*"),
            "X_count": seq.count("X"),
            "invalid_characters": "".join(invalid),
            "sha256_clean_sequence": digest,
            "final_function_group": meta.loc[sid, "final_function_group"],
        })
    write_tsv(work / "01_qc/sequence_QC.tsv", qc_rows, list(qc_rows[0]))
    dup_rows = [
        {"sha256_clean_sequence": h, "n_sequences": len(ids), "sequence_ids": ";".join(ids)}
        for h, ids in hash_to_ids.items() if len(ids) > 1
    ]
    write_tsv(
        work / "01_qc/duplicate_sequences.tsv", dup_rows,
        ["sha256_clean_sequence", "n_sequences", "sequence_ids"]
    )

    old48 = {r.id for r in SeqIO.parse(faa48, "fasta")}
    audit = []
    for sid in sorted(set(seqs) | old48):
        audit.append({
            "sequence_id": sid,
            "present_in_old48": "yes" if sid in old48 else "no",
            "present_in_latest49": "yes" if sid in seqs else "no",
            "status": "shared" if sid in old48 and sid in seqs else (
                "added_in_latest49" if sid in seqs else "absent_from_latest49"
            ),
        })
    write_tsv(work / "01_qc/old48_vs_latest49_sequence_audit.tsv", audit, list(audit[0]))

    host_counts = Counter(meta.MAG_ID)
    metadata_rows = []
    for sid in meta.Sequence_ID:
        row = meta.loc[sid]
        metadata_rows.append({
            "sequence_id": sid,
            "MAG_id": row.MAG_ID,
            "sample": row.Sample_group,
            "habitat": row.Habitat,
            "host_taxonomy": "",
            "protein_length": len(seqs[sid]),
            "phylogenetic_clade": row.Tree_clade,
            "final_function_group": row.final_function_group,
            "neighbor_category": row.neighborhood_class,
            "strict_training_label": (
                "AioA_negative_anchor" if row.final_function_group == "canonical_AioA"
                else "IdrA_internal_sensitivity" if row.final_function_group == "strict_DIRM_like_IdrA"
                else "blind_test_only"
            ),
            "notes": "MAG contains two candidates" if host_counts[row.MAG_ID] > 1 else "",
        })
    write_tsv(work / "00_input/sequence_group_metadata.tsv", metadata_rows, list(metadata_rows[0]))

    order = []
    for group in GROUP_ORDER:
        order.extend(meta.loc[meta.final_function_group == group, "Sequence_ID"].tolist())
    global_aligner, matrix = make_aligner("global")
    local_aligner, _ = make_aligner("local")
    pair_rows = []
    for sid1, sid2 in itertools.combinations(order, 2):
        s1, s2 = seqs[sid1], seqs[sid2]
        g = global_aligner.align(s1, s2)[0]
        gc = alignment_counts(g, s1, s2, matrix)
        l = local_aligner.align(s1, s2)[0]
        lc = alignment_counts(l, s1, s2, matrix)
        group1 = meta.loc[sid1, "final_function_group"]
        group2 = meta.loc[sid2, "final_function_group"]
        gi, gj = GROUP_ORDER.index(group1), GROUP_ORDER.index(group2)
        category = " vs ".join((group1, group2) if gi <= gj else (group2, group1))
        pair_rows.append({
            "seq1": sid1, "seq2": sid2, "group1": group1, "group2": group2,
            "length1": len(s1), "length2": len(s2), **gc,
            "global_identity_with_gaps": 100 * gc["identical_count"] / gc["alignment_length"],
            "global_identity_no_gaps": 100 * gc["identical_count"] / gc["nongap_aligned_columns"],
            "global_similarity": 100 * gc["similar_count"] / gc["alignment_length"],
            "gap_fraction": 100 * gc["gap_count"] / gc["alignment_length"],
            "local_alignment_length": lc["alignment_length"],
            "local_identical_count": lc["identical_count"],
            "local_identity": 100 * lc["identical_count"] / lc["nongap_aligned_columns"],
            "qcov": 100 * lc["q_aligned_residues"] / len(s1),
            "tcov": 100 * lc["t_aligned_residues"] / len(s2),
            "min_cov": min(100 * lc["q_aligned_residues"] / len(s1), 100 * lc["t_aligned_residues"] / len(s2)),
            "global_alignment_score": float(g.score),
            "local_alignment_score": float(l.score),
            "comparison_category": category,
        })
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(work / "03_global_identity/global_pairwise_identity_long.tsv", sep="\t", index=False)
    local_cols = [
        "seq1", "seq2", "group1", "group2", "length1", "length2",
        "local_alignment_length", "local_identical_count", "local_identity",
        "qcov", "tcov", "min_cov", "local_alignment_score", "comparison_category"
    ]
    pairs[local_cols].to_csv(work / "04_local_identity/local_identity_all_pairs.tsv", sep="\t", index=False)
    pairs.loc[pairs.min_cov >= 80, local_cols].to_csv(
        work / "04_local_identity/local_identity_cov80.tsv", sep="\t", index=False
    )
    pairs.loc[pairs.min_cov >= 90, local_cols].to_csv(
        work / "04_local_identity/local_identity_cov90.tsv", sep="\t", index=False
    )

    global_matrix = save_matrix(
        pairs, order, "global_identity_with_gaps",
        work / "03_global_identity/global_identity_matrix.tsv"
    )
    save_matrix(pairs, order, "global_similarity", work / "03_global_identity/global_similarity_matrix.tsv")
    save_matrix(pairs, order, "gap_fraction", work / "03_global_identity/global_gap_matrix.tsv")
    local_matrix = save_matrix(
        pairs, order, "local_identity", work / "04_local_identity/local_identity_matrix.tsv"
    )

    rng = np.random.default_rng(args.seed)
    metric_names = [
        "global_identity_with_gaps", "global_identity_no_gaps", "global_similarity",
        "gap_fraction", "local_identity", "qcov", "tcov", "min_cov",
        "global_alignment_score", "local_alignment_score"
    ]
    summaries = []
    for category, sub in pairs.groupby("comparison_category", sort=False):
        for metric in metric_names:
            summaries.append({
                "comparison_category": category, "metric": metric,
                **describe(sub[metric].to_numpy(), args.bootstrap, rng),
            })
    summary = pd.DataFrame(summaries)
    summary.to_csv(work / "06_statistics/clade_pairwise_identity_summary.tsv", sep="\t", index=False)
    pairs.to_csv(work / "06_statistics/clade_pairwise_identity_all_values.tsv", sep="\t", index=False)

    priority_categories = [
        "strict_DIRM_like_IdrA vs strict_DIRM_like_IdrA",
        "canonical_AioA vs canonical_AioA",
        "canonical_AioA vs strict_DIRM_like_IdrA",
    ]
    priority = summary.loc[
        (summary.metric == "global_identity_with_gaps")
        & summary.comparison_category.isin(priority_categories)
    ].copy()
    priority.to_csv(work / "09_report/priority_global_identity_ranges.tsv", sep="\t", index=False)

    sns.set_theme(style="white", font="Times New Roman", font_scale=0.75)
    group_breaks = np.cumsum([sum(meta.final_function_group == g) for g in GROUP_ORDER])[:-1]
    for matrix_data, cmap, label, stem in (
        (global_matrix, "viridis", "Global identity (%)", "core49_global_identity_phylogenetic_clade_order"),
        (local_matrix, "mako", "Local identity (%)", "core49_local_identity_phylogenetic_clade_order"),
    ):
        fig, ax = plt.subplots(figsize=(11, 10))
        sns.heatmap(
            matrix_data.loc[order, order], cmap=cmap, vmin=0, vmax=100,
            cbar_kws={"label": label}, square=True, ax=ax,
            xticklabels=False, yticklabels=False
        )
        for boundary in group_breaks:
            ax.axhline(boundary, color="white", linewidth=1.2)
            ax.axvline(boundary, color="white", linewidth=1.2)
        ax.set_xlabel("Protein sequences ordered by curated clade")
        ax.set_ylabel("Protein sequences ordered by curated clade")
        save_fig(fig, work / f"07_figures/{stem}")

    plot_df = pairs.copy()
    category_order = []
    for i, g1 in enumerate(GROUP_ORDER):
        for g2 in GROUP_ORDER[i:]:
            category_order.append(f"{g1} vs {g2}")
    fig, ax = plt.subplots(figsize=(12, 7.5))
    sns.boxplot(
        data=plot_df, x="comparison_category", y="global_identity_with_gaps",
        order=category_order, color="#D9D9D9", fliersize=0, ax=ax
    )
    sns.stripplot(
        data=plot_df, x="comparison_category", y="global_identity_with_gaps",
        order=category_order, color="#222222", alpha=0.45, size=2.4, jitter=0.25, ax=ax
    )
    ax.set_xlabel("")
    ax.set_ylabel("Global amino-acid identity (%)")
    ax.tick_params(axis="x", rotation=42)
    ax.set_title("Within- and between-clade global identity of 49 curated K08356 proteins")
    save_fig(fig, work / "07_figures/core49_clade_pairwise_global_identity")

    counts = Counter(meta.final_function_group)
    with (work / "09_report/README_similarity_threshold_CN.md").open("w") as out:
        out.write("# K08356 四 clade 蛋白序列相似度分析（核心 49 条）\n\n")
        out.write("本分析单位为蛋白序列：49 条蛋白来自 48 个唯一 MAG。")
        out.write("HMM 分数沿用既有结果，不作为本轮 pairwise identity 的替代指标。\n\n")
        out.write("## 固定分组\n\n")
        for group in GROUP_ORDER:
            out.write(f"- `{group}`: {counts[group]}\n")
        out.write("\n## 全局比对口径\n\n")
        out.write("Biopython PairwiseAligner 全局模式，BLOSUM62，gap-open 10，")
        out.write("gap-extension 0.5。主指标为 identical residues / full global alignment length（含单边 gap）。\n")

    report_map = {
        "K08356_4clade_similarity_full.tsv": work / "06_statistics/clade_pairwise_identity_all_values.tsv",
        "K08356_pairwise_global_identity_matrix.tsv": work / "03_global_identity/global_identity_matrix.tsv",
        "K08356_pairwise_local_identity_matrix.tsv": work / "04_local_identity/local_identity_matrix.tsv",
        "K08356_4clade_similarity_summary.tsv": work / "06_statistics/clade_pairwise_identity_summary.tsv",
    }
    for name, source in report_map.items():
        (work / "09_report" / name).write_bytes(source.read_bytes())

    print(f"Completed {len(pairs)} non-redundant pairs")
    print(priority[["comparison_category", "minimum", "median", "maximum"]].to_string(index=False))


if __name__ == "__main__":
    main()

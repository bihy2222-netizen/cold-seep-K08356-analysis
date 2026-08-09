#!/usr/bin/env python3
"""Summarize exact three-model HMM rescoring of the T640 union."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


MODELS = ("IdrA", "AioA", "combined")
GROUP_NAMES = {
    "Canonical AioA": "canonical AioA-associated",
    "DIRM-synteny IdrA": "synteny-supported strict DIRM-like IdrA",
    "IdrA phylogenetic": "partial IdrA-associated",
    "Uncertain DMSOR": "AioA-like or unresolved DMSOR",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--workdir", type=Path, required=True)
    p.add_argument("--core49", type=Path, required=True)
    p.add_argument("--new10", type=Path, required=True)
    return p.parse_args()


def parse_tbl(path: Path):
    rows = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split(maxsplit=18)
            rows[f[0]] = {
                "full_evalue": float(f[4]),
                "bitscore": float(f[5]),
                "full_bias": float(f[6]),
            }
    return rows


def union_length(intervals):
    total = 0
    end = 0
    for start, stop in sorted(intervals):
        if stop <= end:
            continue
        total += stop - max(start, end) + 1
        end = stop
    return total


def parse_domtbl(path: Path):
    domains = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split(maxsplit=22)
            domains[f[0]].append({
                "protein_length": int(f[2]),
                "hmm_length": int(f[5]),
                "domain_c_evalue": float(f[11]),
                "domain_i_evalue": float(f[12]),
                "domain_bitscore": float(f[13]),
                "hmm_from": int(f[15]),
                "hmm_to": int(f[16]),
                "ali_from": int(f[17]),
                "ali_to": int(f[18]),
                "env_from": int(f[19]),
                "env_to": int(f[20]),
                "accuracy": float(f[21]),
            })
    out = {}
    for sid, ds in domains.items():
        hmm_intervals = [(d["hmm_from"], d["hmm_to"]) for d in ds]
        protein_intervals = [(d["ali_from"], d["ali_to"]) for d in ds]
        best = max(ds, key=lambda d: d["domain_bitscore"])
        out[sid] = {
            "domain_count": len(ds),
            "best_domain_bitscore": best["domain_bitscore"],
            "best_domain_i_evalue": best["domain_i_evalue"],
            "hmm_length": best["hmm_length"],
            "protein_length": best["protein_length"],
            "hmm_coverage_pct": 100 * union_length(hmm_intervals) / best["hmm_length"],
            "protein_coverage_pct": 100 * union_length(protein_intervals) / best["protein_length"],
            "best_hmm_from": best["hmm_from"],
            "best_hmm_to": best["hmm_to"],
            "best_ali_from": best["ali_from"],
            "best_ali_to": best["ali_to"],
        }
    return out


def save_fig(fig, stem: Path):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=350, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    work = args.workdir.resolve()
    (work / "05_tables").mkdir(exist_ok=True)
    (work / "06_figures").mkdir(exist_ok=True)
    with args.core49.open() as handle:
        core = {r["Sequence_ID"]: r for r in csv.DictReader(handle, delimiter="\t")}
    with args.new10.open() as handle:
        new = {r["protein_ID"].rstrip("\r"): r for r in csv.DictReader(handle, delimiter="\t")}

    all_tbl = {}
    all_dom = {}
    for model in MODELS:
        all_tbl[model] = parse_tbl(work / f"02_rescore58/58_union.{model}.tbl")
        all_dom[model] = parse_domtbl(work / f"02_rescore58/58_union.{model}.domtbl")
    ids = sorted(set(all_tbl["IdrA"]))
    if len(ids) != 58 or any(set(all_tbl[m]) != set(ids) for m in MODELS):
        raise SystemExit("Expected an identical 58-protein set for all three models")

    rows = []
    for target in ids:
        mag_id, protein_id = target.split("|", 1)
        if protein_id in core:
            source_class = "core49_T640_intersection"
            group = GROUP_NAMES[core[protein_id]["clade_display"]]
        elif protein_id in new:
            source_class = "T640_new_blind_candidate"
            group = "blind_unclassified"
        else:
            source_class = "legacy_core48_not_in_latest49_metadata"
            group = "metadata_check_required"
        row = {
            "prefixed_target_ID": target,
            "protein_ID": protein_id,
            "MAG_ID": mag_id,
            "source_class": source_class,
            "curated_or_blind_group": group,
        }
        for model in MODELS:
            t = all_tbl[model][target]
            d = all_dom[model][target]
            for key, value in t.items():
                row[f"{model}_{key}"] = value
            for key, value in d.items():
                row[f"{model}_{key}"] = value
        length = row["IdrA_protein_length"]
        row["IdrA_minus_AioA_bitscore"] = row["IdrA_bitscore"] - row["AioA_bitscore"]
        row["normalized_IdrA_bitscore"] = row["IdrA_bitscore"] / length
        row["normalized_AioA_bitscore"] = row["AioA_bitscore"] / length
        row["normalized_delta_bitscore"] = row["IdrA_minus_AioA_bitscore"] / length
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(work / "05_tables/union58_three_model_score_coverage_matrix.tsv", sep="\t", index=False)
    df[[
        "protein_ID", "MAG_ID", "source_class", "curated_or_blind_group",
        "IdrA_bitscore", "AioA_bitscore", "combined_bitscore",
        "IdrA_minus_AioA_bitscore", "normalized_delta_bitscore"
    ]].to_csv(work / "05_tables/union58_three_model_score_matrix.tsv", sep="\t", index=False)
    coverage_cols = ["protein_ID", "MAG_ID", "source_class"] + [
        f"{m}_{x}" for m in MODELS
        for x in ("hmm_coverage_pct", "protein_coverage_pct", "domain_count", "best_domain_bitscore")
    ]
    df[coverage_cols].to_csv(work / "05_tables/union58_HMM_coverage_matrix.tsv", sep="\t", index=False)

    counts = (
        df.groupby("source_class")
        .agg(protein_count=("protein_ID", "nunique"), MAG_count=("MAG_ID", "nunique"))
        .reset_index()
    )
    counts.to_csv(work / "05_tables/union58_source_counts.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid", font="Times New Roman")
    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    sns.scatterplot(
        data=df, x="AioA_bitscore", y="IdrA_bitscore", hue="curated_or_blind_group",
        style="source_class", s=70, alpha=0.85, ax=ax
    )
    lim = max(df.AioA_bitscore.max(), df.IdrA_bitscore.max()) * 1.03
    ax.plot([0, lim], [0, lim], linestyle="--", color="#555555", linewidth=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("AioA HMM bitscore")
    ax.set_ylabel("IdrA HMM bitscore")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    save_fig(fig, work / "06_figures/union58_IdrA_vs_AioA_HMM_scores")

    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    order = [
        "canonical AioA-associated", "synteny-supported strict DIRM-like IdrA",
        "partial IdrA-associated", "AioA-like or unresolved DMSOR",
        "blind_unclassified"
    ]
    sns.boxplot(
        data=df, x="curated_or_blind_group", y="IdrA_minus_AioA_bitscore",
        order=order, color="#D9D9D9", fliersize=0, ax=ax
    )
    sns.stripplot(
        data=df, x="curated_or_blind_group", y="IdrA_minus_AioA_bitscore",
        order=order, color="#222222", size=3, alpha=0.65, jitter=0.2, ax=ax
    )
    ax.axhline(0, color="#D55E00", linestyle="--", linewidth=1)
    ax.tick_params(axis="x", rotation=30)
    ax.set_xlabel("")
    ax.set_ylabel("IdrA minus AioA HMM bitscore")
    save_fig(fig, work / "06_figures/union58_delta_bitscore_distribution")

    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Summarize the full-808-MAG T100 discovery union and exact rescoring."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO


MODELS = ("IdrA", "AioA", "combined")
BINS = [
    (640, float("inf"), ">=640"),
    (500, 640, "500-<640"),
    (400, 500, "400-<500"),
    (300, 400, "300-<400"),
    (200, 300, "200-<300"),
    (100, 200, "100-<200"),
    (-float("inf"), 100, "<100"),
]


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
            rows[f[0]] = {"evalue": float(f[4]), "bitscore": float(f[5]), "bias": float(f[6])}
    return rows


def union_len(intervals):
    total = 0
    end = 0
    for start, stop in sorted(intervals):
        if stop > end:
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
                "tlen": int(f[2]), "qlen": int(f[5]),
                "domain_score": float(f[13]), "domain_i_evalue": float(f[12]),
                "hmm_from": int(f[15]), "hmm_to": int(f[16]),
                "ali_from": int(f[17]), "ali_to": int(f[18]),
            })
    out = {}
    for sid, ds in domains.items():
        best = max(ds, key=lambda x: x["domain_score"])
        out[sid] = {
            "protein_length": best["tlen"], "hmm_length": best["qlen"],
            "domain_count": len(ds), "best_domain_bitscore": best["domain_score"],
            "best_domain_i_evalue": best["domain_i_evalue"],
            "hmm_coverage_pct": 100 * union_len([(d["hmm_from"], d["hmm_to"]) for d in ds]) / best["qlen"],
            "protein_coverage_pct": 100 * union_len([(d["ali_from"], d["ali_to"]) for d in ds]) / best["tlen"],
        }
    return out


def score_bin(value):
    for low, high, label in BINS:
        if low <= value < high:
            return label
    raise RuntimeError(value)


def main():
    args = parse_args()
    work = args.workdir.resolve()
    outdir = work / "05_tables"
    outdir.mkdir(exist_ok=True)
    with args.core49.open() as handle:
        core = {r["Sequence_ID"]: r for r in csv.DictReader(handle, delimiter="\t")}
    with args.new10.open() as handle:
        new10 = {r["protein_ID"].rstrip("\r") for r in csv.DictReader(handle, delimiter="\t")}
    k08356 = {r.id.split("~~")[-1] for r in SeqIO.parse(work / "00_input/K08356_METABOLIC_47.faa", "fasta")}

    discovery = {}
    exact = {}
    domains = {}
    for model in MODELS:
        discovery[model] = parse_tbl(work / f"03_full808_T100/full808.{model}.T100.tbl")
        exact[model] = parse_tbl(work / f"04_T100_union_rescore/T100_union.{model}.tbl")
        domains[model] = parse_domtbl(work / f"04_T100_union_rescore/T100_union.{model}.domtbl")
    ids = sorted(set().union(*(set(discovery[m]) for m in MODELS)))
    if len(ids) != 720 or any(set(exact[m]) != set(ids) for m in MODELS):
        raise SystemExit("T100 union or exact-rescore set is inconsistent")

    rows = []
    for target in ids:
        mag, protein = target.split("|", 1)
        if protein in core:
            source = "core49_T640_intersection"
            group = core[protein]["clade_display"]
        elif protein in new10:
            source = "T640_new_blind_candidate"
            group = "blind_unclassified"
        else:
            source = "T100_discovery_new_candidate"
            group = "blind_unclassified"
        row = {
            "prefixed_target_ID": target, "protein_ID": protein, "MAG_ID": mag,
            "source_class": source, "curated_or_blind_group": group,
            "K08356_METABOLIC_annotation": "positive" if protein in k08356 else "negative",
        }
        for model in MODELS:
            t, d = exact[model][target], domains[model][target]
            row[f"{model}_T100_discovery_hit"] = "yes" if target in discovery[model] else "no"
            row[f"{model}_bitscore"] = t["bitscore"]
            row[f"{model}_evalue"] = t["evalue"]
            row[f"{model}_score_bin"] = score_bin(t["bitscore"])
            for key, value in d.items():
                row[f"{model}_{key}"] = value
        length = row["IdrA_protein_length"]
        row["IdrA_minus_AioA_bitscore"] = row["IdrA_bitscore"] - row["AioA_bitscore"]
        row["normalized_delta_bitscore"] = row["IdrA_minus_AioA_bitscore"] / length
        row["possible_truncation_lt700aa"] = "yes" if length < 700 else "no"
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "T100_union720_three_model_exact_score_coverage.tsv", sep="\t", index=False)

    bins = []
    for model in MODELS:
        for label in [x[2] for x in BINS]:
            sub = df[df[f"{model}_score_bin"] == label]
            bins.append({
                "model": model, "score_bin": label,
                "protein_count": sub.protein_ID.nunique(), "MAG_count": sub.MAG_ID.nunique(),
                "K08356_positive": (sub.K08356_METABOLIC_annotation == "positive").sum(),
                "K08356_negative": (sub.K08356_METABOLIC_annotation == "negative").sum(),
            })
    pd.DataFrame(bins).to_csv(outdir / "T100_union_score_bin_summary.tsv", sep="\t", index=False)

    source_summary = (
        df.groupby(["source_class", "K08356_METABOLIC_annotation"])
        .agg(protein_count=("protein_ID", "nunique"), MAG_count=("MAG_ID", "nunique"))
        .reset_index()
    )
    source_summary.to_csv(outdir / "T100_union_source_K08356_overlap.tsv", sep="\t", index=False)
    pd.crosstab(df.source_class, df.K08356_METABOLIC_annotation).to_csv(
        outdir / "T100_union_source_by_K08356_contingency.tsv", sep="\t"
    )

    k_only = sorted(k08356 - set(df.protein_ID))
    pd.DataFrame({"K08356_protein_not_in_T100_union": k_only}).to_csv(
        outdir / "K08356_positive_not_in_T100_union.tsv", sep="\t", index=False
    )
    print(f"T100 union: {len(df)} proteins from {df.MAG_ID.nunique()} MAGs")
    print(source_summary.to_string(index=False))
    print(f"K08356 positives absent from T100 union: {len(k_only)}")


if __name__ == "__main__":
    main()

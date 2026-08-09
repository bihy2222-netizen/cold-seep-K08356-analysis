#!/usr/bin/env python3
"""Summarize three-model exact rescoring for the contig T100 candidate union."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


MODELS = ("IdrA", "AioA", "combined")


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", type=Path, required=True)
    return parser.parse_args()


def parse_tbl(path: Path):
    rows = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(maxsplit=18)
            rows[fields[0]] = {
                "full_evalue": float(fields[4]),
                "bitscore": float(fields[5]),
                "bias": float(fields[6]),
            }
    return rows


def parse_domtbl(path: Path):
    domains = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(maxsplit=22)
            domains[fields[0]].append({
                "protein_length": int(fields[2]),
                "hmm_length": int(fields[5]),
                "domain_i_evalue": float(fields[12]),
                "domain_score": float(fields[13]),
                "hmm_from": int(fields[15]),
                "hmm_to": int(fields[16]),
                "ali_from": int(fields[17]),
                "ali_to": int(fields[18]),
            })
    return domains


def union_length(intervals):
    total = 0
    end = 0
    for start, stop in sorted(intervals):
        if stop > end:
            total += stop - max(start, end) + 1
            end = stop
    return total


def domain_summary(domains):
    if not domains:
        return {}
    best = max(domains, key=lambda row: row["domain_score"])
    return {
        "protein_length": best["protein_length"],
        "hmm_length": best["hmm_length"],
        "domain_count": len(domains),
        "best_domain_score": best["domain_score"],
        "best_domain_i_evalue": best["domain_i_evalue"],
        "hmm_coverage_pct": 100 * union_length(
            [(row["hmm_from"], row["hmm_to"]) for row in domains]
        ) / best["hmm_length"],
        "protein_coverage_pct": 100 * union_length(
            [(row["ali_from"], row["ali_to"]) for row in domains]
        ) / best["protein_length"],
    }


def main() -> None:
    work = arguments().workdir.resolve()
    discovery = {
        model: parse_tbl(work / f"03_T100_discovery/contig.{model}.T100.tbl")
        for model in MODELS
    }
    exact = {
        model: parse_tbl(work / f"05_union_exact_rescore/contig_union.{model}.exact.tbl")
        for model in MODELS
    }
    domains = {
        model: parse_domtbl(work / f"05_union_exact_rescore/contig_union.{model}.exact.domtbl")
        for model in MODELS
    }
    ids = sorted(set().union(*(set(rows) for rows in discovery.values())))
    output_rows = []
    for target in ids:
        sample, protein = target.split("|", 1)
        contig = protein.rsplit("_", 1)[0] if "_" in protein else protein
        row = {
            "target_id": target,
            "sample_id": sample,
            "contig_id": contig,
            "orf_id": protein,
        }
        for model in MODELS:
            row[f"{model}_T100_discovery_hit"] = "yes" if target in discovery[model] else "no"
            score = exact[model].get(target)
            row[f"{model}_full_score"] = score["bitscore"] if score else "NA"
            row[f"{model}_full_evalue"] = score["full_evalue"] if score else "NA"
            summary = domain_summary(domains[model].get(target, []))
            for field in (
                "protein_length", "hmm_length", "domain_count", "best_domain_score",
                "best_domain_i_evalue", "hmm_coverage_pct", "protein_coverage_pct",
            ):
                row[f"{model}_{field}"] = summary.get(field, "NA")
        idra = row["IdrA_full_score"]
        aioa = row["AioA_full_score"]
        row["IdrA_minus_AioA_delta_score"] = (
            idra - aioa if idra != "NA" and aioa != "NA" else "NA"
        )
        output_rows.append(row)

    outdir = work / "06_tables"
    outdir.mkdir(exist_ok=True)
    with (outdir / "contig_T100_union_three_model_score_coverage.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    counts = []
    for model in MODELS:
        model_ids = set(discovery[model])
        counts.append({
            "model": model,
            "T100_protein_count": len(model_ids),
            "sample_count": len({target.split("|", 1)[0] for target in model_ids}),
        })
    counts.append({
        "model": "three_model_union",
        "T100_protein_count": len(ids),
        "sample_count": len({target.split("|", 1)[0] for target in ids}),
    })
    with (outdir / "contig_T100_model_hit_counts.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(counts[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(counts)
    print(counts)


if __name__ == "__main__":
    main()

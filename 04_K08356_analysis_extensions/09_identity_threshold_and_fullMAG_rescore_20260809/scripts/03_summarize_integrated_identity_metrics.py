#!/usr/bin/env python3
"""Generate integrated four-clade summaries and sequence-level outlier tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


METRICS = [
    "global_identity_with_gaps",
    "global_identity_no_gaps",
    "MAFFT_identity_excluding_double_gaps",
    "MAFFT_identity_excluding_all_gap_columns",
    "MAFFT_complete_deletion_identity",
    "MMseqs2_local_identity",
    "MMseqs2_qcov",
    "MMseqs2_tcov",
    "MMseqs2_min_cov",
    "MMseqs2_bitscore",
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--workdir", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20260809)
    return p.parse_args()


def describe(values, bootstrap, rng):
    x = np.asarray(values, dtype=float)
    med = np.median(rng.choice(x, size=(bootstrap, len(x)), replace=True), axis=1)
    return {
        "n_pairs": len(x), "minimum": x.min(), "P05": np.quantile(x, 0.05),
        "Q1": np.quantile(x, 0.25), "median": np.median(x), "mean": x.mean(),
        "Q3": np.quantile(x, 0.75), "P95": np.quantile(x, 0.95), "maximum": x.max(),
        "SD": x.std(ddof=1) if len(x) > 1 else 0.0,
        "IQR": np.quantile(x, 0.75) - np.quantile(x, 0.25),
        "median_bootstrap_CI_low": np.quantile(med, 0.025),
        "median_bootstrap_CI_high": np.quantile(med, 0.975),
    }


def save_fig(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=350, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    work = args.workdir.resolve()
    df = pd.read_csv(work / "06_statistics/global_MAFFT_MMseqs_pairwise_comparison.tsv", sep="\t")
    rng = np.random.default_rng(args.seed)
    rows = []
    for category, sub in df.groupby("comparison_category", sort=False):
        for metric in METRICS:
            rows.append({"comparison_category": category, "metric": metric, **describe(sub[metric], args.bootstrap, rng)})
    summary = pd.DataFrame(rows)
    summary.to_csv(work / "06_statistics/clade_pairwise_identity_summary_integrated.tsv", sep="\t", index=False)
    summary.to_csv(work / "09_report/K08356_4clade_similarity_summary_integrated.tsv", sep="\t", index=False)

    long = []
    for row in df.itertuples():
        long.append({"sequence_id": row.seq1, "group": row.group1, "partner_id": row.seq2, "partner_group": row.group2, "identity": row.global_identity_with_gaps})
        long.append({"sequence_id": row.seq2, "group": row.group2, "partner_id": row.seq1, "partner_group": row.group1, "identity": row.global_identity_with_gaps})
    long = pd.DataFrame(long)
    within = long[long.group == long.partner_group]
    outliers = (
        within.groupby(["sequence_id", "group"])
        .identity.agg(["count", "min", "median", "mean", "max"])
        .reset_index()
        .rename(columns={"count": "within_group_pair_count", "min": "within_group_min", "max": "within_group_max"})
    )
    outliers["within_group_median_percentile"] = outliers.groupby("group")["median"].rank(pct=True) * 100
    outliers["low_within_group_similarity_flag"] = np.where(outliers.within_group_median_percentile <= 10, "yes", "no")
    outliers.to_csv(work / "06_statistics/sequence_within_group_identity_outlier_audit.tsv", sep="\t", index=False)

    strict_within = df[df.comparison_category == "synteny-supported strict DIRM-like IdrA vs synteny-supported strict DIRM-like IdrA"].global_identity_with_gaps
    aio_strict = df[df.comparison_category == "canonical AioA-associated vs synteny-supported strict DIRM-like IdrA"].global_identity_with_gaps
    gap = strict_within.min() - aio_strict.max()
    decision = pd.DataFrame([{
        "dataset": "local_curated_core49_only",
        "max_canonical_AioA_associated_vs_synteny_supported_strict_DIRM_like_IdrA": aio_strict.max(),
        "min_synteny_supported_strict_DIRM_like_IdrA_within_clade": strict_within.min(),
        "identity_nonoverlap_gap_percentage_points": gap,
        "local_empirical_nonoverlap": "yes" if gap > 0 else "no",
        "universal_threshold_recommended": "no_pending_external_validated_references",
    }])
    decision.to_csv(work / "08_threshold/local_core49_identity_nonoverlap_audit.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid", font="Times New Roman")
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    sns.scatterplot(
        data=df, x="global_identity_with_gaps", y="MAFFT_identity_excluding_double_gaps",
        hue="comparison_category", s=25, alpha=0.65, ax=ax, legend=False
    )
    low = min(df.global_identity_with_gaps.min(), df.MAFFT_identity_excluding_double_gaps.min())
    ax.plot([low, 100], [low, 100], linestyle="--", color="#555555", linewidth=1)
    ax.set_xlabel("Global pairwise identity, including gaps (%)")
    ax.set_ylabel("Untrimmed MAFFT identity, excluding double gaps (%)")
    save_fig(fig, work / "07_figures/global_vs_MAFFT_identity_consistency")

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    sns.scatterplot(
        data=df, x="global_identity_with_gaps", y="MMseqs2_local_identity",
        size="MMseqs2_min_cov", hue="comparison_category", sizes=(15, 90),
        alpha=0.55, ax=ax, legend=False
    )
    ax.set_xlabel("Global pairwise identity, including gaps (%)")
    ax.set_ylabel("MMseqs2 local identity (%)")
    save_fig(fig, work / "07_figures/global_vs_MMseqs2_local_identity")

    print(decision.to_string(index=False))


if __name__ == "__main__":
    main()

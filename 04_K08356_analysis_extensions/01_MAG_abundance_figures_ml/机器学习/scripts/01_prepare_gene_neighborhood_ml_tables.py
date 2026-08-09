#!/usr/bin/env python3
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

BASE = Path("/Users/catherine/Downloads/MAG-bin-tpm")
ML_DIR = BASE / "机器学习"
INPUT_DIR = ML_DIR / "00_input"
MODEL_DIR = ML_DIR / "01_gene_neighborhood_model"
OUT_DIR = ML_DIR / "outputs"

COORDS = BASE / "K08356_MAG_gene_island_arrows_20260708" / "K08356_four_branch_gene_neighborhood_coordinates.csv"
REMOTE_ORF_PATH = "/home/ps/ps1/data/bihongyu/cold_seep/illu/all_fa/CDS/cat_result/cat_others"

FUNCTION_CLASSES = [
    "Rieske/AioB",
    "Cytochrome/peroxidase",
    "Arsenic resistance",
    "Transporter",
    "Nitrogen metabolism",
    "Sox/sulfur oxidation",
    "Redox/oxidoreductase",
    "Annotated other",
    "Hypothetical/unknown",
]

SAFE = {
    "Rieske/AioB": "Rieske_AioB",
    "Cytochrome/peroxidase": "cyt_peroxidase",
    "Arsenic resistance": "arsenic_resistance",
    "Transporter": "transporter",
    "Nitrogen metabolism": "nitrogen_metabolism",
    "Sox/sulfur oxidation": "sox_sulfur_oxidation",
    "Redox/oxidoreductase": "redox_oxidoreductase",
    "Annotated other": "annotated_other",
    "Hypothetical/unknown": "hypothetical_unknown",
}


def contig_id_from_target(target_gene_id: str) -> str:
    short = target_gene_id.split("-")[-1]
    if "k141_" in target_gene_id:
        return "k141_" + target_gene_id.split("k141_", 1)[1].rsplit("_", 1)[0]
    return ""


def mag_id_from_target(target_gene_id: str) -> str:
    if "-k141_" in target_gene_id:
        return target_gene_id.split("-k141_", 1)[0]
    return target_gene_id.rsplit("_", 2)[0]


def target_short_from_target(target_gene_id: str) -> str:
    if "k141_" in target_gene_id:
        return "k141_" + target_gene_id.split("k141_", 1)[1]
    return target_gene_id


def main():
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    coords = pd.read_csv(COORDS)
    coords["relative_position"] = coords["Local_order"] - coords["Target_local_order"]
    coords["strand"] = coords["Direction"].map({1: "+", -1: "-"}).fillna(coords["Direction"].astype(str))
    coords["source_remote_orf_path"] = REMOTE_ORF_PATH

    target_rows = coords[coords["Is_target"] == "yes"].copy()
    label_cols = [
        "Original_MAG_ID",
        "MAG",
        "Target_short",
        "Target_contig",
        "Feature",
        "Corrected_habitat",
        "Corrected_sample_group",
        "Sample",
        "KO",
        "Description",
    ]
    label = target_rows[label_cols].drop_duplicates().rename(
        columns={
            "Original_MAG_ID": "target_gene_id",
            "MAG": "MAG_id",
            "Target_short": "target_orf_id",
            "Target_contig": "contig_id",
            "Feature": "branch",
            "Corrected_habitat": "habitat",
            "Corrected_sample_group": "sample_group",
            "Sample": "sample_id",
            "KO": "target_KO",
            "Description": "target_product",
        }
    )
    label["source_derep_strategy"] = "overall_dRep_MAG"
    label["source_remote_orf_path"] = REMOTE_ORF_PATH

    orf = coords.rename(
        columns={
            "Original_MAG_ID": "target_gene_id",
            "Protein_short": "neighbor_orf_id",
            "Feature": "branch",
            "Corrected_habitat": "habitat",
            "Corrected_sample_group": "sample_group",
            "Sample": "sample_id",
            "Start": "start",
            "End": "end",
            "Description": "product",
            "Gene_class": "function_class",
            "Gene_label": "function_label",
        }
    )[
        [
            "target_gene_id",
            "MAG",
            "branch",
            "habitat",
            "sample_group",
            "sample_id",
            "Target_contig",
            "neighbor_orf_id",
            "Contig",
            "relative_position",
            "strand",
            "start",
            "end",
            "KO",
            "product",
            "function_class",
            "function_label",
            "Is_target",
            "source_remote_orf_path",
        ]
    ].rename(
        columns={
            "MAG": "MAG_id",
            "Target_contig": "target_contig",
            "Contig": "contig_id",
            "Is_target": "is_target",
        }
    )
    orf["is_target"] = orf["is_target"].map({"yes": True, "no": False})

    feature_rows = []
    neighbor_only = orf[~orf["is_target"]].copy()
    for target_id, grp in neighbor_only.groupby("target_gene_id", sort=False):
        meta = label[label["target_gene_id"] == target_id].iloc[0].to_dict()
        row = {
            "target_gene_id": target_id,
            "MAG_id": meta["MAG_id"],
            "contig_id": meta["contig_id"],
            "branch": meta["branch"],
            "habitat": meta["habitat"],
            "sample_group": meta["sample_group"],
            "sample_id": meta["sample_id"],
            "n_neighbor_orfs": int(len(grp)),
            "n_upstream_orfs": int((grp["relative_position"] < 0).sum()),
            "n_downstream_orfs": int((grp["relative_position"] > 0).sum()),
        }
        for cls in FUNCTION_CLASSES:
            safe = SAFE[cls]
            cls_grp = grp[grp["function_class"] == cls]
            rel = cls_grp["relative_position"].astype(int).tolist()
            row[f"n_{safe}"] = int(len(cls_grp))
            row[f"has_{safe}"] = int(len(cls_grp) > 0)
            row[f"n_up_{safe}"] = int((cls_grp["relative_position"] < 0).sum())
            row[f"n_down_{safe}"] = int((cls_grp["relative_position"] > 0).sum())
            row[f"min_abs_orf_distance_{safe}"] = min([abs(x) for x in rel], default=math.nan)
            row[f"has_within3_{safe}"] = int(any(abs(x) <= 3 for x in rel))
        feature_rows.append(row)

    feature = pd.DataFrame(feature_rows)
    branch_order = [
        "IdrA-associated",
        "canonical aioA-associated",
        "aioA-like-associated",
        "unknown AioA-like / uncertain DMSOR",
    ]
    feature["branch"] = pd.Categorical(feature["branch"], categories=branch_order, ordered=True)
    feature = feature.sort_values(["branch", "habitat", "target_gene_id"]).reset_index(drop=True)

    prevalence_rows = []
    for cls in FUNCTION_CLASSES:
        safe = SAFE[cls]
        for branch, grp in feature.groupby("branch", observed=False):
            prevalence_rows.append(
                {
                    "branch": branch,
                    "function_class": cls,
                    "n_targets": int(len(grp)),
                    "n_present": int(grp[f"has_{safe}"].sum()),
                    "prevalence": float(grp[f"has_{safe}"].mean()) if len(grp) else math.nan,
                    "mean_count": float(grp[f"n_{safe}"].mean()) if len(grp) else math.nan,
                }
            )
    prevalence = pd.DataFrame(prevalence_rows)

    label_path = INPUT_DIR / "K08356_branch_label_table.csv"
    orf_path = INPUT_DIR / "K08356_gene_neighborhood_ORF_table.csv"
    feature_path = MODEL_DIR / "gene_neighborhood_feature_matrix.csv"
    prevalence_path = MODEL_DIR / "gene_neighborhood_feature_prevalence_by_branch.csv"

    label.to_csv(label_path, index=False)
    orf.to_csv(orf_path, index=False)
    feature.to_csv(feature_path, index=False)
    prevalence.to_csv(prevalence_path, index=False)

    readme = ML_DIR / "README_ML_validation_v1.md"
    readme.write_text(
        "\n".join(
            [
                "# K08356 ML-assisted validation v1",
                "",
                "This is a small local test focused on MAG-level K08356 gene-neighborhood features.",
                "",
                "## Source",
                f"- ORF coordinates and annotations were derived from `{COORDS}`.",
                f"- The underlying sample-level ORF/GFF source path provided by the user is `{REMOTE_ORF_PATH}`.",
                "- No new files were written to the supercomputer source directory.",
                "",
                "## Feature design",
                "- Each row in the feature matrix is one target K08356 neighborhood.",
                "- Predictors use only neighboring ORFs; the target K08356 ORF itself is excluded from model features to avoid label leakage.",
                "- Both count and binary features were generated for functional classes used in the gene-island figure.",
                "",
                "## First-pass modeling",
                "- PCA and heatmap are descriptive.",
                "- PERMANOVA tests whether gene-neighborhood feature composition differs among branches.",
                "- Leave-one-out nearest-centroid classification is used as a lightweight validation suitable for small n.",
                "- Decision-tree feature importance is exploratory only.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Wrote {len(label)} label rows")
    print(f"Wrote {len(orf)} ORF rows")
    print(f"Wrote {len(feature)} feature rows")
    print(label_path)
    print(orf_path)
    print(feature_path)
    print(prevalence_path)


if __name__ == "__main__":
    main()

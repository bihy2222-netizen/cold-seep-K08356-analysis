#!/usr/bin/env python3
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests


BASE = Path("/Users/catherine/Downloads/MAG-bin-tpm")
OUT = BASE / "Rhodobacteraceae_Clade4_vs_K08356_negative_20260805"
TAXONOMY = Path("/Users/catherine/Downloads/MAG 宿主是样品里有哪些 MAG，这些 MAG 的相对丰度么/merged_tpm_matrix 处理哦吼.xlsx")
KO_ANNOTATION = Path("/Users/catherine/Downloads/MAG KEGG/all_bin_kofamscan_filtered_split.20260121200414199.xlsx")
MODULE_DEFINITION = Path("/Users/catherine/Downloads/MAG 代谢图/metabolism_summary-师兄给的原始数据.xlsx")
CLADE_MAPPING = Path("/Users/catherine/Downloads/MAG 代谢图/K08356_final_clade_mapping_48.csv")

ANALYSIS_SHEETS = [
    "MISC", "carbon utilization", "carbon utilization (Woodcroft)",
    "Transporters", "N", "Organic Nitrogen",
]

CATEGORY_RULES = {
    "Carbon": r"Methane|Methanogenesis|Acetyl-CoA pathway, CO2|Wood-Ljungdahl|Citrate cycle|Glycolysis|Pyruvate oxidation|TCA /Reductive TCA",
    "Nitrogen": r"Nitrogen fixation|Nitrification|Denitrification|nitrate reduction|Nitrate assimilation|nitrite \+ ammonia|Urea cycle|Urea transport",
    "Sulfur": r"sulfate reduction|Thiosulfate oxidation|Sulfate/thiosulfate transport|Sulfonate transport",
    "Arsenic": r"^Arsenate$|Arsenic|arsenite|arsenate",
    "Vitamin B12": r"Cobalamin|Vitamin B12|cobinamide",
    "Motility": r"^Flagellar Assembly$|Chemotaxis",
}


def habitat_from_mag(mag):
    if mag.startswith("SY"):
        return "IS"
    if re.match(r"^(S1_|S2_|S4_|SQ_)", mag):
        return "AS"
    if re.match(r"^(C1_|C2_|C3_|S13_|S14_|S15_|ES_)", mag):
        return "ES"
    if re.match(r"^(S3_|NS_|R2111_)", mag):
        return "NS"
    return "Unknown"


def clean_taxon(value, prefix):
    value = "" if pd.isna(value) else str(value)
    value = re.sub(rf"^{prefix}__", "", value)
    return value if value else "Unclassified"


def read_module_catalog():
    pieces = []
    for sheet in ANALYSIS_SHEETS:
        frame = pd.read_excel(MODULE_DEFINITION, sheet_name=sheet, usecols=lambda column: column in {"gene_id", "module"})
        frame = frame.dropna(subset=["gene_id"]).copy()
        frame["gene_id"] = frame["gene_id"].astype(str).str.strip()
        frame["module"] = frame["module"].fillna("Unclassified").astype(str).str.strip().replace("", "Unclassified")
        frame["source_sheet"] = sheet
        pieces.append(frame[["source_sheet", "module", "gene_id"]].drop_duplicates())
    return pd.concat(pieces, ignore_index=True).drop_duplicates()


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    adjusted = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    if valid.any():
        adjusted[valid] = multipletests(values[valid], method="fdr_bh")[1]
    return adjusted


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    taxonomy = pd.read_excel(TAXONOMY, sheet_name="Sheet3")
    taxonomy = taxonomy.loc[taxonomy["family"].astype(str).str.contains("Rhodobacteraceae", case=False, na=False)].copy()
    taxonomy = taxonomy.rename(columns={"user_genome": "MAG"})
    taxonomy["genus_clean"] = taxonomy["genus"].map(lambda value: clean_taxon(value, "g"))
    taxonomy["species_clean"] = taxonomy["species"].map(lambda value: clean_taxon(value, "s"))
    taxonomy["habitat"] = taxonomy["MAG"].map(habitat_from_mag)

    clades = pd.read_csv(CLADE_MAPPING)
    clade4 = set(clades.loc[clades["final_clade"].eq("Clade 4"), "MAG"])
    all_mapped_k08356 = set(clades["MAG"])

    raw_ko = pd.read_excel(KO_ANNOTATION, header=None, names=["MAG", "gene", "KO", "description"])
    raw_ko = raw_ko.dropna(subset=["MAG", "KO"])
    raw_ko["MAG"] = raw_ko["MAG"].astype(str).str.strip()
    raw_ko["KO"] = raw_ko["KO"].astype(str).str.strip()
    rhodo_mags = set(taxonomy["MAG"])
    raw_rhodo = raw_ko.loc[raw_ko["MAG"].isin(rhodo_mags)].copy()
    raw_k08356_mags = set(raw_rhodo.loc[raw_rhodo["KO"].eq("K08356"), "MAG"])

    taxonomy["group"] = np.select(
        [taxonomy["MAG"].isin(clade4), ~taxonomy["MAG"].isin(all_mapped_k08356 | raw_k08356_mags)],
        ["Clade 4 carrier", "K08356-negative"],
        default="Other K08356 clade",
    )
    comparison = taxonomy.loc[taxonomy["group"].isin(["Clade 4 carrier", "K08356-negative"])].copy()

    module_catalog = read_module_catalog()
    ko_counts = raw_rhodo.groupby(["MAG", "KO"]).size().rename("gene_count").reset_index()
    module_scores = module_catalog.merge(ko_counts, left_on="gene_id", right_on="KO", how="left")
    module_scores = module_scores.loc[module_scores["MAG"].isin(comparison["MAG"])].copy()
    module_sizes = module_catalog.groupby(["source_sheet", "module"])["gene_id"].nunique().rename("gene_set_size")
    detected = module_scores.loc[module_scores["gene_count"].fillna(0).gt(0)].groupby(
        ["source_sheet", "module", "MAG"]
    ).agg(genes_detected=("gene_id", "nunique"), total_gene_count=("gene_count", "sum")).reset_index()

    grid = pd.MultiIndex.from_product(
        [comparison["MAG"], module_sizes.index], names=["MAG", "module_key"]
    ).to_frame(index=False)
    grid[["source_sheet", "module"]] = pd.DataFrame(grid.pop("module_key").tolist(), index=grid.index)
    scores = grid.merge(detected, on=["MAG", "source_sheet", "module"], how="left")
    scores = scores.merge(module_sizes.reset_index(), on=["source_sheet", "module"], how="left")
    scores[["genes_detected", "total_gene_count"]] = scores[["genes_detected", "total_gene_count"]].fillna(0)
    scores["coverage_pct"] = 100 * scores["genes_detected"] / scores["gene_set_size"]
    scores["present"] = scores["genes_detected"].gt(0)
    scores = scores.merge(comparison[["MAG", "group", "habitat", "genus_clean", "species_clean"]], on="MAG", how="left")

    rows = []
    for (source_sheet, module), subset in scores.groupby(["source_sheet", "module"], sort=False):
        carrier = subset.loc[subset["group"].eq("Clade 4 carrier")]
        negative = subset.loc[subset["group"].eq("K08356-negative")]
        table = [[carrier["present"].sum(), (~carrier["present"]).sum()], [negative["present"].sum(), (~negative["present"]).sum()]]
        odds_ratio, fisher_p = fisher_exact(table)
        mw = mannwhitneyu(carrier["coverage_pct"], negative["coverage_pct"], alternative="two-sided")
        rows.append({
            "source_sheet": source_sheet, "module": module,
            "n_clade4": len(carrier), "n_negative": len(negative),
            "clade4_prevalence_pct": 100 * carrier["present"].mean(),
            "negative_prevalence_pct": 100 * negative["present"].mean(),
            "prevalence_difference_pp": 100 * (carrier["present"].mean() - negative["present"].mean()),
            "fisher_odds_ratio": odds_ratio, "fisher_p": fisher_p,
            "clade4_mean_coverage_pct": carrier["coverage_pct"].mean(),
            "negative_mean_coverage_pct": negative["coverage_pct"].mean(),
            "mean_coverage_difference_pp": carrier["coverage_pct"].mean() - negative["coverage_pct"].mean(),
            "mann_whitney_p": mw.pvalue,
        })
    stats = pd.DataFrame(rows)
    stats["fisher_fdr"] = bh_adjust(stats["fisher_p"])
    stats["mann_whitney_fdr"] = bh_adjust(stats["mann_whitney_p"])

    category_rows = []
    for category, pattern in CATEGORY_RULES.items():
        selected = stats["module"].str.contains(pattern, case=False, regex=True, na=False)
        category_rows.append(stats.loc[selected].assign(category=category))
    focused = pd.concat(category_rows, ignore_index=True).drop_duplicates(["source_sheet", "module"])
    focused["rank_score"] = focused["mean_coverage_difference_pp"].abs() + focused["prevalence_difference_pp"].abs() / 2
    plot_stats = focused.sort_values("rank_score", ascending=False).groupby("category", group_keys=False).head(4)
    plot_stats = plot_stats.sort_values(["category", "mean_coverage_difference_pp"], ascending=[True, False])

    plot_keys = list(plot_stats[["source_sheet", "module"]].itertuples(index=False, name=None))
    plot_scores = scores.set_index(["source_sheet", "module"]).loc[plot_keys].reset_index()
    matrix = plot_scores.pivot(index="MAG", columns="module", values="coverage_pct")
    meta = comparison.set_index("MAG").loc[matrix.index]
    order = meta.assign(group_order=meta["group"].map({"Clade 4 carrier": 0, "K08356-negative": 1})).sort_values(
        ["group_order", "genus_clean", "habitat"]
    ).index
    matrix = matrix.loc[order, plot_stats["module"]]
    meta = meta.loc[order]

    plot_scores.to_csv(OUT / "Rhodobacteraceae_heatmap_long_data.csv", index=False)
    meta.reset_index().to_csv(OUT / "Rhodobacteraceae_heatmap_MAG_order.csv", index=False)

    taxonomy.sort_values(["group", "habitat", "genus_clean", "MAG"]).to_csv(OUT / "Rhodobacteraceae_MAG_group_assignments.csv", index=False)
    scores.to_csv(OUT / "Rhodobacteraceae_all_module_scores_by_MAG.csv", index=False)
    stats.sort_values(["mann_whitney_fdr", "fisher_fdr", "module"]).to_csv(OUT / "Rhodobacteraceae_module_difference_statistics.csv", index=False)
    plot_stats.to_csv(OUT / "Rhodobacteraceae_heatmap_selected_modules.csv", index=False)
    pd.DataFrame({
        "metric": ["Rhodobacteraceae MAGs", "Clade 4 carriers", "K08356-negative MAGs", "Other K08356-clade MAGs", "Modules tested", "Focused modules plotted"],
        "value": [len(taxonomy), (taxonomy["group"] == "Clade 4 carrier").sum(), (taxonomy["group"] == "K08356-negative").sum(), (taxonomy["group"] == "Other K08356 clade").sum(), len(stats), len(plot_stats)],
    }).to_csv(OUT / "analysis_summary.csv", index=False)

    print(pd.read_csv(OUT / "analysis_summary.csv").to_string(index=False))
    print("\nTop focused differences:")
    print(plot_stats[["category", "module", "mean_coverage_difference_pp", "mann_whitney_fdr", "prevalence_difference_pp", "fisher_fdr"]].head(15).to_string(index=False))


if __name__ == "__main__":
    main()

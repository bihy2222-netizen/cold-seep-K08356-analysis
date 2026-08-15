#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

from scipy.optimize import brentq
from scipy.stats import fisher_exact, nchypergeom_fisher


STRICT = "strict_synteny-supported_DIRM-like_IdrA-associated"
PARTIAL = "partial_IdrA-associated"


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else [])
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def taxonomy_parts(classification):
    values = {part.split("__", 1)[0]: part for part in classification.split(";") if "__" in part}
    return {"GTDB_domain": values.get("d", ""), "GTDB_phylum": values.get("p", ""),
            "GTDB_class": values.get("c", ""), "GTDB_order": values.get("o", ""),
            "GTDB_family": values.get("f", ""), "GTDB_genus": values.get("g", "")}


def load_gtdb(paths):
    result = {}
    for path in paths:
        if not path.is_file():
            continue
        for row in read_tsv(path):
            if row["user_genome"] in result:
                raise SystemExit(f'Duplicate GTDB user genome: {row["user_genome"]}')
            result[row["user_genome"]] = row
    return result


def exact_or_ci(table, alpha=0.05):
    a, b = table[0]
    c, d = table[1]
    odds, p_value = fisher_exact(table, alternative="two-sided")
    total, row1, col1 = a + b + c + d, a + b, a + c
    low_support = max(0, col1 - (total - row1))
    high_support = min(row1, col1)

    def survival(log_odds):
        theta = math.exp(log_odds)
        return nchypergeom_fisher.sf(a - 1, total, row1, col1, theta) - alpha / 2

    def cumulative(log_odds):
        theta = math.exp(log_odds)
        return nchypergeom_fisher.cdf(a, total, row1, col1, theta) - alpha / 2

    lower = 0.0 if a == low_support else math.exp(brentq(survival, -40, 40))
    upper = math.inf if a == high_support else math.exp(brentq(cumulative, -40, 40))
    return odds, p_value, lower, upper


def classified_family(row):
    family = row["GTDB_family"]
    return family if family and family != "f__" else "unclassified"


def fasta_metrics(path):
    lengths, current = [], 0
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    total = sum(lengths)
    cumulative, n50 = 0, 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= total / 2:
            n50 = length
            break
    return total, n50


def enrichment_rows(background, carrier_ids, analysis_id, samples=None, unclassified="exclude"):
    selected = [row for row in background if samples is None or row["sample"] in samples]
    if unclassified == "exclude":
        selected = [row for row in selected if classified_family(row) != "unclassified"]
    rhodo = [row for row in selected if classified_family(row) == "f__Rhodobacteraceae"]
    other = [row for row in selected if classified_family(row) != "f__Rhodobacteraceae"]
    table = [[sum(row["MAG_key"] in carrier_ids for row in rhodo),
              sum(row["MAG_key"] not in carrier_ids for row in rhodo)],
             [sum(row["MAG_key"] in carrier_ids for row in other),
              sum(row["MAG_key"] not in carrier_ids for row in other)]]
    odds, p_value, ci_low, ci_high = exact_or_ci(table)
    rhodo_total, other_total = sum(table[0]), sum(table[1])
    rhodo_prev = table[0][0] / rhodo_total if rhodo_total else math.nan
    other_prev = table[1][0] / other_total if other_total else math.nan
    return {
        "analysis_id": analysis_id, "unclassified_handling": unclassified,
        "samples_included": ";".join(sorted(samples)) if samples else "ALL",
        "background_MAGs_used": len(selected),
        "Rhodobacteraceae_carrier": table[0][0], "Rhodobacteraceae_noncarrier": table[0][1],
        "other_carrier": table[1][0], "other_noncarrier": table[1][1],
        "Fisher_two_sided_P": p_value, "odds_ratio": odds,
        "exact_CI95_low": ci_low, "exact_CI95_high": ci_high,
        "Rhodobacteraceae_carrier_prevalence": rhodo_prev,
        "outside_Rhodobacteraceae_carrier_prevalence": other_prev,
        "prevalence_difference": rhodo_prev - other_prev,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-crosswalk", required=True, type=Path)
    parser.add_argument("--gtdb-summary", required=True, nargs="+", type=Path)
    parser.add_argument("--candidate-classification", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    manifest = read_tsv(args.input_crosswalk)
    gtdb = load_gtdb(args.gtdb_summary)
    if len(gtdb) != 255:
        raise SystemExit(f"Expected 255 GTDB classifications, found {len(gtdb)}")

    candidates = read_tsv(args.candidate_classification)
    by_mag = defaultdict(list)
    for row in candidates:
        if row["bin_id"]:
            by_mag[(row["sample"], row["bin_id"])].append(row)

    background = []
    for source in manifest:
        gtdb_row = gtdb.get(source["gtdb_user_genome"])
        if not gtdb_row:
            raise SystemExit(f'Missing GTDB result: {source["gtdb_user_genome"]}')
        key = (source["sample"], source["bin_id"])
        linked = by_mag.get(key, [])
        taxonomy = taxonomy_parts(gtdb_row["classification"])
        genome_size, contig_n50 = fasta_metrics(source["bin_fasta"])
        background.append({
            "sample": source["sample"], "MAG_id": source["bin_id"],
            "MAG_key": f'{source["sample"]}|{source["bin_id"]}',
            "quality_class": source["quality_class"],
            "completeness": source["Completeness"], "contamination": source["Contamination"],
            "genome_size": genome_size, "contig_N50": contig_n50, **taxonomy,
            "classification_method": gtdb_row.get("classification_method", ""),
            "warnings": gtdb_row.get("warnings", ""),
            "strict_IdrA_carrier": "yes" if any(row["final_integrated_class"] == STRICT for row in linked) else "no",
            "partial_IdrA_carrier": "yes" if any(row["final_integrated_class"] == PARTIAL for row in linked) else "no",
            "candidate_ids": ";".join(row["protein_id"] for row in linked),
        })
    write_tsv(args.outdir / "01_all_MQHQ_MAG_GTDB_taxonomy.tsv", background)

    crosswalk_rows = []
    background_lookup = {(row["sample"], row["MAG_id"]): row for row in background}
    for row in candidates:
        mag = background_lookup.get((row["sample"], row["bin_id"]), {})
        crosswalk_rows.append({
            "candidate_id": row["protein_id"], "final_class": row["final_integrated_class"],
            "tree_review_status": "pending" if row["tree_manual_review_required"] == "yes" else "resolved_conservative",
            "sample": row["sample"], "contig_id": row["contig_id"], "bin_id": row["bin_id"],
            "MAG_id": row["bin_id"], "MAG_quality": mag.get("quality_class", ""),
            "GTDB_family": mag.get("GTDB_family", ""), "GTDB_genus": mag.get("GTDB_genus", ""),
        })
    write_tsv(args.outdir / "02_candidate_MAG_taxonomy_crosswalk.tsv", crosswalk_rows)

    strict_ids = {row["MAG_key"] for row in background if row["strict_IdrA_carrier"] == "yes"}
    extended_ids = {row["MAG_key"] for row in background
                    if row["strict_IdrA_carrier"] == "yes" or row["partial_IdrA_carrier"] == "yes"}
    main_rows = [enrichment_rows(background, strict_ids, "strict_bin_level", unclassified=mode)
                 for mode in ("exclude", "include_as_other")]
    write_tsv(args.outdir / "03_strict_carrier_2x2_table.tsv", main_rows)
    write_tsv(args.outdir / "04_strict_Fisher_OR_CI.tsv", main_rows)

    samples = sorted({row["sample"] for row in background})
    per_sample = [enrichment_rows(background, strict_ids, f"strict_{sample}", {sample}, mode)
                  for sample in samples for mode in ("exclude", "include_as_other")]
    write_tsv(args.outdir / "08_per_sample_enrichment.tsv", per_sample)
    leave_one_out = [enrichment_rows(background, strict_ids, f"strict_leave_out_{excluded}",
                                     set(samples) - {excluded}, mode)
                     for excluded in samples for mode in ("exclude", "include_as_other")]
    write_tsv(args.outdir / "09_leave_one_sample_out_enrichment.tsv", leave_one_out)
    extended = [enrichment_rows(background, extended_ids, "strict_plus_partial_bin_level", unclassified=mode)
                for mode in ("exclude", "include_as_other")]
    write_tsv(args.outdir / "10_strict_plus_partial_sensitivity.tsv", extended)

    report = [
        "# Methods and result boundary", "",
        "- Main carrier definition: strict synteny-supported DIRM-like IdrA-associated.",
        "- The 34 tree-review candidates are excluded from the strict carrier set.",
        "- Strict plus partial is a sensitivity analysis only.",
        "- These results are bin-level and must not be treated as independent species-level replicates.",
        "- ANI95 dereplicated statistics remain pending until cluster membership is supplied.",
        "- Family assignments are parsed only from GTDB-Tk r226 output.", "",
    ]
    (args.outdir / "11_methods_and_result_boundary.md").write_text("\n".join(report))


if __name__ == "__main__":
    main()

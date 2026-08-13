#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


FIELDS = ["qseqid", "qlen", "sseqid", "slen", "pident", "length", "mismatch", "gaps",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qframe"]
AMBIGUOUS_DELTA = 10.0


def reference_class(reference_id):
    prefixes = {
        "IdrA_HC": "strict_IdrA", "IdrA_PART": "partial_IdrA", "AioA_CAN": "AioA",
        "DMSOR_UNK": "unknown_DMSOR", "EXT_NapA": "NapA", "EXT_ArrA": "ArrA",
        "EXT_DmsA": "DmsA",
    }
    return next((label for prefix, label in prefixes.items() if reference_id.startswith(prefix)), "other")


def pair_id(read_id):
    return re.sub(r"(?:/|\s)[12]$", "", read_id.split()[0])


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def interval_text(intervals):
    return ";".join(f"{start}-{end}" for start, end in merge_intervals(intervals))


def breadth(intervals, length):
    return sum(end - start + 1 for start, end in merge_intervals(intervals)) / length if length else 0.0


def classify(best_family, delta, strict_intervals, strict_length):
    if not best_family:
        return "no IdrA protein-level support"
    if delta is not None and delta < AMBIGUOUS_DELTA:
        return "ambiguous DMSOR-family hit"
    if best_family != "strict_IdrA":
        return "reassigned to another DMSOR family"
    strict_breadth = breadth(strict_intervals, strict_length)
    regions = len(merge_intervals(strict_intervals))
    if strict_breadth >= 0.50 and regions >= 3:
        return "robust strict IdrA protein-level support"
    return "tentative strict IdrA protein-level support"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--total-reads", type=Path,
                        help="Optional TSV with sample, read_end, total_reads_searched")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    total_reads = {}
    if args.total_reads and args.total_reads.exists():
        with args.total_reads.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                total_reads[(row["sample"], row["read_end"])] = int(row["total_reads_searched"])

    all_hits = defaultdict(list)
    sample_end_hit_ids = defaultdict(set)
    per_reference = defaultdict(lambda: {"reads": set(), "pairs": set(), "alignments": 0,
                                         "intervals": [], "slen": 0, "bitscores": [], "identities": []})
    for path in sorted(args.hits_dir.glob("*.blastx.tsv.gz")):
        sample, mate = path.name.split(".", 2)[:2]
        with gzip.open(path, "rt") as handle:
            for values in csv.reader(handle, delimiter="\t"):
                if len(values) != len(FIELDS):
                    continue
                row = dict(zip(FIELDS, values))
                row.update(sample=sample, read_end=mate, pair_id=pair_id(row["qseqid"]),
                           family=reference_class(row["sseqid"]), bitscore_f=float(row["bitscore"]),
                           pident_f=float(row["pident"]), length_i=int(row["length"]),
                           slen_i=int(row["slen"]))
                row["sint"] = tuple(sorted((int(row["sstart"]), int(row["send"]))))
                all_hits[(sample, mate, row["qseqid"])].append(row)
                sample_end_hit_ids[(sample, mate)].add(row["qseqid"])

    best_rows = []
    pair_families = defaultdict(dict)
    for (sample, mate, read_id), hits in sorted(all_hits.items()):
        family_best = {}
        for hit in hits:
            old = family_best.get(hit["family"])
            if old is None or hit["bitscore_f"] > old["bitscore_f"]:
                family_best[hit["family"]] = hit
        ranked = sorted(family_best.values(), key=lambda x: (-x["bitscore_f"], x["family"], x["sseqid"]))
        best = ranked[0]
        second = ranked[1] if len(ranked) > 1 else None
        delta = best["bitscore_f"] - second["bitscore_f"] if second else None
        assigned = "ambiguous_DMSOR" if delta is not None and delta < AMBIGUOUS_DELTA else best["family"]
        pair_families[(sample, best["pair_id"])][mate] = assigned
        best_rows.append({
            "sample": sample, "read_end": mate, "read_id": read_id, "pair_id": best["pair_id"],
            "best_family": best["family"], "best_reference": best["sseqid"],
            "bitscore": f"{best['bitscore_f']:.1f}",
            "second_best_family": second["family"] if second else "",
            "second_best_reference": second["sseqid"] if second else "",
            "second_best_bitscore": f"{second['bitscore_f']:.1f}" if second else "",
            "bitscore_delta": f"{delta:.1f}" if delta is not None else "",
            "percent_identity": f"{best['pident_f']:.3f}", "alignment_length_aa": best["length_i"],
            "reference_start": best["sint"][0], "reference_end": best["sint"][1],
            "assigned_family": assigned,
            "conserved_region_only": "not_assessed_no_domain_coordinates",
        })
        key = (sample, best["sseqid"], best["family"])
        stats = per_reference[key]
        stats["reads"].add((mate, read_id)); stats["pairs"].add(best["pair_id"])
        stats["alignments"] += 1; stats["intervals"].append(best["sint"]); stats["slen"] = best["slen_i"]
        stats["bitscores"].append(best["bitscore_f"]); stats["identities"].append(best["pident_f"])

    read_out = args.out_dir / "diamond_blastx_best_hit_per_read.tsv"
    read_fields = ["sample", "read_end", "read_id", "pair_id", "best_family", "best_reference", "bitscore",
                   "second_best_family", "second_best_reference", "second_best_bitscore", "bitscore_delta",
                   "percent_identity", "alignment_length_aa", "reference_start", "reference_end",
                   "assigned_family", "conserved_region_only"]
    with read_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=read_fields, delimiter="\t")
        writer.writeheader(); writer.writerows(best_rows)

    ref_out = args.out_dir / "diamond_blastx_per_reference.tsv"
    ref_fields = ["sample", "reference_id", "reference_class", "unique_reads", "read_pairs", "alignments",
                  "protein_length_aa", "reference_coverage_breadth", "covered_intervals", "max_bitscore",
                  "mean_identity_pct", "conserved_region_only"]
    with ref_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ref_fields, delimiter="\t"); writer.writeheader()
        for (sample, ref, family), stats in sorted(per_reference.items()):
            writer.writerow({
                "sample": sample, "reference_id": ref, "reference_class": family,
                "unique_reads": len(stats["reads"]), "read_pairs": len(stats["pairs"]),
                "alignments": stats["alignments"], "protein_length_aa": stats["slen"],
                "reference_coverage_breadth": f"{breadth(stats['intervals'], stats['slen']):.6f}",
                "covered_intervals": interval_text(stats["intervals"]),
                "max_bitscore": f"{max(stats['bitscores']):.1f}",
                "mean_identity_pct": f"{sum(stats['identities']) / len(stats['identities']):.3f}",
                "conserved_region_only": "not_assessed_no_domain_coordinates",
            })

    summary_out = args.out_dir / "diamond_blastx_sample_summary.tsv"
    samples = sorted({sample for sample, _, _ in per_reference} | {sample for sample, _ in sample_end_hit_ids})
    summary_fields = ["sample", "read_end", "total_reads_searched", "DMSOR_hit_reads", "unique_read_ids",
                      "DMSOR_hit_read_pairs", "concordant_family_read_pairs", "best_family", "best_reference",
                      "bitscore", "second_best_family", "bitscore_delta", "percent_identity",
                      "alignment_length_aa", "reference_coverage_breadth", "covered_intervals",
                      "conserved_region_only", "final_classification"]
    with summary_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t"); writer.writeheader()
        for sample in samples:
            sample_rows = [r for r in best_rows if r["sample"] == sample]
            ranked = sorted(sample_rows, key=lambda r: -float(r["bitscore"]))
            best = ranked[0] if ranked else None
            family_scores = {}
            for row in sample_rows:
                family_scores[row["best_family"]] = max(family_scores.get(row["best_family"], 0.0), float(row["bitscore"]))
            fam_rank = sorted(family_scores.items(), key=lambda x: -x[1])
            best_family = fam_rank[0][0] if fam_rank else ""
            delta = fam_rank[0][1] - fam_rank[1][1] if len(fam_rank) > 1 else None
            strict_stats = [s for (ss, _, fam), s in per_reference.items() if ss == sample and fam == "strict_IdrA"]
            strict_intervals = [iv for s in strict_stats for iv in s["intervals"]]
            strict_length = max((s["slen"] for s in strict_stats), default=0)
            pairs = {r["pair_id"] for r in sample_rows}
            concordant = sum(1 for (ss, _), ends in pair_families.items()
                             if ss == sample and {"R1", "R2"} <= set(ends) and ends["R1"] == ends["R2"])
            ref_stats = per_reference.get((sample, best["best_reference"], best["best_family"])) if best else None
            for mate in ("R1", "R2", "combined"):
                mate_rows = sample_rows if mate == "combined" else [r for r in sample_rows if r["read_end"] == mate]
                writer.writerow({
                    "sample": sample, "read_end": mate,
                    "total_reads_searched": (sum(total_reads.get((sample, end), 0) for end in ("R1", "R2"))
                                               if mate == "combined" else total_reads.get((sample, mate), "")),
                    "DMSOR_hit_reads": len(mate_rows), "unique_read_ids": len({r["read_id"] for r in mate_rows}),
                    "DMSOR_hit_read_pairs": len(pairs) if mate == "combined" else "",
                    "concordant_family_read_pairs": concordant if mate == "combined" else "",
                    "best_family": best_family, "best_reference": best["best_reference"] if best else "",
                    "bitscore": best["bitscore"] if best else "",
                    "second_best_family": fam_rank[1][0] if len(fam_rank) > 1 else "",
                    "bitscore_delta": f"{delta:.1f}" if delta is not None else "",
                    "percent_identity": best["percent_identity"] if best else "",
                    "alignment_length_aa": best["alignment_length_aa"] if best else "",
                    "reference_coverage_breadth": (f"{breadth(ref_stats['intervals'], ref_stats['slen']):.6f}"
                                                   if ref_stats else ""),
                    "covered_intervals": interval_text(ref_stats["intervals"]) if ref_stats else "",
                    "conserved_region_only": "not_assessed_no_domain_coordinates",
                    "final_classification": classify(best_family, delta, strict_intervals, strict_length),
                })


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Annotate reviewed neighborhood ORFs using validated IdrB/P-like references."""

import argparse
import csv
from collections import defaultdict


def family(reference):
    if reference.startswith(("IdrP_like", "IdrP1", "IdrP2")): return "P_like"
    if reference.startswith(("IdrB", "AioB")): return "IdrB_related"
    return "other"


def accepted(group, evalue, bitscore, query_coverage, subject_coverage):
    if group == "IdrB_related":
        return evalue <= 1e-5 and bitscore >= 50 and query_coverage >= 50 and subject_coverage >= 35
    if group == "P_like":
        return evalue <= 1e-10 and bitscore >= 80 and query_coverage >= 45 and subject_coverage >= 35
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--neighborhood", required=True); parser.add_argument("--diamond", required=True)
    parser.add_argument("--out-neighborhood", required=True); parser.add_argument("--out-hits", required=True)
    args = parser.parse_args()
    with open(args.neighborhood, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    best = {}
    hit_rows = []
    with open(args.diamond, newline="") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if len(fields) < 14: continue
            query, subject = fields[0], fields[1]; alignment_length = float(fields[3])
            evalue, bitscore = float(fields[10]), float(fields[11]); qlen, slen = float(fields[12]), float(fields[13])
            group = family(subject); qcov = 100*alignment_length/qlen if qlen else 0
            scov = 100*alignment_length/slen if slen else 0
            keep = accepted(group, evalue, bitscore, qcov, scov)
            hit = {"query_id": query, "reference_id": subject, "family": group, "evalue": evalue,
                   "bitscore": bitscore, "query_coverage_pct": round(qcov, 3),
                   "subject_coverage_pct": round(scov, 3), "accepted": "yes" if keep else "no"}
            hit_rows.append(hit)
            if keep and (query not in best or bitscore > best[query]["bitscore"]): best[query] = hit
    for row in rows:
        hit = best.get(row["neighbor_id"])
        if hit:
            row["neighbor_family"] = hit["family"]; row["neighbor_annotation_accepted"] = "yes"
            row["neighbor_bitscore"] = hit["bitscore"]; row["neighbor_reference"] = hit["reference_id"]
    with open(args.out_neighborhood, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    fields = list(hit_rows[0]) if hit_rows else ["query_id","reference_id","family","evalue","bitscore","query_coverage_pct","subject_coverage_pct","accepted"]
    with open(args.out_hits, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(hit_rows)


if __name__ == "__main__": main()


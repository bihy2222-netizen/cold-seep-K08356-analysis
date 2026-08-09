#!/usr/bin/env python3
import csv
import os
from pathlib import Path


ROOT = Path(os.environ.get("ROOT", "/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra"))
HMM_DIR = ROOT / "09_hmmsearch"
SUMMARY = ROOT / "11_summary"


def parse_tblout(path):
    rows = []
    run = path.name.split(".")[0]
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 19:
                continue
            rows.append({
                "run": run,
                "protein_id": parts[0],
                "hmm_name": parts[2],
                "full_evalue": parts[4],
                "full_bitscore": parts[5],
                "full_bias": parts[6],
                "best_domain_evalue": parts[7],
                "best_domain_bitscore": parts[8],
                "description": " ".join(parts[18:]) if len(parts) > 18 else "",
                "candidate_tier": tier(float(parts[5]), float(parts[4])),
            })
    return rows


def tier(bitscore, evalue):
    threshold = os.environ.get("IDRA_BITSCORE_THRESHOLD")
    if threshold:
        try:
            if bitscore >= float(threshold):
                return "Tier_A_calibrated_threshold"
        except ValueError:
            pass
    if evalue <= 1e-20 and bitscore >= 100:
        return "Tier_B_significant_HMM_no_calibrated_threshold"
    if evalue <= 1e-5:
        return "Tier_C_weak_manual_review"
    return "below_reporting_threshold"


def main():
    SUMMARY.mkdir(parents=True, exist_ok=True)
    rows = []
    for tbl in sorted(HMM_DIR.glob("*.idra.tbl")):
        rows.extend(parse_tblout(tbl))
    fields = [
        "run", "protein_id", "hmm_name", "full_evalue", "full_bitscore",
        "full_bias", "best_domain_evalue", "best_domain_bitscore",
        "candidate_tier", "description",
    ]
    out = SUMMARY / "idra_HMM_transcript_candidates.tsv"
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out} rows={len(rows)}")


if __name__ == "__main__":
    main()

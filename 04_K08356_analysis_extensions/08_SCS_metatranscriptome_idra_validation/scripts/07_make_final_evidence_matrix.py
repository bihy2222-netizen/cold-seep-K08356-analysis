#!/usr/bin/env python3
import csv
import os
from collections import defaultdict
from pathlib import Path


ROOT = Path(os.environ.get("ROOT", "/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra"))
SUMMARY = ROOT / "11_summary"


def read_tsv(path):
    if not path.exists():
        return []
    with open(path, newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    expr = read_tsv(SUMMARY / "idra_bowtie_expression.tsv")
    hmm = read_tsv(SUMMARY / "idra_HMM_transcript_candidates.tsv")
    hmm_by_run = defaultdict(list)
    for row in hmm:
        hmm_by_run[row["run"]].append(row)

    out_rows = []
    for row in expr:
        if row.get("reference_panel") not in {"strict", "extended"}:
            continue
        level = row.get("expression_evidence_level", "Not detected")
        has_hmm = "YES" if hmm_by_run.get(row["run"]) else "NO"
        final = level
        if level in {"High confidence", "Moderate confidence"} and has_hmm == "YES":
            final = level + " with assembly-HMM support"
        out = {
            "sample": row.get("sample", ""),
            "run": row.get("run", ""),
            "reference_panel": row.get("reference_panel", ""),
            "reference_gene": row.get("reference_gene", ""),
            "IdrA_FPM": row.get("IdrA_FPM", ""),
            "IdrA_RPKM": row.get("IdrA_RPKM", ""),
            "MAPQ20_fragments_depth_approx": row.get("MAPQ20_fragments_depth_approx", ""),
            "mean_depth": row.get("mean_depth", ""),
            "breadth_1x": row.get("breadth_1x", ""),
            "breadth_3x": row.get("breadth_3x", ""),
            "regions_with_coverage": row.get("regions_with_coverage", ""),
            "assembly_HMM_hit_in_sample": has_hmm,
            "final_evidence_level": final,
            "notes": "Bowtie2 mapped to nucleotide IdrA CDS panel; not HMM-to-read mapping; FPM/RPKM are panel-target expression metrics.",
        }
        out_rows.append(out)

    fields = [
        "sample", "run", "reference_panel", "reference_gene", "IdrA_FPM",
        "IdrA_RPKM", "MAPQ20_fragments_depth_approx", "mean_depth",
        "breadth_1x", "breadth_3x", "regions_with_coverage",
        "assembly_HMM_hit_in_sample", "final_evidence_level", "notes",
    ]
    out = SUMMARY / "idra_final_evidence_matrix.tsv"
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    positives = [r for r in out_rows if not r["final_evidence_level"].startswith("Not detected")]
    with open(SUMMARY / "idra_positive_samples.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(positives)
    with open(SUMMARY / "README_results_CN.md", "w") as handle:
        handle.write("# 南海冷泉宏转录组 IdrA 表达验证结果说明\n\n")
        handle.write("本流程使用两条证据路线：RNA reads 回贴到核酸 IdrA CDS 参考，以及宏转录组组装蛋白后的 idra.hmm 搜索。\n\n")
        handle.write("注意：idra.hmm 是蛋白 profile-HMM，不能直接作为 RNA reads 的 Bowtie2 mapping 参考。\n\n")
        handle.write(f"- final evidence rows: {len(out_rows)}\n")
        handle.write(f"- positive rows: {len(positives)}\n")
        handle.write("- Bowtie2 面板结果报告 IdrA_FPM/RPKM/覆盖度，不等同于全转录组 TPM。\n")
        handle.write("- 所有阳性仍需系统发育、竞争 DMSOR 参考和邻域证据确认后才能称为 DIRM-like IdrA。\n")
    print(f"wrote {out} rows={len(out_rows)} positives={len(positives)}")


if __name__ == "__main__":
    main()

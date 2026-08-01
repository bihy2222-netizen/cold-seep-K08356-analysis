#!/usr/bin/env python3
"""Validate cardinalities and completion markers in a METABOLIC-G output."""

import argparse
import csv
from pathlib import Path


def count_files(path, pattern):
    return sum(item.is_file() for item in path.glob(pattern)) if path.is_dir() else 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-mag-count", type=int)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    input_count = count_files(args.input_dir, "*.faa")
    expected = args.expected_mag_count if args.expected_mag_count is not None else input_count
    run_log = args.output_dir / "METABOLIC_run.log"
    metabolic_log = args.output_dir / "METABOLIC_log.log"
    workbook = args.output_dir / "METABOLIC_result.xlsx"
    kegg_count = count_files(args.output_dir / "KEGG_identifier_result", "*")
    pdf_count = count_files(
        args.output_dir / "METABOLIC_Figures" / "Nutrient_Cycling_Diagrams",
        "*.pdf",
    )
    spreadsheet_count = count_files(args.output_dir / "METABOLIC_result_each_spreadsheet", "*")
    completed = metabolic_log.is_file() and "METABOLIC-G was done" in metabolic_log.read_text(errors="ignore")

    checks = {
        "input_MAG_count": (input_count, input_count == expected),
        "KEGG_identifier_file_count": (kegg_count, kegg_count == expected * 2),
        "nutrient_cycle_PDF_count": (pdf_count, pdf_count == expected * 4),
        "per_category_spreadsheet_count": (spreadsheet_count, spreadsheet_count == 6),
        "final_workbook_bytes": (workbook.stat().st_size if workbook.is_file() else 0, workbook.is_file() and workbook.stat().st_size > 0),
        "METABOLIC_run_log_present": (str(run_log.is_file()).lower(), run_log.is_file()),
        "completion_marker_present": (str(completed).lower(), completed),
    }
    status = "pass" if all(passed for _, passed in checks.values()) else "fail"

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "observed", "expected_or_rule", "status"])
        expected_rules = {
            "input_MAG_count": str(expected),
            "KEGG_identifier_file_count": str(expected * 2),
            "nutrient_cycle_PDF_count": str(expected * 4),
            "per_category_spreadsheet_count": "6",
            "final_workbook_bytes": ">0",
            "METABOLIC_run_log_present": "true",
            "completion_marker_present": "true",
        }
        for metric, (observed, passed) in checks.items():
            writer.writerow([metric, observed, expected_rules[metric], "pass" if passed else "fail"])
        writer.writerow(["overall_validation", status, "pass", status])

    for metric, (observed, passed) in checks.items():
        print(f"{metric}={observed}\t{'pass' if passed else 'fail'}")
    print(f"overall_validation={status}")
    print(f"summary={args.summary}")
    if args.strict and status != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

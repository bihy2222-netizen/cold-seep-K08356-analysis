#!/usr/bin/env python3
import csv
import os
import re
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(os.environ.get("ROOT", "/home/ps/ps1/data/bihongyu/cold_seep/public_validation/SCS_metatranscriptome_idra"))
METADATA = ROOT / "00_metadata"
PROJECTS = ["PRJNA738468", "PRJNA739005", "PRJNA831433"]

TARGET_RE = re.compile(r"QDN|Qiongdongnan|W01B|W03B|W04B|W07|Lingshui|JL[_-]?|Jiaolong|Haima|HM[_-]?|S11[_-]?", re.I)
QDN_RE = re.compile(r"QDN|Qiongdongnan|W01B|W03B|W04B|W07|Lingshui", re.I)
JL_RE = re.compile(r"JL[_-]?|Jiaolong", re.I)
HAIMA_RE = re.compile(r"Haima|HM[_-]?|S11[_-]?", re.I)

ENA_FIELDS = [
    "run_accession", "sample_accession", "experiment_accession", "study_accession",
    "secondary_study_accession", "sample_title", "scientific_name",
    "library_strategy", "library_source", "library_selection", "library_layout",
    "instrument_platform", "instrument_model", "fastq_ftp", "fastq_md5",
    "fastq_bytes", "read_count", "base_count",
]


def runinfo_with_edirect(project, out_path):
    esearch = os.environ.get("ESEARCH", "esearch")
    efetch = os.environ.get("EFETCH", "efetch")
    cmd = f'{esearch} -db sra -query "{project}" | {efetch} -format runinfo'
    with open(out_path, "w") as handle:
        subprocess.run(cmd, shell=True, check=True, stdout=handle)


def read_csv(path):
    with open(path, newline="", errors="replace") as handle:
        return list(csv.DictReader(handle))


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fetch_ena(project):
    query = {
        "accession": project,
        "result": "read_run",
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "download": "true",
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(query)
    text = urllib.request.urlopen(url, timeout=180).read().decode("utf-8", errors="replace")
    (METADATA / f"{project}.ena_read_run.tsv").write_text(text)
    return {row["run_accession"]: row for row in csv.DictReader(text.splitlines(), delimiter="\t")}


def split_semicolon(value):
    return [x for x in str(value or "").split(";") if x]


def total_bytes(ena_row):
    total = 0
    for item in split_semicolon(ena_row.get("fastq_bytes", "")):
        try:
            total += int(item)
        except ValueError:
            pass
    return total


def joined(row):
    return " ".join(str(v) for v in row.values())


def is_rna(row):
    text = " ".join(str(row.get(k, "")) for k in ["LibraryStrategy", "LibrarySource", "LibrarySelection", "ENA_library_strategy", "ENA_library_source", "ENA_library_selection"])
    return bool(re.search(r"RNA|TRANSCRIPTOMIC|METATRANSCRIPTOMIC", text, re.I))


def target_group(row):
    text = joined(row)
    if QDN_RE.search(text):
        return "Priority_1_QDN"
    if JL_RE.search(text):
        return "Priority_3_Jiaolong"
    if HAIMA_RE.search(text):
        return "Priority_4_Haima"
    if TARGET_RE.search(text):
        return "Priority_5_other_SCS"
    return ""


def exclusion_reason(row):
    if not is_rna(row):
        return "not_RNA_metatranscriptome"
    if not target_group(row):
        return "not_SCS_target_keyword"
    text = joined(row)
    if re.search(r"16S|18S|ITS|AMPLICON", text, re.I):
        return "amplicon"
    if re.search(r"WGS|METAGENOMIC", str(row.get("LibraryStrategy", "")), re.I) and not re.search(r"RNA|TRANSCRIPTOMIC", text, re.I):
        return "DNA_WGS_metagenome"
    return ""


def main():
    METADATA.mkdir(parents=True, exist_ok=True)
    all_rows, all_rna, targets, excluded = [], [], [], []

    for project in PROJECTS:
        runinfo_path = METADATA / f"{project}.runinfo.csv"
        if not runinfo_path.exists() or runinfo_path.stat().st_size == 0:
            runinfo_with_edirect(project, runinfo_path)
        ena = fetch_ena(project)
        project_rows = read_csv(runinfo_path)
        for row in project_rows:
            run = row.get("Run", "")
            e = ena.get(run, {})
            row["project"] = project
            for key, value in e.items():
                row[f"ENA_{key}"] = value
            urls = ["https://" + x for x in split_semicolon(e.get("fastq_ftp", ""))]
            row["FASTQ_URL"] = ";".join(urls)
            row["FASTQ_MD5"] = e.get("fastq_md5", "")
            row["FASTQ_BYTES"] = e.get("fastq_bytes", "")
            row["FASTQ_TOTAL_BYTES"] = str(total_bytes(e))
            row["FASTQ_TOTAL_GiB"] = f"{total_bytes(e) / 1024**3:.3f}"
            row["target_group"] = target_group(row)
            row["inferred_is_RNA"] = "YES" if is_rna(row) else "NO"
            reason = exclusion_reason(row)
            if reason:
                ex = dict(row)
                ex["exclude_reason"] = reason
                excluded.append(ex)
            else:
                targets.append(row)
            if is_rna(row):
                all_rna.append(row)
            all_rows.append(row)

        fields = list(dict.fromkeys(k for r in project_rows for k in r.keys()))
        write_tsv(METADATA / f"{project}.enriched_runinfo.tsv", project_rows, fields)

    fields = list(dict.fromkeys(k for r in all_rows for k in r.keys()))
    write_tsv(METADATA / "all_runs.tsv", all_rows, fields)
    write_tsv(METADATA / "all_RNA_runs.tsv", all_rna, fields)
    write_tsv(METADATA / "SCS_sediment_metatranscriptome_runs.tsv", targets, fields)
    write_tsv(METADATA / "excluded_runs.tsv", excluded, list(dict.fromkeys(k for r in excluded for k in r.keys())) if excluded else fields + ["exclude_reason"])

    qdn = [r for r in targets if r["target_group"] == "Priority_1_QDN"]
    qdn.sort(key=lambda r: (r.get("SampleName", ""), r.get("Run", "")))
    write_tsv(METADATA / "QDN_priority_runs.tsv", qdn, fields)
    write_tsv(METADATA / "download_manifest.tsv", qdn, fields)
    write_tsv(METADATA / "PRJNA739005_QDN_RNA_fastq_manifest.tsv", [r for r in qdn if r["project"] == "PRJNA739005"], fields)

    with open(METADATA / "runinfo_summary.tsv", "w") as handle:
        handle.write("project\tall_runs\tRNA_runs\ttarget_runs\tQDN_runs\ttarget_FASTQ_GiB\tQDN_FASTQ_GiB\n")
        for project in PROJECTS:
            pr = [r for r in all_rows if r["project"] == project]
            rr = [r for r in all_rna if r["project"] == project]
            tr = [r for r in targets if r["project"] == project]
            qr = [r for r in qdn if r["project"] == project]
            tg = sum(int(r.get("FASTQ_TOTAL_BYTES") or 0) for r in tr) / 1024**3
            qg = sum(int(r.get("FASTQ_TOTAL_BYTES") or 0) for r in qr) / 1024**3
            handle.write(f"{project}\t{len(pr)}\t{len(rr)}\t{len(tr)}\t{len(qr)}\t{tg:.3f}\t{qg:.3f}\n")

    print((METADATA / "runinfo_summary.tsv").read_text())
    print("QDN priority runs:")
    for row in qdn:
        print("\t".join([row.get("project", ""), row.get("Run", ""), row.get("SampleName", ""), row.get("LibraryStrategy", ""), row.get("LibrarySource", ""), row.get("LibraryLayout", ""), row.get("FASTQ_TOTAL_GiB", "")]))


if __name__ == "__main__":
    main()

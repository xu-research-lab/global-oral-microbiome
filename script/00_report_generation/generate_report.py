#!/usr/bin/env python3
"""Generate the Global Oral Microbiome analysis report.

This script does not rerun any biological analysis. It renders the Markdown
report from a maintained template, updates final metadata counts from the
current TSV files, and checks whether embedded figure paths resolve locally.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path("/beegfs/huangxiaochang/global_oral")
REPORT_DIR = ROOT / "analysis_report"
TEMPLATE = REPORT_DIR / "global_oral_analysis_report.template.md"
OUTPUT = REPORT_DIR / "global_oral_analysis_report.md"

TABLES = {
    "raw": ROOT / "new_analysis/metadata_all_for_qc_raw.tsv",
    "all_qc": ROOT / "new_analysis/metadata_all_qc.tsv",
    "oral_qc": ROOT / "new_analysis/metadata_oral_qc.tsv",
    "oral_dedup": ROOT / "new_analysis/metadata_oral_qc_dedup.tsv",
}

LAYER5_TOP10 = (
    ROOT
    / "ASV_vs_HOMD_analysis/containment_derep/gg2_backbone_suspicious_filter/funnel_layer5_top10.tsv"
)


def fmt_int(value: int) -> str:
    return f"{value:,}"


def count_table(path: Path) -> dict[str, int]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = 0
        studies: set[str] = set()
        subjects: set[str] = set()
        for row in reader:
            rows += 1
            bp = (row.get("bioproject_id") or row.get("bioproject") or "").strip()
            subject = (row.get("subject_id") or "").strip()
            if bp and bp != "NA":
                studies.add(bp)
            if subject and subject != "NA":
                subjects.add(subject)
    return {"samples": rows, "studies": len(studies), "subjects": len(subjects)}


def markdown_table_from_tsv(path: Path, max_rows: int = 10) -> str:
    keep_cols = [
        "derep_id",
        "subject_prevalence",
        "n_bioprojects",
        "total_reads",
        "sites",
        "subject_geo_regions",
        "gg2_phylum",
        "gg2_genus",
        "gg2_species",
    ]
    labels = {
        "derep_id": "ASV",
        "subject_prevalence": "Subjects",
        "n_bioprojects": "BioProjects",
        "total_reads": "Reads",
        "sites": "Sites",
        "subject_geo_regions": "Geographies",
        "gg2_phylum": "GG2 phylum",
        "gg2_genus": "GG2 genus",
        "gg2_species": "GG2 species",
    }
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)[:max_rows]

    def clean(value: str) -> str:
        value = (value or "").replace("|", "\|")
        if ";" in value:
            value = value.replace(";", "; ")
        return value

    header = "| " + " | ".join(labels[col] for col in keep_cols) + " |"
    sep = "|" + "|".join("---" for _ in keep_cols) + "|"
    body = []
    for row in rows:
        vals = []
        for col in keep_cols:
            val = clean(row.get(col, ""))
            if col in {"subject_prevalence", "n_bioprojects", "total_reads"}:
                try:
                    val = fmt_int(int(float(val)))
                except ValueError:
                    pass
            vals.append(val)
        body.append("| " + " | ".join(vals) + " |")
    return "\n".join([header, sep, *body])


def render() -> str:
    text = TEMPLATE.read_text()
    replacements: dict[str, str] = {}
    for key, path in TABLES.items():
        counts = count_table(path)
        for metric, value in counts.items():
            replacements[f"{{{{{key}_{metric}}}}}"] = fmt_int(value)
    replacements["{{layer5_top10_table}}"] = markdown_table_from_tsv(LAYER5_TOP10)

    for token, value in replacements.items():
        text = text.replace(token, value)
    return text


def check_embedded_figures(markdown: str) -> list[str]:
    missing: list[str] = []
    for match in re.finditer(r"!\[[^\]]*\]\(([^)]+)\)", markdown):
        fig_path = match.group(1)
        if fig_path.startswith(("http://", "https://")):
            continue
        resolved = (REPORT_DIR / fig_path).resolve()
        if not resolved.exists():
            missing.append(fig_path)
    return missing


def main() -> None:
    markdown = render()
    missing = check_embedded_figures(markdown)
    OUTPUT.write_text(markdown)
    print(f"Wrote {OUTPUT}")
    if missing:
        print("Missing embedded figure paths:")
        for path in missing:
            print(f"  - {path}")
    else:
        print("All embedded local figure paths exist.")


if __name__ == "__main__":
    main()

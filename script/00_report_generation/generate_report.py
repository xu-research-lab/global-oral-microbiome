#!/usr/bin/env python3
"""Generate the Global Oral Microbiome analysis report.

This script does not rerun any biological analysis. It renders the Markdown
report from a maintained template, updates final metadata counts from the
current TSV files, and checks whether embedded figure paths resolve locally.

Two output files are produced:
  - global_oral_analysis_report.md  — path-referenced figures (compact, for local use)
  - global_oral_analysis_report_self_contained.md — base64-embedded figures (self-contained, for sharing)
"""
from __future__ import annotations

import base64
import csv
import mimetypes
import re
from pathlib import Path

ROOT = Path("/beegfs/huangxiaochang/global_oral")
REPORT_DIR = ROOT / "analysis_report"
TEMPLATE = REPORT_DIR / "global_oral_analysis_report.template.md"
OUTPUT = REPORT_DIR / "global_oral_analysis_report.md"
OUTPUT_SELF_CONTAINED = REPORT_DIR / "global_oral_analysis_report_self_contained.md"

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
        if fig_path.startswith(("http://", "https://", "data:")):
            continue
        resolved = (REPORT_DIR / fig_path).resolve()
        if not resolved.exists():
            missing.append(fig_path)
    return missing


# ---------------------------------------------------------------------------
# Base64 embedding for self-contained output
# ---------------------------------------------------------------------------

# Only cache path → base64 data URI; images can be referenced multiple times.
_EMBED_CACHE: dict[str, str] = {}


def _resolve_figure_path(relative_path: str) -> Path:
    """Resolve a figure path (relative to REPORT_DIR) to an absolute Path."""
    return (REPORT_DIR / relative_path).resolve()


def _encode_image_to_data_uri(abs_path: Path) -> str:
    """Read an image file and return a 'data:image/<type>;base64,…' string."""
    mime_type, _ = mimetypes.guess_type(str(abs_path))
    if mime_type is None:
        mime_type = "image/png"
    data = abs_path.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    return f"data:{mime_type};base64,{b64}"


def embed_figures(markdown: str) -> tuple[str, int, int]:
    """Replace local image paths with base64 data URIs.

    Returns (transformed_markdown, embed_count, missing_count).
    """
    embedded_count = 0
    missing_count = 0

    def _replacer(match: re.Match) -> str:
        nonlocal embedded_count, missing_count
        alt = match.group(1)
        rel_path = match.group(2)

        # Skip already-embedded or remote URLs
        if rel_path.startswith(("http://", "https://", "data:")):
            return match.group(0)

        abs_path = _resolve_figure_path(rel_path)

        if rel_path in _EMBED_CACHE:
            data_uri = _EMBED_CACHE[rel_path]
        elif abs_path.exists():
            data_uri = _encode_image_to_data_uri(abs_path)
            _EMBED_CACHE[rel_path] = data_uri
        else:
            missing_count += 1
            return match.group(0)  # keep original reference

        embedded_count += 1
        return f"![{alt}]({data_uri})"

    result = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", _replacer, markdown)
    return result, embedded_count, missing_count


# ---------------------------------------------------------------------------
def main() -> None:
    markdown = render()

    # --- path-referenced version ---
    missing = check_embedded_figures(markdown)
    OUTPUT.write_text(markdown)
    print(f"Wrote {OUTPUT}")
    if missing:
        print("Missing embedded figure paths:")
        for path in missing:
            print(f"  - {path}")
    else:
        print("All embedded local figure paths exist.")

    # --- self-contained (base64-embedded) version ---
    self_contained, embed_n, embed_missing = embed_figures(markdown)
    OUTPUT_SELF_CONTAINED.write_text(self_contained)
    size_mb = OUTPUT_SELF_CONTAINED.stat().st_size / (1024 * 1024)
    print(f"Wrote {OUTPUT_SELF_CONTAINED}  ({size_mb:.1f} MB, {embed_n} images embedded)")
    if embed_missing:
        print(f"  {embed_missing} images could not be found; paths kept as-is.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit raw per-project SRA metadata for negative / blank / mock control samples.

Curation (report section 1.1) removed control and non-human samples from the analysis
tables, so blanks are absent from the curated metadata but usually remain in the raw
`SraRunTable.csv` / original `metadata.csv`. Knowing which BioProjects carry blanks tells
us where blank-based decontamination (e.g. `decontam` prevalence mode) is possible, and
which blanks could be re-fetched to build a data-driven contaminant profile for the
low-HOMD "novel oral taxa" analysis.

This scanner is intentionally higher-precision than a bare "control" grep (which would
match healthy/case-control language): it categorizes each match (extraction blank, blank,
NTC/no-template, negative control, mock community, reagent/water) and reports per-project
counts. Matches still warrant manual confirmation — in particular, a column literally
valued "blank" (a missing-field placeholder) can inflate the plain `blank` category.

Input : data/projects/<BioProject>/{SraRunTable.csv, metadata.csv}
Output: results/tables/05_homd_gg2_closed_reference/negative_control_bioproject_audit.tsv
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PROJECTS = REPO_ROOT / "data/projects"
DEFAULT_OUT = REPO_ROOT / "results/tables/05_homd_gg2_closed_reference/negative_control_bioproject_audit.tsv"
RAW_NAMES = ("SraRunTable.csv", "metadata.csv")

CATEGORIES = [
    ("extraction_blank", re.compile(r"\bext(raction)?[\s_\-\.]*(blank|b\d|ctrl)|extblank|extractionblank", re.I)),
    ("blank",            re.compile(r"\bblank\b|\bblnk\b", re.I)),
    ("ntc_no_template",  re.compile(r"\bntc\b|no[\s\-_]*template", re.I)),
    ("negative_control", re.compile(r"negative[\s_\-]*control|neg[\s_\-\.]*ctrl|neg[\s_\-\.]*control", re.I)),
    ("mock_community",   re.compile(r"mock[\s_\-]*(community|comm|control)", re.I)),
    ("reagent_water",    re.compile(r"reagent[\s_\-]*(control|blank)|sterile[\s_\-]*water|water[\s_\-]*(control|blank)|no[\s\-_]*dna", re.I)),
]
STRONG = {"extraction_blank", "blank", "ntc_no_template"}


def scan_project(pdir: Path) -> Counter:
    files = [pdir / n for n in RAW_NAMES if (pdir / n).exists()]
    files += [Path(f) for f in glob.glob(str(pdir / "*[Ss]ra*un*able*.csv")) if Path(f) not in files]
    counts: Counter = Counter()
    for f in files:
        try:
            with open(f, newline="", encoding="utf-8", errors="ignore") as fh:
                for row in csv.reader(fh):
                    line = " ".join(row)
                    for cat, pat in CATEGORIES:
                        if pat.search(line):
                            counts[cat] += 1
                            break  # one category per row
        except Exception:
            continue
    return counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--projects-dir", type=Path, default=DEFAULT_PROJECTS)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    rows = []
    for p in sorted(os.listdir(args.projects_dir)):
        pdir = args.projects_dir / p
        if not pdir.is_dir():
            continue
        cc = scan_project(pdir)
        if cc:
            rows.append((p, sum(cc.values()), cc))
    rows.sort(key=lambda x: -x[1])

    cats = [c for c, _ in CATEGORIES]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["bioproject", "n_candidate_control_rows"] + cats)
        for p, n, cc in rows:
            w.writerow([p, n] + [cc.get(c, 0) for c in cats])

    strong = [r for r in rows if any(r[2].get(c, 0) for c in STRONG)]
    print(f"projects with candidate controls: {len(rows)}")
    print(f"  with explicit extraction-blank / blank / NTC (high confidence): {len(strong)}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
jaccard_heatmap.py — Genus-overlap (Jaccard) similarity between the healthy
network's modules (HC) and each disease network's modules (PD, DC, OSCC),
for the >30-threshold all-site networks built by this pipeline.

For every (HC module, disease module) pair,
  J = |genera in both| / |genera in either|
computed on taxon_id (greengene2 genus-level feature IDs). All four module
tables come from the SAME genus table, so the IDs are directly comparable;
the g__ prefix is stripped defensively anyway.

Outputs -> result_paired/jaccard/:
  jaccard_<dis>_vs_HC.csv     (full matrix, rows = HC modules)
  jaccard_summary.md          (best match per disease module + honest caveats)

Plotting is done by script_hc_pd_dc_oscc/plot_jaccard_facets.R (ggplot facet
style, mirroring result_disease/modularity/jaccard_vs_allnet_facets.png).

Usage (from network/):
  python script_hc_pd_dc_oscc/jaccard_heatmap.py
"""

import argparse
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MOD_DIR = os.path.join(ROOT, "result_paired", "modularity")
OUT_DIR = os.path.join(ROOT, "result_paired", "jaccard")
DISEASES = ["PD", "DC", "OSCC"]
DIS_LABEL = {"PD": "Periodontal disease", "DC": "Dental caries",
             "OSCC": "OSCC"}


def norm(ids):
    """genus-level feature IDs: strip a g__ prefix defensively."""
    return {str(t)[3:] if str(t).startswith("g__") else str(t) for t in ids}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diseases", default=",".join(DISEASES),
                    help="comma-separated subset of disease networks "
                         "(default: PD,DC,OSCC)")
    ap.add_argument("--suffix", default="",
                    help="suffix for the summary file when a subset is used "
                         "(e.g. noDC)")
    args = ap.parse_args()
    diseases = [d.strip().upper() for d in args.diseases.split(",") if d.strip()]
    bad = [d for d in diseases if d not in DISEASES]
    if bad:
        raise SystemExit(f"unknown diseases: {bad} (choose from {DISEASES})")
    os.makedirs(OUT_DIR, exist_ok=True)
    hm = pd.read_csv(os.path.join(MOD_DIR, "modules_HC.csv"),
                     dtype={"taxon_id": str})
    healthy = {m: norm(g.taxon_id) for m, g in hm.groupby("module")}

    dis_txt = " / ".join(diseases)
    summary = [f"# Module Jaccard similarity, HC / {dis_txt} networks",
               "",
               "Healthy network = all-site HC network from this pipeline; "
               f"disease networks = the {dis_txt} all-site networks.",
               "Jaccard = shared genera between two modules / union genera, "
               "computed on genus-level feature IDs (taxon_id).",
               ""]
    for dis in diseases:
        dm = pd.read_csv(os.path.join(MOD_DIR, f"modules_{dis}.csv"),
                         dtype={"taxon_id": str})
        disease = {m: norm(g.taxon_id) for m, g in dm.groupby("module")}

        rows = []
        for hm_m, hs in sorted(healthy.items()):
            rows.append([len(hs & ds) / len(hs | ds)
                         for _, ds in sorted(disease.items())])
        jac = pd.DataFrame(rows, index=[f"HC-M{m}" for m in sorted(healthy)],
                           columns=[f"{dis}-M{m}" for m in sorted(disease)])
        jac.to_csv(os.path.join(OUT_DIR, f"jaccard_{dis}_vs_HC.csv"))

        # summary lines: best match per disease module
        summary += [f"## {DIS_LABEL[dis]} ({dis})",
                    "",
                    "| Disease module | Best-matching healthy module | Jaccard | shared/union |",
                    "|---|---|---|---|"]
        for m, ds in sorted(disease.items()):
            jvals = {hm_m: len(hs & ds) / len(hs | ds)
                     for hm_m, hs in healthy.items()}
            best = max(jvals, key=jvals.get)
            shared = len(healthy[best] & ds)
            summary.append(
                f"| {dis}-M{m} | HC-M{best} | {jvals[best]:.2f} "
                f"| {shared}/{len(healthy[best] | ds)} |")
        summary.append("")

    # cross-network caveat
    summary += ["---", "",
                "**How to read this (stated as-is)**:",
                "- Jaccard measures genus-composition overlap only; two modules "
                "may share a high Jaccard while their edges and pooled samples",
                "  are entirely different (the four networks are computed "
                "independently), so the same function cannot be inferred.",
                "- HC 113 cohorts vs PD 11 / DC 5 / OSCC 10 cohorts: the "
                "difference in genus-space coverage systematically depresses",
                "  Jaccard (the HC module genus pool is far larger).",
                "- DC has only 5 cohorts; its network has limited robustness "
                "and the overlap values are indicative only.",
                "- The best match is the per-row maximum; no permutation test "
                "was run, values are for manual reference."]
    summary_name = "jaccard_summary.md" if not args.suffix \
        else f"jaccard_summary_{args.suffix}.md"
    with open(os.path.join(OUT_DIR, summary_name),
              "w", encoding="utf-8") as fh:
        fh.write("\n".join(summary))
    print(f"[done] -> {os.path.join(OUT_DIR, summary_name)}")


if __name__ == "__main__":
    main()
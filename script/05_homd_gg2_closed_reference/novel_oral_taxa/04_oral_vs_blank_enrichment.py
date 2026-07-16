#!/usr/bin/env python3
"""Phase 2/3: oral-vs-blank enrichment and residence flag.

Residence-vs-contamination is tested against BLANKS, not gut (oral-gut co-abundant taxa
are real biology; oral-gut co-abundance argues FOR authenticity, so it must not gate
residence). For each candidate feature this computes, across the BioProjects that contain
negative controls (see 25_negative_control_audit.py), its prevalence and mean relative
abundance in oral samples vs blanks, an enrichment ratio, and a residence flag. This
complements the decontam prevalence model (step 03) with an explicit oral-over-blank
signal, and it also supplies within-oral abundance for the ecology criterion.

Inputs
------
--abundance   long-format TSV: feature_id, sample_id, count   (feature_id in candidate
              space, e.g. derep_id; caller maps source ASVs -> derep upstream)
--metadata    TSV: sample_id, bioproject, is_blank (bool), is_oral (bool)
--candidates  optional TSV with a 'derep_id' column to restrict to Phase-1 candidates

A feature passes residence if it is prevalent and abundant in oral samples AND clearly
depleted in blanks (either low blank prevalence or a high oral/blank enrichment ratio and
higher oral abundance). Thresholds are conservative defaults; tune per data.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C  # noqa: E402


def as_bool(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin({"true", "1", "yes", "t", "y"})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--abundance", required=True, help="long-format feature_id,sample_id,count TSV")
    ap.add_argument("--metadata", required=True, help="sample_id,bioproject,is_blank,is_oral TSV")
    ap.add_argument("--candidates", default=None, help="optional TSV with derep_id column to restrict")
    ap.add_argument("--outdir", default=os.path.join(C.OUT_ROOT, "04_oral_vs_blank"))
    ap.add_argument("--min-oral-prev", type=float, default=0.05, help="min prevalence in oral samples")
    ap.add_argument("--max-blank-prev", type=float, default=0.10, help="blank prevalence at/below which a taxon is clean")
    ap.add_argument("--min-enrich-ratio", type=float, default=5.0, help="oral/blank prevalence ratio if blank prev is high")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    md = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    md["is_blank"] = as_bool(md["is_blank"])
    md["is_oral"] = as_bool(md["is_oral"])
    # Only BioProjects that actually contain a blank inform this test.
    bp_with_blank = set(md.loc[md["is_blank"], "bioproject"].unique())
    md = md[md["bioproject"].isin(bp_with_blank)].copy()
    oral_samples = set(md.loc[md["is_oral"] & ~md["is_blank"], "sample_id"])
    blank_samples = set(md.loc[md["is_blank"], "sample_id"])
    n_oral, n_blank = len(oral_samples), len(blank_samples)
    if n_oral == 0 or n_blank == 0:
        raise SystemExit("no oral and/or blank samples in blank-bearing BioProjects")

    ab = pd.read_csv(args.abundance, sep="\t", low_memory=False)
    ab.columns = ["feature_id", "sample_id", "count"][:len(ab.columns)]
    ab["count"] = pd.to_numeric(ab["count"], errors="coerce").fillna(0)
    if args.candidates:
        keep = set(pd.read_csv(args.candidates, sep="\t")["derep_id"])
        ab = ab[ab["feature_id"].isin(keep)]

    # per-sample total for relative abundance
    totals = ab.groupby("sample_id")["count"].sum()
    ab = ab[ab["count"] > 0].copy()
    ab["relab"] = ab["count"] / ab["sample_id"].map(totals)
    ab["grp"] = np.where(ab["sample_id"].isin(blank_samples), "blank",
                         np.where(ab["sample_id"].isin(oral_samples), "oral", "other"))
    ab = ab[ab["grp"].isin(["oral", "blank"])]

    rows = []
    for fid, sub in ab.groupby("feature_id"):
        o = sub[sub["grp"] == "oral"]; b = sub[sub["grp"] == "blank"]
        oral_prev = o["sample_id"].nunique() / n_oral
        blank_prev = b["sample_id"].nunique() / n_blank
        oral_ab = float(o["relab"].mean()) if len(o) else 0.0
        blank_ab = float(b["relab"].mean()) if len(b) else 0.0
        ratio = oral_prev / blank_prev if blank_prev > 0 else np.inf
        residence = (
            oral_prev >= args.min_oral_prev
            and oral_ab > blank_ab
            and (blank_prev <= args.max_blank_prev or ratio >= args.min_enrich_ratio)
        )
        rows.append(dict(
            feature_id=fid, oral_prevalence=round(oral_prev, 4), blank_prevalence=round(blank_prev, 4),
            oral_mean_relab=round(oral_ab, 6), blank_mean_relab=round(blank_ab, 6),
            oral_over_blank_prev_ratio=round(ratio, 3) if np.isfinite(ratio) else "inf",
            residence_pass=bool(residence),
        ))
    res = pd.DataFrame(rows).sort_values("oral_prevalence", ascending=False)
    out = os.path.join(args.outdir, "oral_vs_blank_enrichment.tsv")
    res.to_csv(out, sep="\t", index=False)
    npass = int(res["residence_pass"].sum()) if len(res) else 0
    print(f"blank-bearing BioProjects: {len(bp_with_blank)}  oral samples: {n_oral}  blanks: {n_blank}")
    print(f"features evaluated: {len(res)}  residence_pass: {npass}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

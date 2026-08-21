#!/usr/bin/env python3
"""
pipeline_v2 Step 06b: Phylotype aggregation — report, don't adjudicate.

For every phylotype, aggregate the member ASV contamination annotations into
per-phylotype composition columns. NO phylotype is excluded here. Two derived
markers are added for the figures:

  clean_phylotype : n_flagged == 0  (no member ASV was contamination-flagged)
  support_score   : max_prevalence * n_bioprojects_total * n_regions
                    (the user's "strongest evidence" ranking; thresholds are set
                     later from the actual distribution, see diagnostic summary)

Outputs:
  - 06b_phylotypes_aggregated.tsv : all 23,088 phylotypes + composition + markers
  - 06b_diagnostic_clean.tsv      : distribution of clean phylotypes (for threshold setting)
  - 06b_summary.json
"""
import json
import os
from pathlib import Path

import pandas as pd
import numpy as np

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
log = lambda s: print(f"[06b_agg] {s}", flush=True)

# ============================================================
# 1. Load novelty-classified phylotypes + member-level annotations
# ============================================================
log("Loading novelty-classified phylotypes...")
phy = pd.read_csv(PIPE_V2 / "06_phylotypes_novelty_classified.tsv", sep="\t", low_memory=False)
log(f"  {len(phy):,} phylotypes")

log("Loading support-floor ASVs (member annotations)...")
asv = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False)
asv["derep_id"] = asv["derep_id"].astype(str)
for c in ["contamination_flag"]:
    if c in asv.columns:
        asv[c] = asv[c].astype(str).str.lower().isin(["true", "1", "t"])
asv_meta = asv[["derep_id", "contamination_flag"]]
asv_lookup = asv_meta.set_index("derep_id").to_dict("index")

# ============================================================
# 2. Aggregate member composition per phylotype
# ============================================================
log("Aggregating member composition...")

def parse_members(x):
    if pd.isna(x):
        return []
    return [m for m in str(x).split(";") if m]

rows = []
for _, r in phy.iterrows():
    members = parse_members(r.get("member_derep_ids"))
    n = len(members)
    n_flag = 0
    for m in members:
        info = asv_lookup.get(m, {})
        n_flag += int(info.get("contamination_flag", False))
    rows.append({
        "phylotype_id": r["phylotype_id"],
        "novelty_class": r.get("novelty_class", ""),
        "oral_clade": bool(r.get("oral_clade", False)),
        "candidate_novel": bool(r.get("candidate_novel", False)),
        "n_members": n,
        "n_flagged": n_flag,
        "n_clean": n - n_flag,
        "max_prevalence": r.get("max_prevalence", 0),
        "max_subject_prevalence": r.get("max_subject_prevalence", 0),
        "n_bioprojects_total": r.get("n_bioprojects_total", 0),
        "n_continents": r.get("n_continents", 0),
        "n_regions": r.get("regions", "").count(";") + 1 if pd.notna(r.get("regions")) and str(r.get("regions")).strip() else 0,
        "total_reads_sum": r.get("total_reads_sum", 0),
        "clean_phylotype": n_flag == 0,
    })

agg = pd.DataFrame(rows)
for c in ["max_prevalence", "max_subject_prevalence", "n_bioprojects_total", "n_continents", "n_regions", "total_reads_sum"]:
    agg[c] = pd.to_numeric(agg[c], errors="coerce").fillna(0).astype(int)

# support_score: subject-level prevalence (distinct people, not samples). 2026-08-13.
agg["support_score"] = (agg["max_subject_prevalence"] * agg["n_bioprojects_total"] * agg["n_regions"])

# ============================================================
# 3. Save full table
# ============================================================
out = PIPE_V2 / "06b_phylotypes_aggregated.tsv"
agg.to_csv(out, sep="\t", index=False)
log(f"  saved: {out} ({len(agg):,} rows)")

# ============================================================
# 4. Diagnostic — clean phylotype distribution (for threshold setting)
# ============================================================
log("\n=== Diagnostic: clean phylotypes ===")
clean = agg[agg["clean_phylotype"]]
log(f"  clean phylotypes: {len(clean):,} / {len(agg):,}")
log(f"  clean candidates (≥10 subjects, not homd): {int((clean['candidate_novel']).sum()):,}")

diag = clean[["phylotype_id", "candidate_novel", "n_members", "n_flagged",
              "max_subject_prevalence", "n_bioprojects_total", "n_continents", "n_regions",
              "total_reads_sum", "support_score"]]
diag.to_csv(PIPE_V2 / "06b_diagnostic_clean.tsv", sep="\t", index=False)

# support distribution among clean candidates
cc = clean[clean["candidate_novel"]]
log(f"\n  clean candidates by support (for choosing 'strongest' cutoffs):")
for col, cuts in [("max_subject_prevalence", [10, 20, 50, 100, 200, 500]),
                  ("n_bioprojects_total", [2, 3, 5, 10, 20]),
                  ("n_regions", [1, 2, 3, 4, 5])]:
    log(f"  {col}:")
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        n = int(((cc[col] >= lo) & (cc[col] < hi)).sum())
        log(f"    [{lo},{hi}): {n:,}")
    n = int((cc[col] >= cuts[-1]).sum())
    log(f"    ≥{cuts[-1]}: {n:,}")

# how many would pass if thresholds were e.g. subj≥50 AND bp≥5 AND regions≥2
n_strict = int(((cc["max_subject_prevalence"] >= 50) & (cc["n_bioprojects_total"] >= 5) & (cc["n_regions"] >= 2)).sum())
n_med = int(((cc["max_subject_prevalence"] >= 30) & (cc["n_bioprojects_total"] >= 3) & (cc["n_regions"] >= 2)).sum())
n_relax = int(((cc["max_subject_prevalence"] >= 10) & (cc["n_bioprojects_total"] >= 2) & (cc["n_regions"] >= 1)).sum())
log(f"\n  example cutoffs (clean candidates only):")
log(f"    strict  (subj≥50 & bp≥5 & reg≥2): {n_strict:,}")
log(f"    medium  (subj≥30 & bp≥3 & reg≥2): {n_med:,}")
log(f"    relaxed (subj≥10 & bp≥2 & reg≥1): {n_relax:,}")

# ============================================================
# 5. Summary
# ============================================================
summary = {
    "total_phylotypes": int(len(agg)),
    "clean_phylotypes": int(clean.shape[0]),
    "candidate_novel": int(agg["candidate_novel"].sum()),
    "clean_candidates": int(clean["candidate_novel"].sum()),
    "example_cutoffs": {
        "strict_subj50_bp5_reg2": n_strict,
        "medium_subj30_bp3_reg2": n_med,
        "relaxed_subj10_bp2_reg1": n_relax,
    },
    "flagged_distribution": {int(k): int(v) for k, v in agg["n_flagged"].value_counts().sort_index().items()},
}
with open(PIPE_V2 / "06b_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
log(f"\nSaved: {PIPE_V2 / '06b_summary.json'}")
log("\nDone!")

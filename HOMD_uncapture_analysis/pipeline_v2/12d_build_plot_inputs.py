#!/usr/bin/env python3
"""
pipeline_v2 Step 12d: build plot inputs for the two new result-layer figures.

  1. anchors_supported186.tsv   — one row per SUPPORTED-CATALOG phylotype (186):
     rep tip + metrics + label, for the CIRCULAR phylotype-circle tree.
  2. headline_profile53.tsv     — one row per HEADLINE phylotype (53):
     rep-level support + identity + regions/sites/continents, for the
     prevalence-profile figure (adapted from region_aware_pipeline/17_*).

Outputs (pipeline_v2/):
  - 12d_anchors_supported.tsv
  - 12d_headline_profile.tsv
"""
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"

log = lambda s: print(f"[12d_inputs] {s}", flush=True)


def parse_set(val):
    if pd.isna(val) or val == "":
        return []
    s = str(val).strip("{}")
    for sep in ["; ", ";", ", ", ","]:
        if sep in s:
            return [x.strip().strip("'\"") for x in s.split(sep)]
    return [s.strip().strip("'\"")]


def main():
    sup = pd.read_csv(PIPE_V2 / "11_supported_clean.tsv", sep="\t", low_memory=False)
    hl = pd.read_csv(PIPE_V2 / "11_headline_candidates.tsv", sep="\t", low_memory=False)
    log(f"supported: {len(sup):,} | headline: {len(hl):,}")

    mp = pd.read_csv(PIPE_V2 / "10_asv_phylotype_map.tsv", sep="\t", low_memory=False,
                     usecols=["derep_id", "phylotype_id"])
    s04 = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False,
                      usecols=["derep_id", "continents", "countries"])
    s01 = pd.read_csv(PIPE_V2 / "01_chimera_filtered.tsv", sep="\t", low_memory=False,
                      usecols=["derep_id", "sites"])

    def agg_set(pid, col_df, col):
        members = mp.loc[mp["phylotype_id"] == pid, "derep_id"]
        vals = col_df.loc[col_df["derep_id"].isin(members), col]
        out = set()
        for v in vals.dropna():
            out.update(parse_set(v))
        return ";".join(sorted(out))

    # ---- 1. anchors for the supported-catalog circles tree ----
    rows = []
    for _, r in sup.iterrows():
        label = r["rep_species"] if pd.notna(r["rep_species"]) and str(r["rep_species"]).strip() \
            else (r["rep_genus"] if pd.notna(r["rep_genus"]) and str(r["rep_genus"]).strip()
                  else r["phylotype_id"])
        rows.append({
            "anchor_tip": r["rep_derep_id"],
            "phylotype_id": r["phylotype_id"],
            "n_members": int(r["n_asvs"]),
            "max_subject_prevalence": r["max_subject_prevalence"],
            "n_bioprojects_total": r["n_bioprojects_total"],
            "n_continents": len(agg_set(r["phylotype_id"], s04, "continents").split(";")) if agg_set(r["phylotype_id"], s04, "continents") else 0,
            "n_regions": r["n_regions"],
            "oral_clade": bool(str(r["oral_clade"]).lower() in ("true", "1")),
            "tree_label": str(label),
        })
    anchors = pd.DataFrame(rows)
    anchors.to_csv(PIPE_V2 / "12d_anchors_supported.tsv", sep="\t", index=False)
    log(f"saved 12d_anchors_supported.tsv ({len(anchors):,} rows)")

    # ---- 2. headline prevalence profile ----
    prof_rows = []
    for _, r in hl.iterrows():
        regions = r["regions"] if pd.notna(r["regions"]) else ""
        prof_rows.append({
            "phylotype_id": r["phylotype_id"],
            "derep_id": r["rep_derep_id"],
            "n_asvs": int(r["n_asvs"]),
            "subject_prevalence": r["max_subject_prevalence"],
            "n_bioprojects": r["n_bioprojects_total"],
            "total_reads": r["total_reads_sum"],
            "regions": regions,
            "sites": agg_set(r["phylotype_id"], s01, "sites"),
            "continents": agg_set(r["phylotype_id"], s04, "continents"),
            "countries": agg_set(r["phylotype_id"], s04, "countries"),
            "gg2_identity": r["rep_identity"],
            "homd_identity": r["homd_best_id"],
            "gg2_species": r["rep_species"] if pd.notna(r["rep_species"]) else "",
            "gg2_genus": r["rep_genus"] if pd.notna(r["rep_genus"]) else "",
            "gg2_phylum": "",
        })
    prof = pd.DataFrame(prof_rows)
    prof.to_csv(PIPE_V2 / "12d_headline_profile.tsv", sep="\t", index=False)
    log(f"saved 12d_headline_profile.tsv ({len(prof):,} rows)")
    log("Done!")


if __name__ == "__main__":
    main()

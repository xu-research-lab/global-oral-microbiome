#!/usr/bin/env python3
"""
pipeline_v2 Step 07: Build the two candidate tables from the aggregated phylotypes.

Sources:
  - 06b_phylotypes_aggregated.tsv    : per-phylotype composition + markers
  - 06_phylotypes_novelty_classified.tsv : rep ASV, genus, novelty hits

Outputs:
  - 07_strict_candidates.tsv      : the strongest candidates (276-equivalent)
  - 07_supplementary_table.tsv    : all phylotypes for the paper supplement

strict_candidate cutoffs (set from the data distribution):
  clean_phylotype (no contamination-flagged member)
  & candidate_novel (>=10 subjects, not HOMD-captured)
  & max_subject_prevalence >= 50  & n_bioprojects_total >= 5  & n_regions >= 2
  (max_subject_prevalence = distinct (bioproject, subject_id) individuals)

Residence columns were removed on 2026-08-05 (oral-vs-blank test dropped).
"""
import json
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
log = lambda s: print(f"[07_tables] {s}", flush=True)

# ---- load ----
agg = pd.read_csv(PIPE_V2 / "06b_phylotypes_aggregated.tsv", sep="\t", low_memory=False)
nov = pd.read_csv(PIPE_V2 / "06_phylotypes_novelty_classified.tsv", sep="\t", low_memory=False)
log(f"  aggregated: {len(agg):,} | novelty-classified: {len(nov):,}")

# ---- merge ----
# GG2/SILVA vsearch columns dropped on 2026-08-13 (novelty = placement-clade vs HOMD).
nov_cols = ["phylotype_id", "rep_derep_id", "rep_gg2_genus", "regions",
            "homd_best_id", "homd_target"]
tbl = agg.merge(nov[nov_cols], on="phylotype_id", how="left")

# ---- strict_candidate flag (subject-level prevalence, 2026-08-13) ----
strict_mask = (
    tbl["clean_phylotype"]
    & tbl["candidate_novel"]
    & (tbl["max_subject_prevalence"] >= 50)
    & (tbl["n_bioprojects_total"] >= 5)
    & (tbl["n_regions"] >= 2)
)
tbl["strict_candidate"] = strict_mask
log(f"  strict candidates: {int(strict_mask.sum()):,}")

# ---- 07_strict_candidates.tsv ----
strict_cols = ["phylotype_id", "novelty_class", "candidate_novel", "n_members",
               "n_flagged", "n_clean",
               "max_subject_prevalence", "max_prevalence",
               "n_bioprojects_total", "n_continents", "n_regions",
               "total_reads_sum", "clean_phylotype", "support_score"]
strict = tbl.loc[strict_mask, strict_cols].sort_values("support_score", ascending=False)
strict.to_csv(PIPE_V2 / "07_strict_candidates.tsv", sep="\t", index=False)
log(f"  saved 07_strict_candidates.tsv ({len(strict):,} rows)")

# ---- 07_supplementary_table.tsv ----
supp_cols = ["phylotype_id", "strict_candidate", "novelty_class", "candidate_novel",
             "rep_derep_id", "rep_gg2_genus",
             "n_members", "n_clean", "n_flagged",
             "max_subject_prevalence", "max_prevalence",
             "n_bioprojects_total", "n_continents", "n_regions",
             "regions", "total_reads_sum", "clean_phylotype", "support_score",
             "homd_best_id", "homd_target"]
supp = tbl[supp_cols].sort_values(["strict_candidate", "support_score"], ascending=[False, False])
supp.to_csv(PIPE_V2 / "07_supplementary_table.tsv", sep="\t", index=False)
log(f"  saved 07_supplementary_table.tsv ({len(supp):,} rows)")

# ---- summary ----
summary = {
    "total_phylotypes": int(len(tbl)),
    "strict_candidates": int(strict_mask.sum()),
    "strict_cutoffs": {"clean_phylotype": True, "candidate_novel": True,
                       "max_subject_prevalence_min": 50, "n_bioprojects_min": 5, "n_regions_min": 2},
    "strict_by_novelty": {k: int(v) for k, v in tbl.loc[strict_mask, "novelty_class"].value_counts().items()},
}
with open(PIPE_V2 / "07_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
log(f"  saved 07_summary.json")
log("Done!")

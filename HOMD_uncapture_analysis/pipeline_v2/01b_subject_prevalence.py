#!/usr/bin/env python3
"""
pipeline_v2 Step 01b: Subject-level prevalence (distinct subjects, not samples).

Motivation (decision 2026-08-13):
  Sample-level prevalence over-counts when a subject contributes multiple
  (matched / time-series) samples — up to 983 samples for one subject in this
  dataset. Novel-oral-taxon calls should rest on how many DISTINCT PEOPLE carry
  a phylotype, not how many tubes.

  subject identity = (bioproject, subject_id). subject_id strings are NOT
  globally unique — the same label ("B1") is reused for different people across
  datasets, so the composite is required. Total distinct subjects = 52,298.

  For each support-floor derep sequence, subject_prevalence = the number of
  distinct (bioproject, subject_id) individuals in whose sample(s) the sequence
  is detected. Because the composite identity is project-scoped, the per-project
  distinct-subject counts are disjoint across projects and simply sum.

Inputs:
  - 04_support_floor_passed.tsv   : target derep_ids + sequences (54,291)
  - /beegfs/dongbiao/oral/*/results/DADA2.biom : per-project ASV × sample tables
  - data/metadata-v2-merge-oral.csv : sample_id → subject_id per bioproject

Outputs:
  - 01b_derep_subject_prevalence.tsv : derep_id, sequence, subject_prevalence
  - 01b_subject_stats.json           : total subjects, coverage diagnostics
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
ASV_ROOT = "/beegfs/dongbiao/oral"
METADATA = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/data/metadata-v2-merge-oral.csv"
log = lambda s: print(f"[01b_subj] {s}", flush=True)


def main(limit_projects=None):
    t0 = time.time()

    # ============================================================
    # 1. Target derep set (support floor) — derep_id + DNA sequence
    # ============================================================
    log("Loading 04_support_floor_passed.tsv (target dereps)...")
    sf = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False)
    sf["derep_id"] = sf["derep_id"].astype(str)
    sf["sequence"] = sf["sequence"].astype(str)
    seq_to_derep = dict(zip(sf["sequence"], sf["derep_id"]))
    target_seqs = set(seq_to_derep.keys())
    log(f"  target dereps: {len(target_seqs):,}")

    # ============================================================
    # 2. Metadata: sample_id → (bioproject, subject_id)
    # ============================================================
    log("Loading metadata (sample → subject)...")
    md = pd.read_csv(METADATA, low_memory=False)
    md["sample_id"] = md["sample_id"].astype(str)
    md["bioproject"] = md["bioproject"].astype(str)
    md["subject_id"] = md["subject_id"].astype(str).fillna("NA")
    sample_subject = dict(zip(md["sample_id"], md["bioproject"] + "|" + md["subject_id"]))
    n_subjects_total = md.groupby(["bioproject", "subject_id"]).ngroups
    log(f"  samples: {len(md):,} | distinct subjects (bioproject,subject_id): {n_subjects_total:,}")

    # ============================================================
    # 3. Per-project: count distinct subjects per target sequence
    # ============================================================
    from biom import load_table

    projects = sorted(md["bioproject"].unique())
    if limit_projects:
        projects = projects[:limit_projects]
        log(f"  TEST MODE: first {len(projects)} projects")

    subj_count = {}   # derep_id → subject_prevalence (running sum)
    n_missing_biom = 0
    n_sample_unmapped = 0
    n_sample_total = 0
    n_proj_target_rows = 0

    for k, bp in enumerate(projects, 1):
        biom_path = os.path.join(ASV_ROOT, bp, "results", "DADA2.biom")
        if not os.path.exists(biom_path):
            n_missing_biom += 1
            continue
        table = load_table(biom_path)
        obs_ids = list(map(str, table.ids(axis="observation")))
        sample_ids = list(map(str, table.ids(axis="sample")))
        mat = table.matrix_data.tocsr()
        indptr, indices = mat.indptr, mat.indices

        # sample index → composite subject (or None if not in metadata)
        n_sample_total += len(sample_ids)
        subj_of = [sample_subject.get(s) for s in sample_ids]
        n_unmapped = sum(1 for x in subj_of if x is None)
        n_sample_unmapped += n_unmapped

        # rows that are target sequences
        idx_obs = {oid: i for i, oid in enumerate(obs_ids) if oid in target_seqs}
        if not idx_obs:
            continue

        for oid, i in idx_obs.items():
            cols = indices[indptr[i]:indptr[i + 1]]
            if cols.size == 0:
                continue
            subs = {subj_of[c] for c in cols if subj_of[c] is not None}
            n_proj_target_rows += 1
            if subs:
                derep = seq_to_derep[oid]
                subj_count[derep] = subj_count.get(derep, 0) + len(subs)

        if k % 50 == 0:
            log(f"  {k}/{len(projects)} projects | covered {len(subj_count):,}/{len(target_seqs):,} dereps | elapsed {time.time()-t0:.0f}s")

    # ============================================================
    # 4. Output
    # ============================================================
    rows = [{"derep_id": d, "subject_prevalence": subj_count.get(d, 0)} for d in sf["derep_id"]]
    out = pd.DataFrame(rows)
    out["sequence"] = sf["sequence"].values
    out.to_csv(PIPE_V2 / "01b_derep_subject_prevalence.tsv", sep="\t", index=False)
    log(f"  saved 01b_derep_subject_prevalence.tsv ({len(out):,} rows)")

    n_detected = int((out["subject_prevalence"] > 0).sum())
    stats = {
        "total_distinct_subjects": n_subjects_total,
        "target_dereps": len(target_seqs),
        "dereps_detected_in_any_subject": n_detected,
        "subject_prevalence_range": [int(out["subject_prevalence"].min()),
                                     int(out["subject_prevalence"].max())],
        "projects_processed": len(projects),
        "projects_missing_biom": n_missing_biom,
        "samples_unmapped_to_subject": n_sample_unmapped,
        "samples_total": n_sample_total,
        "target_obs_rows_found_across_projects": n_proj_target_rows,
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(PIPE_V2 / "01b_subject_stats.json", "w") as f:
        json.dump(stats, f, indent=2)
    log(f"  saved 01b_subject_stats.json")
    log(f"  total subjects: {n_subjects_total:,} | dereps with ≥1 subject: {n_detected:,} | elapsed {time.time()-t0:.0f}s")

    # preview: max subject-prevalence dereps
    top = out.sort_values("subject_prevalence", ascending=False).head(10)
    log("  top subject-prevalence dereps:")
    for _, r in top.iterrows():
        log(f"    {r['derep_id']:>12} | {int(r['subject_prevalence'])}")
    log("Done!")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else None
    main(limit_projects=n)

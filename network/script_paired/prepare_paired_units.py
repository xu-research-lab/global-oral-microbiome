#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_paired_units.py — Superset of (project, status) units for the
paired HC/PD/DC/OSCC networks.

Unit selection (decided 2026-10-07):
  A project (bioproject_id) enters the pool if it has >= 1 sample of at
  least one disease group (PD / DC / OSCC) AND >= 1 Healthy Control sample
  (body site ignored — all-site). It is then split into one unit per status
  present, so a paired project contributes up to 4 units and every disease
  network has a matched HC unit available at the pooling layer.

The SparCC layer is a superset: inclusion rules (which units enter NetMoss)
are applied later in backbone_pooling.py, where they cost nothing to vary
because per-unit SparCC matrices are cached and reused.

Floors (unit-level sanity checks only, NOT project inclusion):
  MIN_UNIT_SAMPLES = 3   a unit needs >= 3 samples
  MIN_FEAT         = 5   genus kept if present in >= 5 samples of the unit
                         (units with < 5 samples therefore yield no genera
                         and are dropped, same as prepare_hc_pd_dc_oscc.py)

Outputs -> result_paired/
  paired_projects.csv        per project x status raw sample counts (audit)
  biom_paired/<grp>_genus/<study>.tsv   per-unit genus count matrices
  units.csv                  unit table for backbone_pooling.py
  study_list_genus.txt       "<grp>_genus <study>" lines for the slurm array

Usage (from network/):
  python script_hc_pd_dc_oscc/prepare_paired_units.py
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ──────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent                                  # network/
METADATA_TSV = ROOT / "data" / "metadata_after_cleaning.tsv"
COUNTS_TSV = ROOT / "data" / "genus_table_counts.tsv"
OUT_BASE = ROOT / "result_paired"

# ──────────────────────────────────────────────────────────────────────
# Criteria
# ──────────────────────────────────────────────────────────────────────
STUDY_COL = "bioproject_id"
MIN_UNIT_SAMPLES = 3      # a SparCC unit needs >= 3 samples
MIN_FEAT = 5              # genus kept if present in >= 5 samples of the unit

DISEASES = ["Periodontal Disease", "Dental Caries", "OSCC"]
GROUPS = {
    "HC": ["Healthy Control"],
    "PD": ["Periodontal Disease"],
    "DC": ["Dental Caries"],
    "OSCC": ["OSCC"],
}


def select_paired_projects(meta):
    """Projects with >=1 disease sample AND >=1 HC sample (all sites)."""
    counts = (meta.loc[meta["site_cleaned"].notna()]
              .groupby(STUDY_COL)["status_cleaned"]
              .value_counts().unstack(fill_value=0))
    for c in GROUPS.values():
        col = c[0]
        if col not in counts.columns:
            counts[col] = 0
    cols = ["Healthy Control"] + DISEASES
    counts = counts[cols]

    paired = counts[(counts[DISEASES].sum(axis=1) > 0) &
                    (counts["Healthy Control"] > 0)]
    paired.to_csv(OUT_BASE / "paired_projects.csv", sep="\t")
    return paired.index, counts


def write_unit_tsv(counts, grp, study_id, samples, out_dir):
    """Subset one unit's samples, prevalence-filter genera, write counts TSV."""
    cols = [s for s in samples if s in counts.columns]
    if len(cols) < MIN_UNIT_SAMPLES:
        print(f"  [{grp}/{study_id}] too few samples ({len(cols)}), skipped")
        return None
    sub = counts.loc[:, cols]
    present = (sub > 0).sum(axis=1)
    sub = sub.loc[(present >= MIN_FEAT).to_numpy()]
    if sub.shape[0] < 2:
        print(f"  [{grp}/{study_id}] too few genera ({sub.shape[0]}), skipped")
        return None
    if (sub.sum(axis=0) == 0).any():
        # all-zero samples break SparCC; drop them (and re-check)
        sub = sub.loc[:, sub.sum(axis=0) > 0]
        if sub.shape[0] < 2 or sub.shape[1] < MIN_UNIT_SAMPLES:
            print(f"  [{grp}/{study_id}] degenerate after zero-sample drop, "
                  f"skipped")
            return None
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{study_id}.tsv"
    sub.astype(np.int64).to_csv(out_path, sep="\t")
    return dict(n_feat=sub.shape[0], n_samp=sub.shape[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    args = ap.parse_args()

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(METADATA_TSV, sep="\t", dtype=str, low_memory=False)
    print(f"[meta] global metadata: {len(meta)} samples")

    projects, raw_counts = select_paired_projects(meta)
    print(f"[meta] paired projects (>=1 disease sample AND >=1 HC sample): "
          f"{len(projects)}")

    print("\n[genus] loading global genus count table ...")
    counts = pd.read_csv(COUNTS_TSV, sep="\t", index_col=0).astype(np.float32)
    print(f"[genus] genus table: {counts.shape[0]} genera x "
          f"{counts.shape[1]} samples")

    unit_rows, lines = [], []
    for study_id in sorted(projects):
        for grp, diseases in GROUPS.items():
            samples = meta.loc[(meta[STUDY_COL] == study_id) &
                               (meta["status_cleaned"].isin(diseases)) &
                               meta["site_cleaned"].notna(),
                               "sample_id"].astype(str).tolist()
            if not samples:
                continue
            n_raw = len(samples)
            out_dir = OUT_BASE / "biom_paired" / f"{grp}_genus"
            got = write_unit_tsv(counts, grp, study_id, samples, out_dir)
            if got is None:
                got = dict(n_feat=0, n_samp=0)
            unit_rows.append(dict(unit_id=f"{grp}_genus", category=grp,
                                  study=str(study_id), n_samples=len(samples),
                                  n_genera=got["n_feat"],
                                  kept=bool(got["n_feat"] >= 2),
                                  path=(str(out_dir / f"{study_id}.tsv")
                                        if got["n_feat"] >= 2 else "")))
            if got["n_feat"] >= 2:
                lines.append(f"{grp}_genus {study_id}")
                print(f"  [{grp}/{study_id}] {got['n_feat']} genera x "
                      f"{got['n_samp']}/{n_raw} samples")

    pd.DataFrame(unit_rows).to_csv(OUT_BASE / "units.csv", index=False)
    (OUT_BASE / "study_list_genus.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")

    kept = [r for r in unit_rows if r["kept"]]
    print(f"\n[units] {len(kept)}/{len(unit_rows)} (project,status) units "
          f"kept: "
          f"HC={sum(1 for r in kept if r['category']=='HC')} "
          f"PD={sum(1 for r in kept if r['category']=='PD')} "
          f"DC={sum(1 for r in kept if r['category']=='DC')} "
          f"OSCC={sum(1 for r in kept if r['category']=='OSCC')} "
          f"-> {OUT_BASE / 'units.csv'}")
    print(f"[units] {len(lines)} units -> study_list_genus.txt")
    print("[units] next: sbatch script_hc_pd_dc_oscc/run_slurm_paired.sh")


if __name__ == "__main__":
    main()
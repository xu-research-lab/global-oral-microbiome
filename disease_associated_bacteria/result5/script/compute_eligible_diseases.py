#!/usr/bin/env python3
"""
Compute eligible diseases for da_v4 DA analysis.

Rules (identical to the runners):
  - matched BP: the same bioproject has Case > 0 and HC > 0
  - across all matched BPs of that disease: Case >= 30 and HC >= 30
  - number of matched BPs >= 2

Basis: samples passing the shared cleaning filter total reads >= MIN_READS
(1000). This is the starting point of all three depth tracks, so a single
eligible_diseases.txt is used by raw / rare1k / rare5k.

eligible_diseases.csv additionally records the counts under each rarefaction
depth (1000 / 5000) so the effect of dropping shallow samples is documented.
"""
from __future__ import print_function
import os
from pathlib import Path
import numpy as np
import pandas as pd
import biom

# Repository root defaults to this script's folder; raw inputs (not
# distributed, see README) default to data/ under it.
WORK_DIR = Path(os.environ.get("DA_WORK_DIR", Path(__file__).resolve().parent.parent))
DATA_DIR = Path(os.environ.get("DA_DATA_DIR", WORK_DIR / "data"))
BIOM_PATH = DATA_DIR / "biom_oral_qc_dedup.biom"
META_PATH = DATA_DIR / "metadata_oral_qc_dedup.tsv"

MIN_READS = 1000          # shared sample filter, defines the eligible list
TRACK_DEPTHS = [1000, 5000]   # extra thresholds reported for reference

SITE_MAP = {
    "saliva": "Saliva", "tongue": "Tongue", "tooth": "Tooth",
    "oral swab": "Oral_swab", "oral mucosa": "Oral_Mucosa",
    "gingiva": "Gingiva", "other": "Others",
}
STATUS_MAP = {"HC": "Healthy Control", "Healthy": "Healthy Control"}
HC = "Healthy Control"


def load_metadata():
    meta = pd.read_csv(META_PATH, sep="\t", low_memory=False)
    str_cols = [
        "sample_id", "bioproject_id", "site_cleaned", "continent_group",
        "status_cleaned", "sex_cleaned", "age_cleaned_group",
        "sequencing_region_cleaned",
    ]
    for col in str_cols:
        if col in meta.columns:
            meta[col] = meta[col].astype(str).str.strip()
            meta.loc[meta[col].isin(["", "nan", "None", "NA"]), col] = pd.NA
    if "site_cleaned" in meta.columns:
        meta["site_cleaned"] = meta["site_cleaned"].replace(SITE_MAP)
    if "status_cleaned" in meta.columns:
        meta["status_cleaned"] = meta["status_cleaned"].replace(STATUS_MAP)
    if "age_cleaned_group" in meta.columns:
        meta["age_cleaned_group"] = meta["age_cleaned_group"].str.replace(
            "＞", ">", regex=False
        )
    return meta


def evaluate_diseases(meta, sample_set, track_name, verbose=True):
    """Stats for every disease label within sample_set (pass flag included)."""
    sub_meta = meta[meta["sample_id"].isin(sample_set)].copy()
    valid = (
        sub_meta["bioproject_id"].notna()
        & sub_meta["status_cleaned"].notna()
        & (sub_meta["status_cleaned"] != "")
        & (sub_meta["status_cleaned"] != "Other")
    )
    sub_meta = sub_meta[valid]
    diseases = sorted(d for d in sub_meta["status_cleaned"].unique() if d != HC)
    rows = []
    for disease in diseases:
        disease_bp = sub_meta.loc[
            sub_meta["status_cleaned"] == disease, "bioproject_id"
        ].unique()
        cohort = sub_meta[
            sub_meta["status_cleaned"].isin([HC, disease])
            & sub_meta["bioproject_id"].isin(disease_bp)
        ]
        if cohort.empty:
            continue
        bp_tbl = pd.crosstab(cohort["bioproject_id"], cohort["status_cleaned"])
        if disease not in bp_tbl.columns or HC not in bp_tbl.columns:
            continue
        valid_bp = bp_tbl.index[(bp_tbl[disease] > 0) & (bp_tbl[HC] > 0)]
        cohort = cohort[cohort["bioproject_id"].isin(valid_bp)]
        n_case = int((cohort["status_cleaned"] == disease).sum())
        n_hc = int((cohort["status_cleaned"] == HC).sum())
        n_bp = int(cohort["bioproject_id"].nunique())
        passes = n_case >= 30 and n_hc >= 30 and n_bp >= 2
        per_bp = [
            "%s(%d/%d)" % (bp, int(bp_tbl.loc[bp, disease]), int(bp_tbl.loc[bp, HC]))
            for bp in sorted(valid_bp.astype(str))
        ]
        rows.append({
            "Disease": disease,
            "track": track_name,
            "n_case": n_case,
            "n_hc": n_hc,
            "n_bp": n_bp,
            "passes_filter": passes,
            "matched_bps": "; ".join(per_bp),
        })
        if passes and verbose:
            print(
                "  [%s] %-34s case=%4d HC=%4d BP=%2d"
                % (track_name, disease, n_case, n_hc, n_bp)
            )
    rows.sort(key=lambda r: (-r["n_case"], r["Disease"]))
    return rows


def main():
    WORK_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading metadata:", META_PATH)
    meta = load_metadata()
    meta = meta.dropna(subset=["sample_id"]).copy()
    meta["sample_id"] = meta["sample_id"].astype(str)

    print("Loading BIOM for sample ids / library sizes:", BIOM_PATH)
    table = biom.load_table(str(BIOM_PATH))
    sample_ids = np.array(table.ids(axis="sample"), dtype=str)
    totals = np.asarray(table.sum(axis="sample")).ravel()
    print("BIOM samples: %d" % len(sample_ids))

    sets = {}
    sets["min%d" % MIN_READS] = set(sample_ids[totals >= MIN_READS])
    for depth in TRACK_DEPTHS:
        sets["rare%dk" % (depth // 1000)] = set(sample_ids[totals >= depth])
    for name, sset in sets.items():
        print("  %-8s samples: %d" % (name, len(sset)))

    base_key = "min%d" % MIN_READS
    print("\n=== Eligible diseases (basis: total reads >= %d) ===" % MIN_READS)
    rows_base = evaluate_diseases(meta, sets[base_key], base_key)
    diseases_base = [r["Disease"] for r in rows_base if r["passes_filter"]]

    stats_by_track = {base_key: {r["Disease"]: r for r in rows_base}}
    for name in sets:
        if name == base_key:
            continue
        print("\n=== Reference: %s ===" % name)
        rows = evaluate_diseases(meta, sets[name], name)
        stats_by_track[name] = {r["Disease"]: r for r in rows}

    list_path = WORK_DIR / "eligible_diseases.txt"
    with list_path.open("w") as fh:
        for d in diseases_base:
            fh.write(d + "\n")
    print("\nWrote %d diseases -> %s" % (len(diseases_base), list_path))

    track_order = [base_key] + [k for k in sets if k != base_key]
    all_names = sorted(
        set().union(*[set(s) for s in stats_by_track.values()]),
        key=lambda d: (
            -stats_by_track[base_key].get(d, {}).get("n_case", 0),
            d,
        ),
    )
    combined = []
    for d in all_names:
        row = {"Disease": d, "pass_%s" % base_key: d in diseases_base}
        for name in track_order:
            s = stats_by_track[name].get(d, {})
            row["n_case_%s" % name] = s.get("n_case", "")
            row["n_hc_%s" % name] = s.get("n_hc", "")
            row["n_bp_%s" % name] = s.get("n_bp", "")
            row["passes_%s" % name] = s.get("passes_filter", False)
        row["matched_bps"] = stats_by_track[base_key].get(d, {}).get("matched_bps", "")
        combined.append(row)
    csv_path = WORK_DIR / "eligible_diseases.csv"
    pd.DataFrame(combined).to_csv(csv_path, index=False)
    print("Wrote:", csv_path)

    for name in track_order:
        n_pass = sum(
            1 for s in stats_by_track[name].values() if s.get("passes_filter")
        )
        print("Summary: %-8s passing diseases = %d" % (name, n_pass))


if __name__ == "__main__":
    main()

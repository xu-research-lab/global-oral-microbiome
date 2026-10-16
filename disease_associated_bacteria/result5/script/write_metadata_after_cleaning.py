#!/usr/bin/env python3
"""
Write post–data-cleaning metadata for the da_v4 raw track.

Keeps only samples present in genus_table_relabund.tsv (reads >= 1000 after
BIOM filtering in prepare_genus_tables.py) and applies the same light
field-cleaning used by the DA runners.
"""
from __future__ import print_function
import os
from pathlib import Path
import pandas as pd

# Repository root defaults to this script's folder; raw metadata (not
# distributed, see README) defaults to data/ under it.
WORK = Path(os.environ.get("DA_WORK_DIR", Path(__file__).resolve().parent.parent))
META_SRC = Path(os.environ.get("DA_DATA_DIR", WORK / "data")) / "metadata_oral_qc_dedup.tsv"
RELAB = WORK / "genus_table_relabund.tsv"
OUT = WORK / "metadata_after_cleaning.tsv"

SITE_MAP = {
    "saliva": "Saliva",
    "tongue": "Tongue",
    "tooth": "Tooth",
    "oral swab": "Oral_swab",
    "oral mucosa": "Oral_Mucosa",
    "gingiva": "Gingiva",
    "other": "Others",
}
STATUS_MAP = {"HC": "Healthy Control", "Healthy": "Healthy Control"}
STR_COLS = [
    "sample_id",
    "bioproject_id",
    "site_cleaned",
    "continent_group",
    "status_cleaned",
    "sex_cleaned",
    "age_cleaned_group",
    "sequencing_region_cleaned",
]


def main():
    # Only the header is needed (sample IDs); do not load the abundance matrix.
    with RELAB.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
    sample_ids = header[1:]
    print("genus table samples:", len(sample_ids))

    meta = pd.read_csv(META_SRC, sep="\t", low_memory=False)
    print("source metadata rows:", len(meta))

    for col in STR_COLS:
        if col in meta.columns:
            meta[col] = meta[col].astype(str).str.strip()
            meta.loc[meta[col].isin(["", "nan", "None", "NA"]), col] = pd.NA

    if "site_cleaned" in meta.columns:
        meta["site_cleaned"] = meta["site_cleaned"].replace(SITE_MAP)
    if "status_cleaned" in meta.columns:
        meta["status_cleaned"] = meta["status_cleaned"].replace(STATUS_MAP)
    if "age_cleaned_group" in meta.columns:
        meta["age_cleaned_group"] = (
            meta["age_cleaned_group"]
            .astype(str)
            .str.replace("\uff1e", ">", regex=False)  # full-width ＞
        )
        meta.loc[
            meta["age_cleaned_group"].isin(["nan", "None", "NA", ""]),
            "age_cleaned_group",
        ] = pd.NA

    meta["sample_id"] = meta["sample_id"].astype(str).str.strip()
    meta = meta[
        meta["sample_id"].notna()
        & (meta["sample_id"] != "")
        & (meta["sample_id"] != "nan")
    ]
    meta = meta.drop_duplicates(subset=["sample_id"], keep="first")

    kept = set(sample_ids)
    sub = meta[meta["sample_id"].isin(kept)].copy()
    order = {sid: i for i, sid in enumerate(sample_ids)}
    sub["_ord"] = sub["sample_id"].map(order)
    sub = sub.sort_values("_ord").drop(columns=["_ord"])

    missing = [s for s in sample_ids if s not in set(sub["sample_id"])]
    print("matched metadata rows:", len(sub))
    print("genus samples missing in metadata:", len(missing))
    if missing:
        print("  e.g.", missing[:5])

    sub.insert(0, "da_v4_track", "raw")
    sub.insert(1, "da_v4_min_reads", 1000)
    sub.insert(2, "da_v4_in_genus_table", True)

    sub.to_csv(OUT, sep="\t", index=False)
    print("wrote:", OUT)
    print("columns:", len(sub.columns), "rows:", len(sub))
    print("status_cleaned top:")
    print(sub["status_cleaned"].value_counts().head(12).to_string())


if __name__ == "__main__":
    main()

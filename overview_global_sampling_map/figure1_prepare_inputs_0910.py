#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# figure1_prepare_inputs_0910.py
#
# Prepares the Figure 1 plotting input files from the raw metadata and the
# oral-gut paired-sample table. The R plotting script
# (figure1_render_panels_0910.R) reads these pre-computed files and only draws
# the figure; it does not re-derive any counts.
#
# Outputs (written to <script_dir>/input/):
#   cohort_site.tsv     label, n   (7 oral site categories, canonical order)
#   cohort_age.tsv       label, n   (6 age groups, canonical order)
#   cohort_sex.tsv       label, n   (F / M)
#   cohort_country.tsv   label, n   (top 10 countries by sample count)
#   cohort_region.tsv    label, n   (top 8 16S variable regions)
#   cohort_pairs.tsv     n_pairs, n_pair_ids  (oral-gut paired samples)
#   top20_diseases.tsv   status, n, category  (top 20 disease statuses)
#
# Inputs (paths unchanged):
#   <repo>/03-patch-pipeline-map/input/metadata_oral_qc.tsv   (80,788 samples)
#   <script_dir>/metadata_oral_gut_pair_v2_merge_qc.csv       (oral-gut pairs)
# ----------------------------------------------------------------------------

import csv
import os
import sys
from collections import Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Raw metadata lives in the 03-patch-pipeline-map/input/ folder (two levels up
# from this script's directory).
METADATA = os.path.join(SCRIPT_DIR, "..", "..", "input", "metadata_oral_qc.tsv")
PAIR_CSV = os.path.join(SCRIPT_DIR, "metadata_oral_gut_pair_v2_merge_qc.csv")
OUT_DIR  = os.path.join(SCRIPT_DIR, "input")

SITE_ORDER = ["Saliva", "Tooth", "Oral Mucosa", "Oral Swab",
              "Tongue", "Gingiva", "Others"]
AGE_ORDER  = ["0-6 months", "6-36 months", "3-12 years", "12-18 years",
              "18-60 years", ">60 years"]

# Disease categories for the top-20 panel (panel E).
ORAL_DISEASES = {"Dental Caries", "Periodontal Disease", "OSCC",
                 "Oral Mucositis", "Aphthous"}
CANCER_DISEASES = {"Lung Cancer", "Head and Neck Cancer", "Colorectal Cancer",
                   "Breast Cancer", "Acute Myeloid Leukemia", "Pancreatic Cancer"}


def is_na(v):
    return v is None or str(v).strip() == "" or str(v).strip().upper() == "NA"


def read_metadata(path):
    """Return a list of dicts for the 80,788-sample metadata table."""
    rows = []
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def count_column(rows, col, drop_na=True):
    """Count non-NA values of a metadata column, returning {value: n}."""
    counts = Counter()
    for r in rows:
        v = r.get(col, "")
        if drop_na and is_na(v):
            continue
        counts[v] += 1
    return counts


def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print("wrote", path, f"({len(rows)} rows)")


def main():
    if not os.path.exists(METADATA):
        sys.exit(f"metadata not found: {METADATA}")
    print("reading", METADATA)
    md = read_metadata(METADATA)
    print("n samples:", len(md))

    # ---- cohort composition: body site (7 oral categories) ----
    site_counts = count_column(md, "site_cleaned")
    site_rows = []
    for label in SITE_ORDER:
        # raw labels are lowercase; the "Others" category is stored as "Other"
        raw = "Other" if label == "Others" else label.lower()
        n = site_counts.get(raw, 0)
        site_rows.append((label, n))
    write_tsv(os.path.join(OUT_DIR, "cohort_site.tsv"),
              ["label", "n"], site_rows)

    # ---- cohort composition: age group (6 groups) ----
    age_counts = count_column(md, "age_cleaned_group")
    # normalize full-width "＞" to ">"
    age_counts = Counter({k.replace("＞", ">"): v for k, v in age_counts.items()})
    age_rows = [(label, age_counts.get(label, 0)) for label in AGE_ORDER]
    write_tsv(os.path.join(OUT_DIR, "cohort_age.tsv"),
              ["label", "n"], age_rows)

    # ---- cohort composition: sex ----
    sex_counts = count_column(md, "sex_cleaned")
    sex_rows = sorted(sex_counts.items(), key=lambda kv: -kv[1])
    write_tsv(os.path.join(OUT_DIR, "cohort_sex.tsv"),
              ["label", "n"], sex_rows)

    # ---- cohort composition: country (top 10) ----
    country_counts = count_column(md, "country_cleaned")
    country_rows = sorted(country_counts.items(), key=lambda kv: -kv[1])[:10]
    write_tsv(os.path.join(OUT_DIR, "cohort_country.tsv"),
              ["label", "n"], country_rows)

    # ---- cohort composition: 16S variable region (top 8) ----
    region_counts = count_column(md, "sequencing_region_cleaned")
    region_rows = sorted(region_counts.items(), key=lambda kv: -kv[1])[:8]
    write_tsv(os.path.join(OUT_DIR, "cohort_region.tsv"),
              ["label", "n"], region_rows)

    # ---- oral-gut paired samples ----
    if not os.path.exists(PAIR_CSV):
        print("WARNING: pair csv not found, writing NA:", PAIR_CSV)
        write_tsv(os.path.join(OUT_DIR, "cohort_pairs.tsv"),
                  ["n_pairs", "n_pair_ids"], [["NA", "NA"]])
    else:
        n_pairs = 0
        pair_ids = set()
        with open(PAIR_CSV, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                n_pairs += 1
                pair_ids.add(row.get("pair_id", ""))
        write_tsv(os.path.join(OUT_DIR, "cohort_pairs.tsv"),
                  ["n_pairs", "n_pair_ids"], [[n_pairs, len(pair_ids)]])
        print("  oral-gut paired samples:", n_pairs, "pairs:", len(pair_ids))

    # ---- top-20 disease statuses (healthy controls excluded) ----
    status_counts = count_column(md, "status_cleaned")
    status_counts.pop("HC", None)
    top20 = sorted(status_counts.items(), key=lambda kv: -kv[1])[:20]
    disease_rows = []
    for status, n in top20:
        if status in ORAL_DISEASES:
            cat = "Oral diseases"
        elif status in CANCER_DISEASES:
            cat = "Cancers"
        else:
            cat = "Systemic & other"
        disease_rows.append((status, n, cat))
    write_tsv(os.path.join(OUT_DIR, "top20_diseases.tsv"),
              ["status", "n", "category"], disease_rows)

    print("done. inputs written to", OUT_DIR)


if __name__ == "__main__":
    main()

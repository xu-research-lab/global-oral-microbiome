#!/usr/bin/env python3
"""Prepare HOMD tree tip annotation for ggtree rings."""

from __future__ import annotations

import argparse
import json
import math
import os
from collections import defaultdict

import pandas as pd

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_TAX = "/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03.taxonomy"


def log(msg):
    print(f"[tree-annot] {msg}", flush=True)


def load_homd_tax(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            target, tax = parts[0], parts[1]
            levels = tax.split(";")
            while len(levels) < 7:
                levels.append("NA")
            rows.append({
                "homd_target": target,
                "domain": levels[0] or "NA",
                "phylum": levels[1] or "NA",
                "class": levels[2] or "NA",
                "order": levels[3] or "NA",
                "family": levels[4] or "NA",
                "genus": levels[5] or "NA",
                "species_epithet": levels[6] or "NA",
                "taxonomy": tax,
            })
    return pd.DataFrame(rows)


def load_hits(path):
    cols = [
        "derep_id", "homd_target", "identity", "alnlen", "mismatches", "gaps", "qlo", "qhi", "tlo", "thi",
        "query_length", "target_length", "query_cov", "target_cov",
    ]
    hits = pd.read_csv(path, sep="\t", names=cols)
    hits = hits.sort_values(["derep_id", "identity", "query_cov", "alnlen", "homd_target"], ascending=[True, False, False, False, True])
    return hits.drop_duplicates("derep_id", keep="first")


def split_values(value):
    if pd.isna(value):
        return []
    return [v for v in str(value).split(";") if v and v != "nan"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--taxonomy", default=DEFAULT_TAX)
    args = parser.parse_args()

    burden_path = os.path.join(args.outdir, "homd_reference_hit_burden_qcov99.tsv")
    derep_path = os.path.join(args.outdir, "exact_derep_summary.tsv")
    hit_path = os.path.join(args.outdir, "homd_vsearch_best_hits.tsv")
    out_path = os.path.join(args.outdir, "homd_tree_tip_annotation.tsv")
    stats_path = os.path.join(args.outdir, "homd_tree_tip_annotation_stats.json")

    log("loading taxonomy and burden")
    tax = load_homd_tax(args.taxonomy)
    burden = pd.read_csv(burden_path, sep="\t")

    log("loading best hits and derep metadata for unique cohort counts")
    derep_cols = ["derep_id", "bioprojects", "countries", "continents"]
    derep = pd.read_csv(derep_path, sep="\t", usecols=derep_cols)
    hits = load_hits(hit_path)
    merged = hits.merge(derep, on="derep_id", how="left")
    merged["identity"] = pd.to_numeric(merged["identity"], errors="coerce").fillna(0)
    merged["query_cov"] = pd.to_numeric(merged["query_cov"], errors="coerce").fillna(0)

    cohort = defaultdict(lambda: {"bioprojects_any": set(), "countries_any": set(), "continents_any": set(),
                                  "bioprojects_id99": set(), "countries_id99": set(), "continents_id99": set()})
    for row in merged.itertuples(index=False):
        rec = cohort[row.homd_target]
        bps = split_values(row.bioprojects)
        countries = split_values(row.countries)
        continents = split_values(row.continents)
        rec["bioprojects_any"].update(bps)
        rec["countries_any"].update(countries)
        rec["continents_any"].update(continents)
        if row.identity >= 99 and row.query_cov >= 99:
            rec["bioprojects_id99"].update(bps)
            rec["countries_id99"].update(countries)
            rec["continents_id99"].update(continents)

    cohort_rows = []
    for target, rec in cohort.items():
        cohort_rows.append({
            "homd_target": target,
            "n_bioprojects_any": len(rec["bioprojects_any"]),
            "n_countries_any": len(rec["countries_any"]),
            "n_continents_any": len(rec["continents_any"]),
            "n_bioprojects_id99_qcov99": len(rec["bioprojects_id99"]),
            "n_countries_id99_qcov99": len(rec["countries_id99"]),
            "n_continents_id99_qcov99": len(rec["continents_id99"]),
        })
    cohort_df = pd.DataFrame(cohort_rows)

    annot = burden.merge(tax, on="homd_target", how="left").merge(cohort_df, on="homd_target", how="left")
    for col in ["phylum", "class", "order", "family", "genus", "species_epithet", "taxonomy"]:
        annot[col] = annot[col].fillna("NA")
    for col in ["n_bioprojects_any", "n_countries_any", "n_continents_any", "n_bioprojects_id99_qcov99", "n_countries_id99_qcov99", "n_continents_id99_qcov99"]:
        annot[col] = annot[col].fillna(0).astype(int)

    annot["log10_any_hit_derep"] = annot["any_hit_derep_sequences"].fillna(0).map(lambda x: math.log10(x + 1))
    annot["log10_id99_derep"] = annot["id99p0_qcov99_derep_sequences"].fillna(0).map(lambda x: math.log10(x + 1))
    annot["median_identity_for_plot"] = annot["any_hit_median_identity"].fillna(0)
    annot["identity_gap_from_100"] = (100 - annot["median_identity_for_plot"]).clip(lower=0)
    annot["fraction_asv_id99_qcov99"] = 0.0
    has_any_asv = annot["any_hit_source_asvs"].fillna(0) > 0
    annot.loc[has_any_asv, "fraction_asv_id99_qcov99"] = (
        annot.loc[has_any_asv, "id99p0_qcov99_source_asvs"].fillna(0)
        / annot.loc[has_any_asv, "any_hit_source_asvs"].fillna(0)
    )
    annot["fraction_asv_id99_qcov99"] = annot["fraction_asv_id99_qcov99"].clip(lower=0, upper=1)
    annot["clipped_identity_gap_to_99"] = (99 - annot["median_identity_for_plot"]).clip(lower=0, upper=10)
    annot["common_low_identity_score"] = annot["log10_any_hit_derep"] * annot["n_bioprojects_any"].map(lambda x: math.log10(x + 1)) * annot["identity_gap_from_100"]
    annot["reference_gap_score"] = (
        annot["log10_any_hit_derep"]
        * annot["n_bioprojects_any"].map(lambda x: math.log10(x + 1))
        * (1 - annot["fraction_asv_id99_qcov99"])
        * annot["clipped_identity_gap_to_99"]
    )

    # Stable phylum labels for plotting.
    annot["phylum_plot"] = annot["phylum"].replace({"": "NA"}).fillna("NA")
    top_phyla = set(annot["phylum_plot"].value_counts().head(12).index)
    annot.loc[~annot["phylum_plot"].isin(top_phyla), "phylum_plot"] = "Other"

    annot.to_csv(out_path, sep="\t", index=False)
    stats = {
        "tips": int(len(annot)),
        "with_any_hit": int((annot["any_hit_derep_sequences"] > 0).sum()),
        "with_id99_qcov99_hit": int((annot["id99p0_qcov99_derep_sequences"] > 0).sum()),
        "top_reference_gap_score": annot.sort_values("reference_gap_score", ascending=False)[
            ["homd_target", "homd_hmt", "homd_species", "phylum", "any_hit_derep_sequences", "any_hit_source_asvs", "id99p0_qcov99_source_asvs", "fraction_asv_id99_qcov99", "any_hit_median_identity", "n_bioprojects_any", "reference_gap_score"]
        ].head(20).to_dict(orient="records"),
        "output": out_path,
    }
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)
    log(f"wrote {out_path}")


if __name__ == "__main__":
    main()

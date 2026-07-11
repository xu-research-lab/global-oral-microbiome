#!/usr/bin/env python3
"""Per-taxon geographic reference gap: which bacterial taxa are missing from the reference
in non-Western oral samples.

This turns the sample-level "non-Western regions map worse to Greengenes2" observation
(see 23_cr_mapping_by_continent_amplicon_adjusted.py) into a per-taxon, amplicon-aware,
contamination-guarded statement: it identifies the specific *bacterial* 16S lineages that
recur in oral samples yet are poorly represented in HOMD/GG2, and asks whether that
reference gap is concentrated in non-Western continents.

Why this is the defensible version of the claim
------------------------------------------------
A low mapping rate is only interesting if the unmapped reads are genuine bacterial 16S
rather than host/plastid/chimera noise. This script therefore consumes the GG2 *bacterial-
only* low-identity ASV table (already screened to d__Bacteria and past the suspicious-
backbone / chloroplast-mito filter by scripts 14-15), so every ASV counted is plausible
novel bacterial diversity. It also reports results per amplicon region, because identity to
the reference is amplicon-dependent.

Inputs (produced by scripts 14-15 of this pipeline; live on the cluster)
------------------------------------------------------------------------
  --bacterial-table  low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_only.tsv
      Expected columns (subset used here):
        derep_id, total_reads, prevalence_samples, subject_prevalence,
        continents            (';'-joined list of continents where the ASV was observed),
        sequencing_regions    (';'-joined list of amplicon regions),
        gg2_best_identity, gg2_phylum, gg2_family, gg2_genus
  --metadata         metadata-v2-merge-oral.csv  (or metadata_all_for_qc_raw.tsv)
      Used only to count sampling effort (unique subjects per continent) so that
      reference-gap burden is normalized per 1000 subjects and not driven by sample count.

Outputs (under --outdir)
------------------------
  reference_gap_by_continent.tsv        burden per continent (n gap ASVs, per-1000-subjects)
  reference_gap_top_genera.tsv          top GG2 genera driving each continent's gap
  reference_gap_nonwestern_specific.tsv gap ASVs seen ONLY in non-Western continents + taxonomy
  reference_gap_by_continent.png/.pdf   figure (burden barplot + genus x continent heatmap)

NOTE ON STATUS: this script targets the cluster intermediates under
`/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline` and has
NOT been executed here (those tables are not in the local checkout). Column-name guards are
included; adjust `--*-col` flags if the upstream schema differs.
"""

from __future__ import annotations

import argparse
import os
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_BACTERIAL = os.path.join(
    DEFAULT_OUTDIR, "low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_only.tsv"
)
DEFAULT_META = "/beegfs/huangxiaochang/global_oral/metadata_v2_merge/metadata-v2-merge-oral.csv"

CONTINENTS = ["North America", "Europe", "Oceania", "Asia", "Africa", "South America"]
WESTERN = {"North America", "Europe", "Oceania"}
NON_WESTERN = {"Asia", "Africa", "South America"}


def log(msg: str) -> None:
    print(f"[geo-refgap] {msg}", flush=True)


def split_values(x):
    if pd.isna(x):
        return []
    return [v.strip() for v in str(x).split(";") if v.strip() and v.strip() not in {"nan", "NA"}]


def subjects_per_continent(meta_path: str) -> pd.Series:
    sep = "\t" if meta_path.endswith((".tsv", ".txt")) else ","
    md = pd.read_csv(meta_path, sep=sep, low_memory=False)
    subj = "subject_id" if "subject_id" in md.columns else "sample_id"
    cont = "continent_group" if "continent_group" in md.columns else "continent"
    n = md.dropna(subset=[cont]).groupby(cont)[subj].nunique()
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bacterial-table", default=DEFAULT_BACTERIAL)
    ap.add_argument("--metadata", default=DEFAULT_META)
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
    ap.add_argument("--continents-col", default="continents")
    ap.add_argument("--regions-col", default="sequencing_regions")
    ap.add_argument("--genus-col", default="gg2_genus")
    ap.add_argument("--weight-col", default="subject_prevalence",
                    help="per-ASV weight for burden (subject_prevalence | total_reads | prevalence_samples)")
    ap.add_argument("--amplicon", default=None,
                    help="restrict to ASVs observed in this amplicon region (e.g. V4) for an amplicon-controlled view")
    ap.add_argument("--top-genera", type=int, default=20)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.bacterial_table, sep="\t", low_memory=False)
    log(f"loaded {len(df):,} bacterial low-identity ASVs")
    for col in (args.continents_col, args.genus_col):
        if col not in df.columns:
            raise SystemExit(f"missing expected column '{col}'; available: {list(df.columns)}")

    df["_conts"] = df[args.continents_col].map(split_values)
    if args.amplicon:
        if args.regions_col not in df.columns:
            raise SystemExit(f"--amplicon set but region column '{args.regions_col}' absent")
        keep = df[args.regions_col].map(lambda x: args.amplicon in split_values(x))
        df = df[keep].copy()
        log(f"amplicon={args.amplicon}: {len(df):,} ASVs retained")

    w = pd.to_numeric(df.get(args.weight_col, 1), errors="coerce").fillna(0.0).to_numpy()
    n_subj = subjects_per_continent(args.metadata)

    # ---- burden per continent -------------------------------------------------
    rows = []
    for c in CONTINENTS:
        mask = df["_conts"].map(lambda s, c=c: c in s).to_numpy(bool)
        n_asv = int(mask.sum())
        wsum = float(w[mask].sum())
        subs = float(n_subj.get(c, np.nan))
        rows.append(dict(
            continent=c, group="non-Western" if c in NON_WESTERN else "Western",
            n_gap_asvs=n_asv, weighted_burden=round(wsum, 1),
            n_subjects=int(subs) if subs == subs else None,
            gap_asvs_per_1000_subjects=round(n_asv / subs * 1000, 2) if subs == subs and subs else None,
        ))
    burden = pd.DataFrame(rows)
    burden.to_csv(os.path.join(args.outdir, "reference_gap_by_continent.tsv"), sep="\t", index=False)
    log("wrote reference_gap_by_continent.tsv")
    print(burden.to_string(index=False))

    # ---- top GG2 genera per continent ----------------------------------------
    genus = df[args.genus_col].fillna("Unassigned").replace("", "Unassigned")
    per_cont_genus = {}
    for c in CONTINENTS:
        mask = df["_conts"].map(lambda s, c=c: c in s).to_numpy(bool)
        per_cont_genus[c] = Counter()
        for g, wi in zip(genus[mask], w[mask]):
            per_cont_genus[c][g] += wi
    top_genera = sorted({g for c in CONTINENTS for g, _ in per_cont_genus[c].most_common(args.top_genera)})
    mat = pd.DataFrame({c: {g: per_cont_genus[c].get(g, 0.0) for g in top_genera} for c in CONTINENTS})
    mat.index.name = "gg2_genus"
    mat.to_csv(os.path.join(args.outdir, "reference_gap_top_genera.tsv"), sep="\t")
    log("wrote reference_gap_top_genera.tsv")

    # ---- ASVs whose entire continent footprint is non-Western -----------------
    nw = df[df["_conts"].map(lambda s: len(s) > 0 and set(s).issubset(NON_WESTERN))].copy()
    keep_cols = [c for c in ["derep_id", "gg2_phylum", "gg2_family", args.genus_col,
                             "gg2_best_identity", args.weight_col, args.continents_col] if c in nw.columns]
    nw[keep_cols].sort_values(args.weight_col, ascending=False).to_csv(
        os.path.join(args.outdir, "reference_gap_nonwestern_specific.tsv"), sep="\t", index=False)
    log(f"wrote reference_gap_nonwestern_specific.tsv ({len(nw):,} non-Western-only gap ASVs)")

    # ---- figure ---------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.4]})
    bcol = ["#C44E52" if g == "non-Western" else "#4C72B0" for g in burden["group"]]
    yvals = burden["gap_asvs_per_1000_subjects"].fillna(0)
    axes[0].bar(burden["continent"], yvals, color=bcol)
    axes[0].set_ylabel("bacterial reference-gap ASVs per 1000 subjects")
    axes[0].set_title("A  Reference-gap burden by continent")
    axes[0].tick_params(axis="x", rotation=30)
    for lbl in axes[0].get_xticklabels():
        lbl.set_ha("right")

    heat = mat.reindex(columns=CONTINENTS)
    heat_norm = heat / heat.max(axis=1).replace(0, np.nan).to_numpy()[:, None]
    im = axes[1].imshow(heat_norm.values, aspect="auto", cmap="magma_r")
    axes[1].set_xticks(range(len(CONTINENTS)))
    axes[1].set_xticklabels(CONTINENTS, rotation=30, ha="right", fontsize=8)
    axes[1].set_yticks(range(len(heat_norm.index)))
    axes[1].set_yticklabels(heat_norm.index, fontsize=7)
    axes[1].set_title("B  Top reference-gap genera (row-normalized burden)")
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04, label="relative burden")
    amp = f" (amplicon={args.amplicon})" if args.amplicon else ""
    fig.suptitle(f"Geographic reference gap in bacterial oral 16S{amp}", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    base = os.path.join(args.outdir, "reference_gap_by_continent")
    fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig)
    log("wrote reference_gap_by_continent.png / .pdf")


if __name__ == "__main__":
    main()

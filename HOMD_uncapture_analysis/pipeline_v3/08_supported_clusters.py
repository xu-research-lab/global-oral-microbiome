#!/usr/bin/env python3
"""
pipeline_v3 Step 08: supported-cluster catalog (v2 14 equivalent).

Definition (user, 2026-08-18, adapted to TreeCluster units): clusters passing
the support gate
    max_subject_prevalence > 50 & n_bioprojects_total >= 2 & n_continents >= 2
AND clean_cluster (no contamination-flagged member) AND no organelle
(mitochondrion/chloroplast) member - the SUPPORTED CATALOG
(cross-CONTINENT evidence). Singleton clusters (size-1 TreeCluster units) are
included - each singleton ASV is its own cluster. The HEADLINE layer =
supported catalog + HOMD-uncaptured (no tree-distance filter, per user
2026-08-29; the headline figure/table shows the whole supported catalog).

Outputs (pipeline_v3/):
  - 08_supported_clusters.tsv   : the supported catalog (one row per cluster)
  - 08_supported_summary.json
  - analysis_report/figs/fig_supported.png : composition figure
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V3 = BASE / "pipeline_v3"
FIG = PIPE_V3 / "analysis_report" / "figs"
FIG.mkdir(parents=True, exist_ok=True)

SUPPORT = {"subject": 50, "bioprojects": 2, "continents": 2, "subject_gt": True,
           "bioprojects_ge": True, "continents_ge": True}
DIST_SPECIES_99 = 0.01   # ~99% identity species boundary on the GG2 tree

BLUE = "#4E79A7"; RED = "#E15759"; GRAY = "#BAB0AC"; TEAL = "#76B7B2"; INK = "#2B2B2B"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.edgecolor": "#BBBBBB", "axes.linewidth": 0.8,
                     "axes.labelcolor": INK, "text.color": INK, "figure.dpi": 150})
log = lambda s: print(f"[08_supported] {s}", flush=True)


def main():
    log("Loading aggregated clusters + per-ASV distances ...")
    agg = pd.read_csv(PIPE_V3 / "02_clusters_aggregated.tsv", sep="\t", low_memory=False)
    asv = pd.read_csv(PIPE_V3 / "asv_clusters.tsv", sep="\t", low_memory=False,
                      usecols=["asv_id", "cluster_id", "nearest_ref_dist", "covered"])
    asv["asv_id"] = asv["asv_id"].astype(str)
    # rep ASV = max subject prevalence member; recover rep distance from asv table
    rep_map = dict(zip(agg["rep_derep_id"].astype(str), agg["cluster_id"]))
    asv["is_rep"] = asv.apply(lambda r: rep_map.get(r["asv_id"]) == r["cluster_id"], axis=1)
    rep_dist = asv[asv["is_rep"]].set_index("cluster_id")["nearest_ref_dist"]
    agg["rep_nearest_ref_dist"] = agg["cluster_id"].map(rep_dist)

    for c in ["max_subject_prevalence", "n_bioprojects_total", "n_regions"]:
        agg[c] = pd.to_numeric(agg[c], errors="coerce")
    agg["clean_cluster"] = agg["clean_cluster"].astype(str).str.lower().isin(["true", "1"])
    agg["cluster_novelty_class"] = agg["cluster_novelty_class"].fillna("homd_uncaptured")

    m = ((agg["max_subject_prevalence"] > 50)
         & (agg["n_bioprojects_total"] >= 2)
         & (agg["n_continents"] >= 2)
         & agg["clean_cluster"])
    sup = agg[m].copy()
    # exclude organelle (mitochondrion / chloroplast) clusters
    if "any_organelle" in agg.columns:
        n_org_excl = int(sup["any_organelle"].fillna(False).sum())
        sup = sup[~sup["any_organelle"].fillna(False)]
        log(f"  organelle-flagged clusters excluded: {n_org_excl}")
    # headline = supported + HOMD-uncaptured (no tree-distance filter)
    sup["headline"] = (sup["cluster_novelty_class"] == "homd_uncaptured")
    sup["tree_tier"] = np.where(sup["rep_nearest_ref_dist"] <= 0.01, "<=0.01 (~99%)",
                       np.where(sup["rep_nearest_ref_dist"] <= 0.03, "0.01-0.03 (97-99%)",
                                ">0.03 (<97%)"))
    sup = sup.sort_values(["max_subject_prevalence", "n_bioprojects_total", "n_regions"],
                          ascending=False)
    sup.to_csv(PIPE_V3 / "08_supported_clusters.tsv", sep="\t", index=False)
    log(f"  saved 08_supported_clusters.tsv ({len(sup):,} rows)")

    # composition figure: tree tier x HOMD novelty
    ctab = pd.crosstab(sup["tree_tier"], sup["cluster_novelty_class"])
    order = [">0.03 (<97%)", "0.01-0.03 (97-99%)", "<=0.01 (~99%)"]
    ctab = ctab.reindex([c for c in order if c in ctab.index])

    fig, ax = plt.subplots(figsize=(8, 4.4))
    y = np.arange(len(ctab))
    left = np.zeros(len(ctab))
    for cls, color in [("homd_captured", GRAY), ("homd_uncaptured", RED)]:
        if cls not in ctab.columns:
            continue
        vals = ctab[cls].to_numpy()
        ax.barh(y, vals, left=left, color=color, height=0.58,
                label="HOMD-captured" if cls == "homd_captured" else "HOMD-uncaptured")
        left += vals
    for yi, (idx, row) in enumerate(ctab.iterrows()):
        total = int(row.sum())
        ax.text(total + 2, yi, f"{total}", va="center", fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels(ctab.index, fontsize=9)
    ax.set_xlim(0, ctab.sum(axis=1).max() * 1.15)
    ax.set_xlabel("Supported clusters")
    ax.set_title(f"{len(sup)} clean supported clusters (subjects >50, BioProjects ≥2, continents ≥2)\n"
                 "by tree distance and HOMD status (all supported clusters, incl. singletons)",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.7)
    fig.tight_layout()
    fig.savefig(FIG / "fig_supported.png", bbox_inches="tight")
    plt.close(fig)
    log("  saved fig_supported.png")

    summary = {
        "support_cutoffs": SUPPORT,
        "n_supported_clean": int(len(sup)),
        "by_tree_tier": sup["tree_tier"].value_counts().to_dict(),
        "by_novelty": sup["cluster_novelty_class"].value_counts().to_dict(),
        "n_headline_among_supported": int(sup["headline"].sum()),
        "headline_cutoffs": {"novelty": "homd_uncaptured",
                             "rep_nearest_ref_dist_max": "no tree-distance filter"},
        "organelle_excluded": "organelle-flagged (mitochondrion/chloroplast) clusters removed",
        "singletons_included": "size-1 TreeCluster units are their own clusters",
    }
    with open(PIPE_V3 / "08_supported_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log(f"  summary: {summary}")
    log("Done!")


if __name__ == "__main__":
    main()

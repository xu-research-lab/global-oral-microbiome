#!/usr/bin/env python3
"""
pipeline_v3 Step 05: Phylum-level distribution of candidate novel clusters
(v2 08 equivalent, cluster units).

Maps each candidate novel cluster to phylum/class via its annotation genus
(nearest-reference genus from Step 01) through the GG2 2022.10 taxonomy.

Input : 02_clusters_aggregated.tsv (candidate_novel subset)
Output: 05_candidate_novel_with_phylum.tsv + analysis_report/figs/fig_cluster_phylum.png
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
PIPE_V2 = BASE / "pipeline_v2"
OUT = PIPE_V3 / "analysis_report" / "figs"
OUT.mkdir(parents=True, exist_ok=True)
TAX = PIPE_V2 / "gappa" / "gg2_2022.10_tree_taxonomy.tsv"

BLUE = "#4E79A7"
RED = "#E15759"
GRAY = "#BAB0AC"
TEAL = "#76B7B2"
INK = "#2B2B2B"
GRID = "#E6E6E6"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.edgecolor": "#BABABA",
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "text.color": INK,
    "figure.dpi": 150,
})

log = lambda s: print(f"[05_geo] {s}", flush=True)


def parse_ranks(taxon_str):
    out = {}
    for part in str(taxon_str).split(";"):
        part = part.strip()
        if len(part) >= 3 and part[1:3] == "__":
            out[part[0]] = part[3:]
    return out


def main():
    log("Loading candidate novel clusters ...")
    agg = pd.read_csv(PIPE_V3 / "02_clusters_aggregated.tsv", sep="\t", low_memory=False)
    cn = agg[agg["candidate_novel"] == True].copy()
    log(f"  {len(cn):,} candidate novel clusters")

    log("Loading GG2 taxonomy ...")
    tax = pd.read_csv(TAX, sep="\t", low_memory=False, dtype=str, header=None,
                      names=["tip", "Taxon"])
    genus_to_phylum, genus_to_class = {}, {}
    base_to_phylum, base_to_class = {}, {}
    for _, row in tax.iterrows():
        ranks = parse_ranks(row["Taxon"])
        g = ranks.get("g", "")
        if not g:
            continue
        ph = ranks.get("p", "Unknown")
        cl = ranks.get("c", "Unknown")
        genus_to_phylum.setdefault(g, ph)
        genus_to_class.setdefault(g, cl)
        base = g.split("_")[0]
        base_to_phylum.setdefault(base, ph)
        base_to_class.setdefault(base, cl)
    log(f"  taxonomy entries: {len(tax):,} | genus lookups: {len(genus_to_phylum):,}")

    def lookup_ranks(genus):
        genus = str(genus)
        if genus in genus_to_phylum:
            return genus_to_phylum[genus], genus_to_class[genus]
        base = genus.split("_")[0]
        if base in base_to_phylum:
            return base_to_phylum[base], base_to_class[base]
        return "Unknown", "Unknown"

    log("Mapping cluster phylum via member ASV nearest-reference taxonomy ...")
    asv = pd.read_csv(PIPE_V3 / "asv_clusters.tsv", sep="\t", low_memory=False,
                      usecols=["asv_id", "cluster_id", "nearest_ref"])
    asv["asv_id"] = asv["asv_id"].astype(str)
    asv["nearest_ref"] = asv["nearest_ref"].astype(str)
    tip_tax = dict(zip(tax["tip"].astype(str), tax["Taxon"]))
    asv["taxon"] = asv["nearest_ref"].map(tip_tax)
    asv["gg2_phylum"] = asv["taxon"].map(lambda t: parse_ranks(t).get("p", "Unknown") if pd.notna(t) else "Unknown")
    asv["gg2_class"] = asv["taxon"].map(lambda t: parse_ranks(t).get("c", "Unknown") if pd.notna(t) else "Unknown")
    phy_maj = asv.groupby("cluster_id")["gg2_phylum"].agg(lambda x: x.value_counts().idxmax())
    cls_maj = asv.groupby("cluster_id")["gg2_class"].agg(lambda x: x.value_counts().idxmax())
    cn["gg2_phylum"] = cn["cluster_id"].map(phy_maj).fillna("Unknown")
    cn["gg2_class"] = cn["cluster_id"].map(cls_maj).fillna("Unknown")

    n_mapped = (cn["gg2_phylum"] != "Unknown").sum()
    log(f"  Mapped to phylum: {n_mapped:,} / {len(cn):,} ({n_mapped/len(cn)*100:.1f}%)")

    log("\nPhylum distribution (top 20):")
    phylum_counts = cn["gg2_phylum"].value_counts()
    for ph, n in phylum_counts.head(20).items():
        log(f"  {ph}: {n:,}")

    # oral clade: annotation genus in HOMD genus set
    homd_genus = set()
    with open("/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta") as f:
        for line in f:
            if line.startswith(">"):
                parts = line.split()
                if len(parts) >= 3:
                    homd_genus.add(parts[2].strip())
    cn["is_oral_clade"] = cn["annotation"].map(
        lambda g: str(g).split("_")[0] in homd_genus if pd.notna(g) and str(g).strip() else False)

    # Fig: phylum distribution split by oral clade
    TOP_N = 15
    top_phyla = phylum_counts.head(TOP_N).index.tolist()
    cn_plot = cn.copy()
    cn_plot["phylum_plot"] = cn_plot["gg2_phylum"].apply(
        lambda p: p if p in top_phyla else "Other")
    phylum_plot_order = cn_plot["phylum_plot"].value_counts().index.tolist()
    phylum_plot_order = [p for p in phylum_plot_order if p != "Other"] + \
                        (["Other"] if "Other" in phylum_plot_order else [])

    out_by_ph = cn_plot[~cn_plot["is_oral_clade"]]["phylum_plot"].value_counts()
    in_by_ph = cn_plot[cn_plot["is_oral_clade"]]["phylum_plot"].value_counts()

    fig, ax = plt.subplots(figsize=(9, 6))
    y = np.arange(len(phylum_plot_order))
    h = 0.38
    ax.barh(y + h / 2, [out_by_ph.get(p, 0) for p in phylum_plot_order],
            height=h, color=RED, label="outside HOMD genus")
    ax.barh(y - h / 2, [in_by_ph.get(p, 0) for p in phylum_plot_order],
            height=h, color=BLUE, label="in HOMD genus (oral clade)")
    ax.set_yticks(y)
    ax.set_yticklabels(phylum_plot_order, fontsize=9)
    ax.set_xlabel("Candidate novel clusters")
    ax.set_title(f"Phylum distribution of {len(cn):,} candidate novel clusters\n"
                 "(TreeCluster t=0.03, HOMD-uncaptured, subjects ≥10)",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_cluster_phylum.png", bbox_inches="tight")
    plt.close(fig)
    log(f"  saved {OUT / 'fig_cluster_phylum.png'}")

    cn.to_csv(PIPE_V3 / "05_candidate_novel_with_phylum.tsv", sep="\t", index=False)
    log("  saved 05_candidate_novel_with_phylum.tsv")

    summary = {
        "n_candidate_novel": int(len(cn)),
        "n_mapped_to_phylum": int(n_mapped),
        "phylum_counts": {k: int(v) for k, v in phylum_counts.head(20).items()},
        "oral_clade": {"in_homd_genus": int(cn["is_oral_clade"].sum()),
                       "outside_homd_genus": int((~cn["is_oral_clade"]).sum())},
    }
    with open(PIPE_V3 / "05_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log("  saved 05_summary.json")
    log("Done!")


if __name__ == "__main__":
    main()

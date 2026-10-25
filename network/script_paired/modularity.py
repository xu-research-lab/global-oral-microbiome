#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
modularity.py — Modularity analysis of the four HC/PD/DC/OSCC pooled genus
backbones, mirroring scripts_disease/modularity_disease.py.

Reads the outputs of script_hc_pd_dc_oscc/backbone_pooling.py:
  result_paired/backbone_v2/<grp>_genus_uw/backbone_edges.csv

Per group it reports:
  Q, Q_sd   mean/sd Louvain modularity over --seeds restarts (stochastic, so a
            single run is not reproducible on its own)
  k         median number of modules
  Q_null, z degree-preserving null (double-edge swap): sparse graphs have high
            Q by construction, so z = (Q - mean Q_null)/sd tells whether the
            module structure exceeds what the degree distribution already forces
  AMI, Q_phy agreement between the Louvain partition and phylum assignment
            (low values mean the modules are NOT a taxonomic artefact)

Outputs -> result_paired/modularity/:
  hc_pd_dc_oscc_modularity.csv  one row per group
  modules_<grp>.csv             node -> module assignment of the best-Q partition
  module_sizes.csv              module size table
  modularity_summary.png        Q vs null + module size distribution per group

Usage (from network/):
  python script_hc_pd_dc_oscc/modularity.py
"""

import argparse
import os

import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms.community import louvain_communities, modularity

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # network/
BACKBONE_DIR = os.path.join(ROOT, "result_paired", "backbone_v2")
OUT_DIR = os.path.join(ROOT, "result_paired", "modularity")
TAXONOMY_DEFAULT = os.environ.get("TAXONOMY_TSV",
                            "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv")

GROUPS = ["HC", "PD", "DC", "OSCC"]


# ---------------------------------------------------------------------------
# taxonomy helpers (same lookup policy as modularity_disease.py, genus level)
# ---------------------------------------------------------------------------

def load_taxonomy(path):
    if not path or not os.path.exists(path):
        return None
    tax = pd.read_csv(path, sep="\t", low_memory=False, index_col=0)
    tax = tax["Taxon"].str.split(";", expand=True).iloc[:, :7]
    tax.columns = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                   "Species"]
    return tax


def get_phylum_map(tax, fids):
    """Map genus ids to phyla: direct Feature-ID match, then the Genus column
    with its g__ prefix stripped (genus-name rows keep an accession id)."""
    _EMPTY = {"", "nan", "none", "p__", "d__", "g__"}
    if tax is None:
        return {f: "unknown" for f in fids}

    genus_to_phylum = {}
    if "Genus" in tax.columns and "Phylum" in tax.columns:
        g = tax["Genus"].dropna().astype(str).str.strip()
        p = tax.loc[g.index, "Phylum"].astype(str).str.strip()
        ok = ~g.isin(_EMPTY) & ~p.isin(_EMPTY)
        genus_to_phylum = dict(zip(
            g[ok].str.replace("^g__", "", regex=True), p[ok]))

    result = {}
    for fid in fids:
        if fid in tax.index:
            row = tax.loc[fid]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            pv = row.get("Phylum", "")
            if pd.notna(pv) and str(pv).strip() not in _EMPTY:
                result[fid] = str(pv).strip()
                continue
        if fid in genus_to_phylum:
            result[fid] = genus_to_phylum[fid]
        else:
            result[fid] = "unknown"
    return result


# ---------------------------------------------------------------------------
# core metrics
# ---------------------------------------------------------------------------

def louvain_stats(G, weight=None, seeds=50):
    """Mean/sd of Q and median k over `seeds` Louvain restarts."""
    qs, ks, parts = [], [], []
    for s in range(seeds):
        comms = louvain_communities(G, weight=weight, seed=s)
        qs.append(modularity(G, comms, weight=weight))
        ks.append(len(comms))
        parts.append(comms)
    best = parts[int(np.argmax(qs))]
    return float(np.mean(qs)), float(np.std(qs)), int(np.median(ks)), best


def null_modularity(G, seeds=30, swap_mult=10):
    """Q of degree-preserving rewired graphs (configuration-model style null)."""
    m = G.number_of_edges()
    qs = []
    for s in range(seeds):
        H = G.copy()
        try:
            nx.double_edge_swap(H, nswap=swap_mult * m,
                                max_tries=20 * swap_mult * m, seed=s)
        except (nx.NetworkXAlgorithmError, nx.NetworkXError):
            pass  # too few swappable edges; keep whatever swaps succeeded
        comms = louvain_communities(H, weight=None, seed=s)
        qs.append(modularity(H, comms, weight=None))
    return float(np.mean(qs)), float(np.std(qs))


def partition_vs_attribute(G, partition, attr_map):
    """AMI(partition, attribute) and the modularity the attribute itself gives."""
    from sklearn.metrics import adjusted_mutual_info_score

    nodes = list(G.nodes())
    lab = {}
    for i, comm in enumerate(partition):
        for n in comm:
            lab[n] = i
    a_lab = [attr_map.get(n, "unknown") for n in nodes]
    ami = float(adjusted_mutual_info_score([lab[n] for n in nodes], a_lab))

    by_attr = {}
    for n, a in zip(nodes, a_lab):
        by_attr.setdefault(a, set()).add(n)
    q_attr = float(modularity(G, list(by_attr.values()), weight=None))
    return ami, q_attr


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def run_group(grp, seeds, null_seeds):
    d = os.path.join(BACKBONE_DIR, f"{grp}_genus_uw")
    edges = pd.read_csv(os.path.join(d, "backbone_edges.csv"))
    G = nx.Graph()
    for _, e in edges.iterrows():
        G.add_edge(str(e.taxon_i), str(e.taxon_j),
                   w=abs(float(e.pooled_obs)), r=float(e.pooled_obs))
    if G.number_of_edges() == 0:
        print(f"[mod] {grp}: empty backbone, skipped")
        return None

    tax = load_taxonomy(TAXONOMY_DEFAULT)
    if tax is None:
        print(f"[warn] taxonomy not found at {TAXONOMY_DEFAULT}; "
              f"phyla = 'unknown'")
    phy = get_phylum_map(tax, sorted(G.nodes()))

    q_uw, sd_uw, k_uw, best = louvain_stats(G, None, seeds)
    q_w, sd_w, k_w, _ = louvain_stats(G, "w", seeds)
    q_null, sd_null = null_modularity(G, null_seeds)
    z = (q_uw - q_null) / sd_null if sd_null > 0 else np.nan
    ami, q_phy = partition_vs_attribute(G, best, phy)

    row = dict(
        group=grp, nodes=G.number_of_nodes(), edges=G.number_of_edges(),
        Q=round(q_uw, 4), Q_sd=round(sd_uw, 4), k=k_uw,
        Q_weighted=round(q_w, 4), k_weighted=k_w,
        Q_null=round(q_null, 4), Q_null_sd=round(sd_null, 4),
        z_vs_null=round(z, 2),
        AMI_module_phylum=round(ami, 4),
        Q_phylum_partition=round(q_phy, 4),
    )
    print(f"[mod] {grp}: {row['nodes']} nodes, {row['edges']} edges  "
          f"Q={q_uw:.3f} k={k_uw} z={z:.1f} AMI={ami:.3f}")

    # per-node module assignment of the best-Q partition (M1 = largest module)
    mod_of = {}
    for i, comm in enumerate(sorted(best, key=len, reverse=True)):
        for n in comm:
            mod_of[n] = i + 1
    deg = dict(G.degree())
    mods = pd.DataFrame([
        dict(taxon_id=n, module=mod_of[n], degree=deg[n],
             phylum=phy.get(n, "unknown"), taxon_name=n)
        for n in sorted(G.nodes(), key=lambda x: (mod_of[x], -deg[x]))
    ])
    mods.to_csv(os.path.join(OUT_DIR, f"modules_{grp}.csv"), index=False)
    sizes = mods.groupby("module").size()
    print(f"[mod] {grp} module sizes:", sizes.to_dict())
    return row, sizes


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--null-seeds", type=int, default=30)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    rows, size_frames = [], []
    for grp in GROUPS:
        if not os.path.exists(os.path.join(BACKBONE_DIR, f"{grp}_genus_uw",
                                           "backbone_edges.csv")):
            print(f"[mod] {grp}: no backbone_edges.csv (pooling not run), "
                  f"skipped")
            continue
        got = run_group(grp, args.seeds, args.null_seeds)
        if got is None:
            continue
        row, sizes = got
        rows.append(row)
        for m, n in sizes.items():
            size_frames.append(dict(group=grp, module=m, n_genera=int(n)))

    pd.DataFrame(rows).to_csv(
        os.path.join(OUT_DIR, "hc_pd_dc_oscc_modularity.csv"), index=False)
    pd.DataFrame(size_frames).to_csv(
        os.path.join(OUT_DIR, "module_sizes.csv"), index=False)

    # summary figure: per group Q vs null + module size distribution
    rows = [r for r in rows if r]
    n = len(rows)
    if n == 0:
        print("[mod] nothing to plot")
        return
    fig, axes = plt.subplots(2, n, figsize=(4.6 * n, 7.6))
    if n == 1:
        axes = axes.reshape(2, 1)
    for i, r in enumerate(rows):
        ax = axes[0, i]
        # this matplotlib cannot take per-bar color lists with yerr (errorbar
        # builds one Line2D), so draw the two bars separately
        ax.bar([0], [r["Q"]], 0.55, yerr=[r["Q_sd"]], capsize=3,
               color="#2c7fb8", ecolor="#08519c")
        ax.bar([1], [r["Q_null"]], 0.55, yerr=[r["Q_null_sd"]], capsize=3,
               color="#bdbdbd", ecolor="#636363")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["observed", "null"], fontsize=8)
        ax.text(0.5, max(r["Q"], r["Q_null"]) + 0.03,
                f"z={r['z_vs_null']:.0f}", ha="center", fontsize=9)
        ax.set_ylabel("modularity Q")
        ax.set_ylim(0, max(r["Q"], r["Q_null"]) * 1.35)
        ax.set_title(f"{r['group']}: Q vs null", fontsize=10)
        ax = axes[1, i]
        sub = [s for s in size_frames if s["group"] == r["group"]]
        sub = sorted(sub, key=lambda s: s["module"])
        ax.bar([f"M{s['module']}" for s in sub],
               [s["n_genera"] for s in sub], color="#2c7fb8", width=0.8)
        ax.set_xlabel("module (M1 = largest)")
        ax.set_ylabel("# genera")
        ax.set_title(f"{r['group']}: module sizes", fontsize=10)
        ax.tick_params(axis="x", labelsize=7, rotation=45)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "modularity_summary.png"), dpi=200)
    plt.close(fig)
    print(f"[done] -> {OUT_DIR}")


if __name__ == "__main__":
    main()
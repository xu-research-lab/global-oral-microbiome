#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
backbone_pooling.py — HC/PD/DC/OSCC consensus backbones (NetMoss uw),
paired-cohort scope (run_paired_downstream.py patches RESULT/SPARCC_BASE
to the paired directories).

Reads the per-(group, cohort) genus SparCC matrices under
  result_paired/per_site_sparcc/<grp>_genus/<study>.csv
and pools them per group (HC, PD, DC, OSCC) onto the union genus space:
  - zero-padded pooling            -> consensus sign reference
    (performed by the OFFICIAL NetMoss2::getNetwork R package, v0.2.0,
    Xiao et al. 2022 Nat Comput Sci, via netmoss_get_network.R)
  - observation-conditional pooling -> edge weight (pooled_obs)

Pooling unit = (group, bioproject) paired cohort scope; units are prepared
by prepare_paired_units.py. Gates applied uniformly to all four groups:
  uw weight   w_k = (n_k - 1) / (1 - r_k^2)          (clip 1 - r^2 at 1e-8)
  three gates coverage >= k_min(5)  AND  repro >= 0.60  AND  |pooled_obs| >= 0.20
  binomial test (report-only)  p0 = 0.025, BH-FDR
  LOO diagnostics (pooled_wo, top_weight_share) kept as columns, not filters

Outputs -> result_paired/backbone_v2/<grp>_genus_uw/:
  candidates.csv, backbone_edges.csv, reproducibility.png, backbone.png,
  pooled_corr_distribution.png, units.csv, counts.json
plus result_paired/backbone_summary/hc_pd_dc_oscc_backbone_counts.csv

Usage (from network/, normally via run_paired_figures.py / run_paired_downstream.py):
  python script_paired/backbone_pooling.py
  python script_paired/backbone_pooling.py --k-min 3   # sensitivity
"""

import argparse
import hashlib
import json
import os
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # network/
RESULT = os.path.join(ROOT, "result_paired")
SPARCC_BASE = os.path.join(RESULT, "per_site_sparcc")
OUT_BASE = os.path.join(RESULT, "backbone_v2")
SUM_DIR = os.path.join(RESULT, "backbone_summary")

# official NetMoss2 getNetwork pooling (zero-padded consensus reference)
RSCRIPT = os.environ.get("RSCRIPT",
                         "/home/zhangyuchen/.conda/envs/r_latest/bin/Rscript")
NETMOSS_WRAPPER = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "netmoss_get_network.R")

GROUPS = ["HC", "PD", "DC", "OSCC"]
WEIGHT = "uw"
EPS = 1e-8


# ---------------------------------------------------------------------------
# helpers (identical semantics to backbone_pooling_disease.py)
# ---------------------------------------------------------------------------

def bh_fdr(p):
    """Benjamini-Hochberg FDR (NaN preserved)."""
    p = np.asarray(p, dtype=float).ravel()
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv, kind="mergesort")
    pv = pv[order]
    adj = pv * (n / np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    inv = np.empty(n, dtype=np.intp)
    inv[order] = np.arange(n)
    out[ok] = adj[inv]
    return out


def weight_block(r, n):
    """uw: w = (n-1) / (1 - r^2), clipped to avoid divergence as |r| -> 1."""
    return (n - 1) / np.clip(1.0 - r ** 2, EPS, 1.0)


def load_unit_matrix(unit_id, sid, pos):
    fp = os.path.join(SPARCC_BASE, unit_id, f"{sid}.csv")
    if not os.path.exists(fp):
        return None
    df = pd.read_csv(fp, sep=",", index_col=0)
    df.index = [str(x) for x in df.index]
    df.columns = [str(x) for x in df.columns]
    inter = [f for f in df.index if f in pos]
    if len(inter) < 2:
        return None
    df = df.loc[inter, inter]
    idx = np.fromiter((pos[f] for f in inter), dtype=np.int64, count=len(inter))
    return idx, df.to_numpy(dtype=np.float64)


# ---------------------------------------------------------------------------
# unit discovery
# ---------------------------------------------------------------------------

def discover_units(grp, units_all):
    """Units of one group with a SparCC matrix present."""
    units = []
    for _, u in units_all[units_all.category == grp].iterrows():
        fp = os.path.join(SPARCC_BASE, u.unit_id, f"{u.study}.csv")
        if os.path.exists(fp):
            units.append(dict(unit_id=u.unit_id, study=str(u.study),
                              n=int(u.n_samples), path=fp))
        else:
            print(f"[warn] no SparCC matrix for {u.unit_id}/{u.study}, skipped")
    return units


def collect_fid(units):
    """Union genus space across the group's units."""
    fid = set()
    for u in units:
        df = pd.read_csv(u["path"], sep=",", index_col=0, nrows=0)
        fid.update(str(x) for x in df.index)
        fid.update(str(x) for x in df.columns)
    return np.array(sorted(fid))


# ---------------------------------------------------------------------------
# core accumulation (ported from backbone_pooling_disease.py, keyed by unit)
# ---------------------------------------------------------------------------

def build_stats(units, fid):
    F = len(fid)
    pos = {f: k for k, f in enumerate(fid)}

    W1_all = np.zeros((F, F))   # zero-padded weight sum        (consensus)
    VR_all = np.zeros((F, F))   # zero-padded weight*value sum  (consensus)
    W1 = np.zeros((F, F))       # observed weight sum
    W2 = np.zeros((F, F))       # observed weight^2 sum
    VR = np.zeros((F, F))       # observed weight*value sum
    VR2 = np.zeros((F, F))      # observed weight*value^2 sum
    cov = np.zeros((F, F), dtype=np.int32)

    cohorts = []
    for u in units:
        n = u["n"]
        if n < 2:
            continue
        got = load_unit_matrix(u["unit_id"], u["study"], pos)
        if got is None:
            continue
        idx, r = got
        w = weight_block(r, n)

        present = np.zeros(F, dtype=bool)
        present[idx] = True
        obs = np.outer(present, present)
        obsb = obs.astype(np.float64)

        v_full = np.zeros((F, F))
        w_full = np.zeros((F, F))
        v_full[np.ix_(idx, idx)] = r            # uw: pooling space == r space
        w_full[np.ix_(idx, idx)] = w

        w_full[~obs] = n - 1                    # zero-pad missing pairs
        np.fill_diagonal(w_full, 1.0)           # self-correlation weight 1
        np.fill_diagonal(v_full, 0.0)

        W1_all += w_full
        VR_all += w_full * v_full
        W1 += w_full * obsb
        W2 += w_full ** 2 * obsb
        VR += w_full * v_full * obsb
        VR2 += w_full * v_full ** 2 * obsb
        cov += obs

        cohorts.append((idx, r, w, n, u))
        print(f"  [unit] {u['study']}  n={n}  genera={len(idx)}")

    return dict(W1_all=W1_all, VR_all=VR_all, W1=W1, W2=W2, VR=VR, VR2=VR2,
                cov=cov), cohorts


def consensus_sign(acc):
    with np.errstate(invalid="ignore", divide="ignore"):
        cons = np.where(acc["W1_all"] > 0, acc["VR_all"] / acc["W1_all"], 0.0)
    cs = np.sign(cons)
    np.fill_diagonal(cs, 1.0)
    return cs


# ---------------------------------------------------------------------------
# NetMoss::getNetwork consensus (zero-padded pooling, official R package)
# ---------------------------------------------------------------------------

def consensus_pool_py(acc):
    """In-house zero-padded pooling — kept ONLY as the patch source for the
    saturated (non-finite) pairs that NetMoss2 reports as NaN; all results
    come from the official package."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(acc["W1_all"] > 0, acc["VR_all"] / acc["W1_all"], 0.0)


def _manifest_text(units):
    rows = pd.DataFrame([{"study": u["study"],
                          "matrix_path": os.path.abspath(u["path"]),
                          "n": int(u["n"])} for u in units])
    return rows.to_csv(index=False)


def netmoss_pooled_matrix(units, fid, rscript, grp):
    """Zero-padded pooled r over fid from the official NetMoss::getNetwork.

    Writes the manifest to its cache dir and reuses the pooled matrix when
    the manifest is unchanged. Returns (cons_matrix F x F, meta dict).
    """
    out_dir = os.path.join(RESULT, "backbone_netmoss", f"{grp}_genus")
    os.makedirs(out_dir, exist_ok=True)
    man = os.path.join(out_dir, "manifest.csv")
    man_text = _manifest_text(units)
    pooled_fp = os.path.join(out_dir, "pooled_union.csv")
    run_meta = {}
    meta_fp = os.path.join(out_dir, "netmoss_run.json")
    if os.path.isfile(meta_fp):
        try:
            with open(meta_fp) as f:
                run_meta = {k: v for k, v in json.load(f).items()}
        except (json.JSONDecodeError, OSError):
            run_meta = {}
    # cache is valid only for the same engine + same manifest
    cached = (os.path.isfile(man) and os.path.isfile(pooled_fp)
              and run_meta.get("engine") == "netmoss2"
              and open(man).read() == man_text)
    if not cached:
        with open(man, "w") as f:
            f.write(man_text)
        cmd = [rscript, NETMOSS_WRAPPER, "--manifest", man, "--out", out_dir]
        try:
            r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                               timeout=7200)
        except FileNotFoundError as e:
            raise RuntimeError(f"Rscript not found: {e}") from None
        if r.returncode != 0 or "NETMOSS_OK" not in r.stdout:
            raise RuntimeError(
                f"NetMoss2::getNetwork failed (exit {r.returncode})\n"
                f"{r.stderr[-4000:]}")
        with open(meta_fp) as f:
            run_meta = {k: v for k, v in json.load(f).items()}
    else:
        print(f"[{grp}] netmoss2 cache hit ({out_dir})")

    cons = pd.read_csv(pooled_fp, index_col=0)
    cons.index = cons.index.astype(str)
    cons.columns = cons.columns.astype(str)
    missing = set(fid) - set(cons.index)
    if missing:
        raise RuntimeError(f"pooled matrix missing {len(missing)} taxa, e.g. "
                           f"{sorted(missing)[:5]}")
    cons = cons.reindex(index=list(fid), columns=list(fid)).to_numpy()
    return cons, run_meta


def resolve_consensus(acc, units, fid, rscript, grp):
    """Consensus pooled matrix for sign reference (official NetMoss2 only).

    Non-finite NetMoss entries (saturated per-cohort r, no eps clip in the
    package) are patched with the in-house clipped value and counted.
    """
    cons_py = consensus_pool_py(acc)   # NaN-patch source (see its docstring)
    cons_r, meta = netmoss_pooled_matrix(units, fid, rscript, grp)
    engine = "netmoss2::getNetwork" \
        if meta.get("engine") == "netmoss2" else "netmoss::getNetwork"
    nf = ~np.isfinite(cons_r) & ~np.eye(len(cons_r), dtype=bool)
    n_sub = int(np.count_nonzero(nf))
    cons_r[nf] = cons_py[nf]          # in-house clipped value, mirrored side too
    off = ~np.eye(len(cons_r), dtype=bool)
    diff = cons_r[off] - cons_py[off]
    meta["py_vs_r_max_abs_diff"] = float(np.max(np.abs(diff))) if diff.size else 0.0
    meta["n_nonfinite_substituted_pairs"] = n_sub // 2
    meta.pop("n_substituted_nonfinite", None)
    print(f"[{grp}] netmoss vs in-house regression check: max|d|="
          f"{meta['py_vs_r_max_abs_diff']:.2e}, non-finite patched={n_sub // 2} pairs")
    return cons_r, engine, meta


def count_confirms(cohorts, F, cons_sign, r_thr):
    n_confirm = np.zeros((F, F), dtype=np.int32)
    for idx, r, w, n, u in cohorts:
        r_full = np.zeros((F, F))
        r_full[np.ix_(idx, idx)] = r
        strong = (np.abs(r_full) >= r_thr)
        confirm = strong & (np.sign(r_full) == cons_sign) & (cons_sign != 0)
        np.fill_diagonal(confirm, False)
        n_confirm += confirm
    return n_confirm


def build_tables(acc, cohorts, F, fid, r_thr, k_min, repro_thr, consensus_thr,
                 p0, cons_sign=None):
    if cons_sign is None:
        cons_sign = consensus_sign(acc)
    n_confirm = count_confirms(cohorts, F, cons_sign, r_thr)

    iu = np.triu_indices(F, k=1)
    cov_u = acc["cov"][iu]
    W1_u = acc["W1"][iu]
    VR_u = acc["VR"][iu]

    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = np.where(W1_u > 0, VR_u / W1_u, np.nan)

    # LOO diagnostics: drop the single highest-weight cohort per pair
    # (the max-|r| cohort carries the max uw weight)
    rmax = np.zeros((F, F))
    w_at = np.zeros((F, F))
    v_at = np.zeros((F, F))
    for idx, r, w, n, u in cohorts:
        ix = np.ix_(idx, idx)
        rabs = np.abs(r)
        blk = rmax[ix]
        newmax = rabs > blk
        rmax[ix] = np.where(newmax, rabs, blk)
        w_at[ix] = np.where(newmax, w, w_at[ix])
        v_at[ix] = np.where(newmax, r, v_at[ix])
    W_at_u = w_at[iu]
    V_at_u = v_at[iu]

    with np.errstate(invalid="ignore", divide="ignore"):
        pooled_wo = np.where(W1_u - W_at_u > 0,
                             (VR_u - W_at_u * V_at_u) / (W1_u - W_at_u), np.nan)

    nconf_u = n_confirm[iu]
    rep_u = np.where(cov_u > 0, nconf_u / cov_u, np.nan)

    cand = cov_u >= k_min
    base = cand & (rep_u >= repro_thr) & (np.abs(pooled) >= consensus_thr)
    final = base

    binom_p = np.full(F * (F - 1) // 2, np.nan)
    binom_p[cand] = stats.binom.sf(nconf_u[cand] - 1, cov_u[cand], p0)
    binom_p_adj = np.full_like(binom_p, np.nan)
    if cand.any():
        binom_p_adj[cand] = bh_fdr(binom_p[cand])

    rows = np.flatnonzero(cand)
    ii, jj = iu[0][rows], iu[1][rows]
    df = pd.DataFrame({
        "taxon_i": fid[ii], "taxon_j": fid[jj],
        "pooled_obs": pooled[rows],
        "coverage": cov_u[rows].astype(int),
        "n_confirm": nconf_u[rows].astype(int),
        "repro": rep_u[rows],
        "binom_p": binom_p[rows],
        "binom_p_adj": binom_p_adj[rows],
        "top_weight_share": np.where(W1_u[rows] > 0,
                                     W_at_u[rows] / W1_u[rows], np.nan),
        "pooled_wo": pooled_wo[rows],
        "baseline": base[rows],
        "final": final[rows],
    })

    counts = dict(baseline=int(base.sum()), final=int(final.sum()))
    return df, counts, pooled, iu


# ---------------------------------------------------------------------------
# plots (same presentation as the disease pipeline)
# ---------------------------------------------------------------------------

def plot_repro(df, out):
    if len(df) == 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "no candidate pairs", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
        return
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].hist(df.repro, bins=40, color="#7f8c8d", edgecolor="black")
    ax[0].axvline(np.nanmedian(df.repro), color="red", ls="--",
                  label=f"median={np.nanmedian(df.repro):.2f}")
    ax[0].set_title("reproducibility (candidate pairs)"); ax[0].legend(frameon=False)
    ax[1].hist(df.coverage, bins=range(1, int(df.coverage.max()) + 2),
               color="#95a5a6", edgecolor="black", align="left")
    ax[1].set_title("coverage distribution")
    fig.tight_layout(); fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_backbone(df, out, grp):
    if len(df) == 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "empty backbone", ha="center"); ax.axis("off")
        fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
        return
    G = nx.Graph()
    for _, e in df.iterrows():
        G.add_edge(e.taxon_i, e.taxon_j, weight=float(e.pooled_obs))
    pos = nx.spring_layout(G, seed=42, k=1.2 / np.sqrt(max(len(G), 1)),
                           iterations=120)
    ecol = ["#c0392b" if G[u][v]["weight"] > 0 else "#2980b9"
            for u, v in G.edges()]
    ew = [0.3 + 2.0 * abs(G[u][v]["weight"]) for u, v in G.edges()]
    fig, ax = plt.subplots(figsize=(11, 8))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=ecol, width=ew, alpha=0.4)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=150, node_color="#2c3e50",
                           linewidths=0)
    legend_handles = [Line2D([0], [0], color="#c0392b", lw=2, label="positive"),
                      Line2D([0], [0], color="#2980b9", lw=2, label="negative")]
    ax.legend(handles=legend_handles, loc="lower right", frameon=False)
    ax.set_title(f"{grp} backbone: {G.number_of_nodes()} taxa, "
                 f"{G.number_of_edges()} edges")
    ax.axis("off"); fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_dist(pooled, out, grp):
    z = pooled[np.isfinite(pooled)]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(z, bins=100, color="skyblue", edgecolor="black", log=True)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_title(f"{grp} pooled correlation (zero-padded, uw)")
    ax.set_xlabel("pooled corr"); ax.set_ylabel("pairs (log)")
    fig.tight_layout(); fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k-min", type=int, default=5)
    ap.add_argument("--r-thr", type=float, default=0.20)
    ap.add_argument("--repro-thr", type=float, default=0.60)
    ap.add_argument("--consensus-thr", type=float, default=0.20)
    ap.add_argument("--p0", type=float, default=0.025, help="binomial null rate")
    ap.add_argument("--tag", default="", help="output subdir suffix (e.g. k3)")
    ap.add_argument("--rscript", default=RSCRIPT,
                    help=f"Rscript binary for the NetMoss wrapper")
    args = ap.parse_args()

    units_all = pd.read_csv(os.path.join(RESULT, "units.csv"),
                            dtype={"study": str})
    summary_rows = []
    for grp in GROUPS:
        print(f"\n===== {grp} =====")
        units = discover_units(grp, units_all)
        if not units:
            print(f"[{grp}] no units with SparCC matrices, skipped")
            continue
        print(f"[{grp}] {len(units)} units "
              f"({sum(u['n'] for u in units)} samples)")

        fid = collect_fid(units)
        print(f"[{grp}] union genus space: {len(fid)} genera")

        acc, cohorts = build_stats(units, fid)
        cons, engine, nm_meta = resolve_consensus(
            acc, units, fid, args.rscript, grp, )
        cons_sign = np.sign(cons)
        np.fill_diagonal(cons_sign, 1.0)
        df, counts, pooled, iu = build_tables(
            acc, cohorts, len(fid), fid,
            args.r_thr, args.k_min, args.repro_thr, args.consensus_thr,
            args.p0, cons_sign=cons_sign)
        print(f"[{grp}] consensus engine: {engine}")
        print(f"[{grp}] baseline={counts['baseline']}  final={counts['final']}")

        d = os.path.join(OUT_BASE, f"{grp}_genus_uw{args.tag}")
        os.makedirs(d, exist_ok=True)
        df.to_csv(os.path.join(d, "candidates.csv"), index=False)

        final_df = df.loc[df.final].copy()
        cols = ["taxon_i", "taxon_j", "pooled_obs", "coverage", "n_confirm",
                "repro", "binom_p", "binom_p_adj", "top_weight_share",
                "pooled_wo"]
        final_df[cols].to_csv(os.path.join(d, "backbone_edges.csv"), index=False)

        plot_repro(df, os.path.join(d, "reproducibility.png"))
        plot_backbone(final_df, os.path.join(d, "backbone.png"), grp)
        plot_dist(pooled, os.path.join(d, "pooled_corr_distribution.png"), grp)

        pd.DataFrame(units).to_csv(os.path.join(d, "units.csv"), index=False)
        with open(os.path.join(d, "counts.json"), "w") as f:
            json.dump({"group": grp, "weight": "uw",
                       "n_units": len(units),
                       "n_samples": int(sum(u["n"] for u in units)),
                       "n_genera": int(len(fid)),
                       "k_min": args.k_min, "repro_thr": args.repro_thr,
                       "consensus_thr": args.consensus_thr,
                       "pooling_engine": engine,
                       "pooling_detail": {k: v for k, v in nm_meta.items()
                                          if isinstance(v, (int, float, str,
                                                            bool))},
                       **counts},
                      f, indent=2)
        n_nodes = len(set(final_df.taxon_i) | set(final_df.taxon_j)) \
            if len(final_df) else 0
        summary_rows.append(dict(group=grp, weight=WEIGHT, pooling_engine=engine,
                                 n_units=len(units),
                                 n_samples=int(sum(u["n"] for u in units)),
                                 n_genera=int(len(fid)), **counts,
                                 n_nodes=n_nodes))

    os.makedirs(SUM_DIR, exist_ok=True)
    pd.DataFrame(summary_rows).to_csv(
        os.path.join(SUM_DIR, "hc_pd_dc_oscc_backbone_counts.csv"), index=False)
    print(f"[done] -> {OUT_BASE}")


if __name__ == "__main__":
    main()
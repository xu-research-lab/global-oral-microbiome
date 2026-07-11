#!/usr/bin/env python3
"""Closed-reference (GG2) mapping rate by continent, controlling for amplicon region.

Motivation
----------
The per-sample GG2 closed-reference mapping rate (`sample_cr_mapping_rate`) is used
as a QC metric. A raw by-continent comparison suggests that understudied / non-Western
regions map worse to Greengenes2 (a potential reference-representation bias). However,
mapping rate is strongly confounded by:
  * amplicon region (V1-V3 / V1-V2 map worse than V4 regardless of geography), and
  * study/protocol (mapping rate is largely a study-level property; continent is nearly
    nested within study, so sample-level tests massively pseudoreplicate).

This script therefore reports the continent effect three ways, from least to most rigorous:
  1. Raw per-sample distributions by continent (all amplicons).
  2. Amplicon-held-constant distributions (V4 only, the largest stratum).
  3. A STUDY-LEVEL weighted regression of logit(mapping rate) on continent + amplicon
     region, so the residual continent effect is quantified net of amplicon and without
     sample-level pseudoreplication.

Outputs (under results/):
  figures/05_homd_gg2_closed_reference_profiling/
      qc_cr_mapping_by_continent.{png,pdf}            # all amplicons, sample + study panels
      qc_cr_mapping_by_continent_V4only.{png,pdf}     # V4 only (amplicon held constant)
  tables/05_homd_gg2_closed_reference/
      cr_mapping_by_continent_summary.csv             # per-continent summary (all + V4)
      cr_mapping_continent_regression.csv             # study-level amplicon-adjusted model

Runnable locally against the committed QC metadata; no cluster inputs required.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_META = REPO_ROOT / "data/metadata-and-biom-v2-qc/metadata_all_for_qc_raw.tsv"
DEFAULT_FIGDIR = REPO_ROOT / "results/figures/05_homd_gg2_closed_reference_profiling"
DEFAULT_TABDIR = REPO_ROOT / "results/tables/05_homd_gg2_closed_reference"

ORAL_SITES = ["saliva", "tooth", "tongue", "oral mucosa", "oral swab", "gingiva"]
# Western-derived populations are well represented in Greengenes2; the contrast of interest
# is whether non-Western regions map worse after holding amplicon constant.
NON_WESTERN = {"Africa", "South America", "Asia"}
COL_WEST = "#4C72B0"
COL_NONWEST = "#C44E52"
QC_FLOOR = 0.30


def log(msg: str) -> None:
    print(f"[cr-continent] {msg}", flush=True)


def load_oral(meta_path: Path) -> pd.DataFrame:
    d = pd.read_csv(meta_path, sep="\t", low_memory=False)
    d["cr"] = pd.to_numeric(d["sample_cr_mapping_rate"], errors="coerce")
    o = d[
        d["site_cleaned"].isin(ORAL_SITES)
        & d["cr"].notna()
        & d["continent_group"].notna()
    ].copy()
    log(f"oral samples with CR mapping + continent: {len(o):,}")
    return o


def study_table(o: pd.DataFrame) -> pd.DataFrame:
    """Collapse to one row per BioProject: mean mapping rate, n, continent, dominant region."""
    def dominant(s: pd.Series):
        m = s.mode()
        return m.iloc[0] if len(m) else np.nan

    st = (
        o.groupby("bioproject_id")
        .agg(
            cr=("cr", "mean"),
            n=("cr", "size"),
            continent=("continent_group", dominant),
            region=("sequencing_region_cleaned", dominant),
        )
        .reset_index()
    )
    return st


def make_figure(o: pd.DataFrame, out_base: Path, title: str, subtitle_note: str) -> None:
    order = o["continent_group"].value_counts().index.tolist()
    colors = [COL_NONWEST if c in NON_WESTERN else COL_WEST for c in order]
    st = study_table(o)
    rng = np.random.RandomState(42)

    by_sample = {c: o[o["continent_group"] == c] for c in order}
    by_study = {c: st[st["continent"] == c] for c in order}

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

    def panel(ax, data_by, valcol, ttl, dot_size_col=None, max_dots=None):
        bp = ax.boxplot(
            [data_by[c][valcol].values for c in order],
            positions=range(len(order)), widths=0.6, showfliers=False,
            patch_artist=True, medianprops=dict(color="black", lw=1.6),
            whiskerprops=dict(color="#555"), capprops=dict(color="#555"),
            boxprops=dict(alpha=0.35),
        )
        for patch, col in zip(bp["boxes"], colors):
            patch.set_facecolor(col)
        for i, c in enumerate(order):
            sub = data_by[c]
            v = sub[valcol].values
            if max_dots and len(v) > max_dots:
                idx = rng.choice(len(v), max_dots, replace=False)
                v = v[idx]
            x = rng.normal(i, 0.07, size=len(v))
            if dot_size_col is not None:
                s = np.clip(np.sqrt(sub[dot_size_col].values) * 3, 8, 300)
                ax.scatter(x, v, s=s, color=colors[i], alpha=0.55,
                           edgecolor="white", linewidth=0.4, zorder=3)
            else:
                ax.scatter(x, v, s=6, color=colors[i], alpha=0.20,
                           edgecolor="none", zorder=3)
        ax.set_xticks(range(len(order)))
        ax.set_title(ttl, fontsize=12, pad=10)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(axis="y", alpha=0.25)
        ax.axhline(QC_FLOOR, ls="--", color="#888", lw=1)

    panel(axes[0], by_sample, "cr", "A  Per-sample mapping rate", max_dots=400)
    axes[0].set_xticklabels(
        [f"{c}\n{len(by_sample[c]):,} samples" for c in order], fontsize=9
    )
    axes[0].set_ylabel("GG2 closed-reference mapping rate", fontsize=11)

    panel(axes[1], by_study, "cr", "B  Per-study mean mapping rate", dot_size_col="n")
    axes[1].set_xticklabels(
        [f"{c}\n{len(by_study[c])} studies" for c in order], fontsize=9
    )
    axes[1].text(len(order) - 0.5, QC_FLOOR + 0.015, "0.30 QC floor",
                 ha="right", va="bottom", fontsize=8, color="#666")

    leg = [
        Patch(facecolor=COL_NONWEST, alpha=0.6, label="Non-Western (Africa, S. America, Asia)"),
        Patch(facecolor=COL_WEST, alpha=0.6, label="Western-referenced (NA, Europe, Oceania)"),
    ]
    axes[0].legend(handles=leg, loc="lower left", fontsize=8, framealpha=0.9)
    fig.suptitle(title, fontsize=13, y=0.99)
    fig.text(0.5, 0.005, subtitle_note, ha="center", fontsize=8, color="#555")
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])

    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_base}.png", dpi=200, bbox_inches="tight")
    fig.savefig(f"{out_base}.pdf", bbox_inches="tight")
    plt.close(fig)
    log(f"wrote {out_base}.png / .pdf")


def summarize(o: pd.DataFrame, label: str) -> pd.DataFrame:
    order = o["continent_group"].value_counts().index.tolist()
    st = study_table(o)
    rows = []
    for c in order:
        s = o.loc[o["continent_group"] == c, "cr"]
        nstud = int((st["continent"] == c).sum())
        rows.append(dict(
            subset=label, continent=c, n_samples=len(s), n_studies=nstud,
            median_cr=round(float(s.median()), 4),
            mean_unmapped_pct=round(float((1 - s.mean()) * 100), 2),
            pct_below_floor=round(float((s < QC_FLOOR).mean() * 100), 2),
        ))
    return pd.DataFrame(rows)


def logit(p: np.ndarray, eps: float = 1e-4) -> np.ndarray:
    p = np.clip(p, eps, 1 - eps)
    return np.log(p / (1 - p))


def _norm_sf(z: np.ndarray) -> np.ndarray:
    """Upper-tail standard-normal probability, via erfc (no scipy dependency)."""
    from math import erfc, sqrt
    return np.array([0.5 * erfc(abs(v) / sqrt(2.0)) for v in np.atleast_1d(z)])


def study_regression(o: pd.DataFrame, ref_continent: str = "North America") -> pd.DataFrame:
    """Study-level weighted least squares: logit(mean mapping) ~ continent + region.

    Fitted at the STUDY level (not sample level) to avoid pseudoreplication, since
    within-study samples share a mapping rate; weighted by n samples per study so
    better-sampled studies count more. Continent coefficients are the amplicon-adjusted
    effect vs the reference continent, on the logit scale. Implemented in numpy so the
    script has no statsmodels/scipy version dependency. Two-sided p-values use a normal
    approximation (n_studies is large relative to the number of parameters).
    """
    st = study_table(o).dropna(subset=["continent", "region"]).copy()
    st = st[st["n"] >= 3]  # require a few samples for a stable study mean
    y = logit(st["cr"].values)
    w = st["n"].to_numpy(dtype=float)

    cont_lvls = [c for c in st["continent"].unique() if c != ref_continent]
    reg_lvls = sorted(st["region"].unique())[1:]  # drop first region as baseline
    cols, names = [np.ones(len(st))], ["Intercept"]
    for c in cont_lvls:
        cols.append((st["continent"] == c).to_numpy(float)); names.append(f"continent[{c}]")
    for r in reg_lvls:
        cols.append((st["region"] == r).to_numpy(float)); names.append(f"region[{r}]")
    X = np.column_stack(cols)

    XtW = X.T * w
    XtWX = XtW @ X
    beta = np.linalg.solve(XtWX, XtW @ y)
    resid = y - X @ beta
    dof = len(st) - X.shape[1]
    sigma2 = float((w * resid**2).sum() / dof)
    cov = sigma2 * np.linalg.inv(XtWX)
    se = np.sqrt(np.diag(cov))
    z = beta / se
    p = 2.0 * _norm_sf(z)

    res = pd.DataFrame({
        "term": names, "coef_logit": beta,
        "ci_low": beta - 1.96 * se, "ci_high": beta + 1.96 * se,
        "p_value": p,
    })
    res["is_continent"] = res["term"].str.startswith("continent")
    res = pd.concat([res[res["is_continent"]], res[~res["is_continent"]]]).drop(columns="is_continent")
    res.attrs["n_studies"] = len(st)
    res.attrs["ref_continent"] = ref_continent
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta", type=Path, default=DEFAULT_META)
    ap.add_argument("--figdir", type=Path, default=DEFAULT_FIGDIR)
    ap.add_argument("--tabdir", type=Path, default=DEFAULT_TABDIR)
    ap.add_argument("--amplicon", default="V4",
                    help="amplicon region used for the held-constant control figure")
    args = ap.parse_args()

    o = load_oral(args.meta)
    args.tabdir.mkdir(parents=True, exist_ok=True)

    # 1 + 2: figures
    make_figure(
        o, args.figdir / "qc_cr_mapping_by_continent",
        "Closed-reference (GG2) mapping rate by continent (all amplicons)",
        "Dots: left = individual samples (subsampled to 400/continent); "
        "right = individual studies (size proportional to sqrt n). "
        "NB: not adjusted for amplicon region — see V4-only companion.",
    )
    o_amp = o[o["sequencing_region_cleaned"] == args.amplicon]
    make_figure(
        o_amp, args.figdir / f"qc_cr_mapping_by_continent_{args.amplicon}only",
        f"CR mapping by continent — {args.amplicon} amplicon only (region held constant)",
        f"Amplicon held constant ({args.amplicon}). Dots: left = samples "
        "(<=400/continent); right = studies (size proportional to sqrt n).",
    )

    # summary tables
    summ = pd.concat([summarize(o, "all_amplicons"),
                      summarize(o_amp, f"{args.amplicon}_only")], ignore_index=True)
    summ_path = args.tabdir / "cr_mapping_by_continent_summary.csv"
    summ.to_csv(summ_path, index=False)
    log(f"wrote {summ_path}")
    print(summ.to_string(index=False))

    # 3: amplicon-adjusted study-level regression
    reg = study_regression(o)
    reg_path = args.tabdir / "cr_mapping_continent_regression.csv"
    reg.to_csv(reg_path, index=False)
    log(f"wrote {reg_path}  (study-level WLS, n_studies={reg.attrs['n_studies']}, "
        f"ref={reg.attrs['ref_continent']})")
    print("\nAmplicon-adjusted continent effect (logit scale, vs North America):")
    print(reg[reg["term"].str.contains("continent")].to_string(index=False))


if __name__ == "__main__":
    main()

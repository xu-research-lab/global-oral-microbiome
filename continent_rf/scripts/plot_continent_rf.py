#!/usr/bin/env python3
"""Plot Saliva/Gut OVO AUC panels from separate or combined summary TSVs."""

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from rf_data import REGIONS

MODULE_DIR = Path(__file__).resolve().parents[1]
FIGURE_NAME = "region_ovo_oral_gut_two_panels_upper_triangle"


def _ovo_pair_frame(summary, site):
    frame = summary[summary["site"] == site]
    classes = [
        r for r in REGIONS if r in set(frame["region_1"]) | set(frame["region_2"])
    ]
    matrix = pd.DataFrame(np.nan, index=classes, columns=classes)
    for row in frame.itertuples(index=False):
        matrix.loc[row.region_1, row.region_2] = row.mean_auc
        matrix.loc[row.region_2, row.region_1] = row.mean_auc
    return matrix


def plot_summary(summary, outdir):
    sites = [s for s in ("Saliva", "Gut") if s in set(summary["site"])]
    if not sites:
        raise ValueError("Summary has no Saliva or Gut rows")
    summary = summary[summary["site"].isin(sites)].copy()
    if summary.duplicated(["site", "region_1", "region_2"]).any():
        raise ValueError("Duplicate site/region pairs in OVO summaries")
    if not np.isfinite(summary["mean_auc"]).all() or not summary["mean_auc"].between(0, 1).all():
        raise ValueError("mean_auc must contain finite values between 0 and 1")
    if not set(summary["region_1"]).union(summary["region_2"]).issubset(REGIONS):
        raise ValueError("Unknown regions in OVO summaries")
    vmin = np.floor(
        min(summary.loc[summary.site == s, "mean_auc"].min() for s in sites) * 20
    ) / 20
    vmin = min(vmin, 0.95)
    fig = plt.figure(figsize=(7.0, 3.3))
    panel = (2.3 / 7.0, 2.3 / 3.3)
    x0 = 0.10
    y0 = 0.22
    gap = 0.121
    cax = fig.add_axes((0.93, y0, 0.012, panel[1]))
    for k, site in enumerate(sites):
        ax = fig.add_axes((x0 + k * (panel[0] + gap), y0, panel[0], panel[1]))
        matrix = _ovo_pair_frame(summary, site)
        classes = [r for r in REGIONS if r in matrix.index]
        matrix = matrix.reindex(index=classes, columns=classes)
        matrix = matrix.iloc[1:, :-1]
        labels = matrix.map(lambda v: f"{v * 100:.1f}" if pd.notna(v) else "")
        mask = np.triu(np.ones(matrix.shape, dtype=bool), k=1)
        sns.heatmap(
            matrix,
            mask=mask,
            annot=labels,
            fmt="",
            annot_kws={"size": 8.5, "color": "white"},
            cmap="Blues",
            vmin=vmin,
            vmax=1,
            linewidths=0.8,
            linecolor="white",
            cbar=False,
            square=True,
            xticklabels=classes[:-1],
            yticklabels=classes[1:],
            ax=ax,
        )
        displayed = "Oral (saliva)" if site == "Saliva" else site
        ax.set_title(displayed, fontsize=11, pad=8)
        ax.tick_params(labelsize=7.5)
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", pad=1)
        for lbl in ax.get_xticklabels():
            lbl.set_rotation_mode("anchor")
            lbl.set_ha("right")
        norm_vals = (matrix.to_numpy() - vmin) / (1 - vmin)
        shown = norm_vals[(~mask) & matrix.notna().to_numpy()]
        for t, frac in zip(ax.texts, shown):
            t.set_color("#333333" if frac < 0.30 else "white")
    sm = plt.cm.ScalarMappable(cmap="Blues", norm=plt.Normalize(vmin, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, label="Mean AUC (%)")
    cbar.ax.yaxis.set_major_formatter(
        matplotlib.ticker.PercentFormatter(xmax=1.0, decimals=0, symbol="")
    )
    fig.savefig(outdir / f"{FIGURE_NAME}.png",
                dpi=300, bbox_inches="tight")
    fig.savefig(outdir / f"{FIGURE_NAME}.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"Plots written to {outdir}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=Path,
                        default=MODULE_DIR / os.environ.get("PLOT_RESULT_DIR", "result_genus_5rep"))
    parser.add_argument("--summaries", type=Path, nargs="+",
                        help="Explicit summary TSVs; also supports old combined Oral/Gut tables")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    paths = args.summaries or [args.result_dir / d / "region_ovo_summary.tsv" for d in ("oral", "gut")]
    summary = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    required = {"site", "region_1", "region_2", "mean_auc"}
    if not required.issubset(summary.columns):
        parser.error(f"Missing summary columns: {sorted(required - set(summary.columns))}")
    outdir = args.output_dir or args.result_dir
    outdir.mkdir(parents=True, exist_ok=True)
    plot_summary(summary, outdir)


if __name__ == "__main__":
    main()

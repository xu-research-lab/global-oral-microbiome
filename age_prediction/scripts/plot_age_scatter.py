#!/usr/bin/env python3
"""Plot pooled HC true-versus-predicted ages for project or random CV."""

import argparse
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import r2_score, mean_absolute_error

from age_common import FIGURE_DIR, ORAL_SITES, RESULT, SCHEMES, load_metadata, site_predictions

OUT_DIR = FIGURE_DIR


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scheme", choices=SCHEMES, default="random")
    args = parser.parse_args(argv)
    directory, token = SCHEMES[args.scheme]
    data_dir = RESULT / directory
    name = "random_split" if args.scheme == "random" else "cross_study"
    prefix = f"age_{name}_scatter_by_site"
    metadata = load_metadata()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # One panel per oral site.
    fig, axes = plt.subplots(1, 4, figsize=(16 * 2 / 3, 4.5 * 2 / 3))
    stat_rows = []

    for ax, site in zip(axes.ravel(), ORAL_SITES):
        site_df = site_predictions(site, metadata, data_dir, token)
        site_df = site_df[site_df.status_cleaned == "HC"]

        true_all = site_df.true_age.values
        pred_all = site_df.pred_age.values
        if len(true_all) == 0:
            ax.set_axis_off()
            ax.set_title(f"{site} — no data")
            continue

        mae = mean_absolute_error(true_all, pred_all)
        r2 = r2_score(true_all, pred_all)
        rho = spearmanr(true_all, pred_all)[0]
        stat_rows.append(dict(site=site, n=len(true_all), MAE=mae, R2=r2,
                              spearman=rho))

        ax.scatter(true_all, pred_all, alpha=0.30, s=16, color="steelblue",
                   edgecolor="none")

        lo = min(true_all.min(), pred_all.min()) - 2
        hi = max(true_all.max(), pred_all.max()) + 2
        ax.plot([lo, hi], [lo, hi], "k--", lw=1, alpha=0.5, label="Ideal (y = x)")

        z = np.polyfit(true_all, pred_all, 1)
        xl = np.linspace(true_all.min(), true_all.max(), 100)
        ax.plot(xl, np.poly1d(z)(xl), color="tomato", lw=2, label="Fit")

        ax.text(0.04, 0.96,
                f"MAE = {mae:.2f} yr\n$R^2$ = {r2:.3f}"
                f"\nSpearman ρ = {rho:.3f}\nn = {len(true_all)}",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.8))

        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xlabel("Chronological age (years)")
        ax.set_ylabel("Predicted age (years)")
        ax.set_title(site, fontsize=10)
        ax.legend(loc="lower right", frameon=False, fontsize=7)

    if stat_rows:
        stats = pd.DataFrame(stat_rows)
        stats.to_csv(OUT_DIR / f"{prefix}_stats.csv",
                     index=False)
        print("\n=== pooled HC stats (HC only) ===")
        print(stats.round(3).to_string(index=False))

    fig.tight_layout()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"{prefix}.{ext}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_DIR / prefix}.png")


if __name__ == "__main__":
    main()

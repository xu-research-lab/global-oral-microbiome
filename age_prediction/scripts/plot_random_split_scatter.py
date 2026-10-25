#!/usr/bin/env python
# =====================================================================
# Plot the random-split (sample-level, studies mixed) age-prediction
# scatter: one panel per oral site, held-out samples pooled over all CV
# folds (HC only, same rule as the existing cross-study figure).
# Style mirrors age_cross_study_scatter_by_site.png so the two versions
# are directly comparable:
#   blue scatter + grey dashed ideal line + red fit + metric box.
# Reads the per-fold CSVs from data/age_prediction/random_split_rf/
# (written by run_age_predict_random_split.py) and saves PNG + PDF to
#   age_prediction/result/figures/
#     age_random_split_scatter_by_site.{png,pdf}
# =====================================================================
from pathlib import Path

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import r2_score, mean_absolute_error

BASE = Path(os.environ.get("GLOBAL_ORAL_ROOT") or Path(__file__).resolve().parents[3])
RESULT = Path(__file__).resolve().parents[1] / "result"
DATA_DIR = RESULT / "random_split_rf"
METADATA_PATH = (BASE / "data" / "metadata-and-biom-v2-qc"
                 / "metadata_oral_qc_dedup.tsv")
OUT_DIR = RESULT / "figures"

ORAL_SITES = ["saliva", "tooth", "tongue", "oral mucosa"]


def main():
    metadata = pd.read_csv(METADATA_PATH, sep="\t", low_memory=False)
    metadata["sample_id"] = metadata["sample_id"].astype(str)

    # same 1x4 geometry as the notebook cell (16*2/3, 4.5*2/3) that made
    # the cross-study figure
    fig, axes = plt.subplots(1, 4, figsize=(16 * 2 / 3, 4.5 * 2 / 3))
    stat_rows = []

    for ax, site in zip(axes.ravel(), ORAL_SITES):
        site_tag = site.replace(" ", "_")
        # CV script writes the raw site name in the filename (may contain
        # a space); accept both spaced and underscored variants.
        files = sorted(DATA_DIR.glob(f"age_predict_test_*_{site_tag}_random_range_1_80.csv"))
        if not files:
            files = sorted(DATA_DIR.glob(
                f"age_predict_test_*_{site}_random_range_1_80.csv"))
        if not files:
            print(f"[{site}] no prediction CSVs found, skipping")
            ax.set_axis_off()
            ax.set_title(f"{site} — no data")
            continue
        site_df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        site_df = site_df.merge(metadata, on="sample_id", how="inner")
        site_df = site_df.loc[site_df.status_cleaned == "HC"]   # HC only, as before

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
        stats.to_csv(OUT_DIR / "age_random_split_scatter_by_site_stats.csv",
                     index=False)
        print("\n=== pooled HC stats (HC only) ===")
        print(stats.round(3).to_string(index=False))

    fig.tight_layout()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"age_random_split_scatter_by_site.{ext}",
                    dpi=300, bbox_inches="tight")
    print(f"Saved: {OUT_DIR / 'age_random_split_scatter_by_site.png'}")


if __name__ == "__main__":
    main()
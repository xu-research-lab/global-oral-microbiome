#!/usr/bin/env python3
"""Quick-look figure: per-site HC delta_age medians (+ IQR), LOSO vs random split.

Reads result/figures/hc_delta_age_site_summary.xlsx (built by
make_hc_delta_age_report.py); colors follow the age-acceleration figure
conventions (#2471a3 / #c0392b / grays). Writes
result/figures/hc_delta_age_median_site_summary.png.

Env used previously: /home/xuhuan/miniconda3/envs/R4/bin/python
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RESULT_ROOT = "/home/xuhuan/global_oral/age_prediction/result"
XLSX = os.path.join(RESULT_ROOT, "figures", "hc_delta_age_site_summary.xlsx")
OUT = os.path.join(RESULT_ROOT, "figures", "hc_delta_age_median_site_summary.png")

BLUE, NEG_RED, GRAY = "#2471a3", "#c0392b", "0.55"
SITE_ORDER = ["saliva", "tooth", "tongue", "oral mucosa"]
RNG = np.random.RandomState(42)


def median_colors(vals):
    """diverging ink: red above 0 (over-estimate), blue below (under-estimate)"""
    return [NEG_RED if v > 0 else BLUE for v in vals]


def main():
    site = pd.read_excel(XLSX, sheet_name="hc_site_summary_lso")
    site_rand = pd.read_excel(XLSX, sheet_name="hc_site_summary_random")
    samples = pd.read_excel(XLSX, sheet_name="hc_delta_samples")
    loso = samples[samples["scheme"] == "leave-studies-out"]

    site = site.set_index("site").reindex(SITE_ORDER).reset_index()
    site_rand = site_rand.set_index("site").reindex(SITE_ORDER).reset_index()

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2),
                             gridspec_kw={"width_ratios": [1.35, 1]})

    # ---- panel A: LOSO per-site sample strip + IQR bar + median tick -----------
    ax = axes[0]
    y = np.arange(len(SITE_ORDER))[::-1]
    for i, s in enumerate(SITE_ORDER):
        yi = y[i]
        v = loso.loc[loso["site"] == s, "delta_age"].values
        row = site[site["site"] == s].iloc[0]
        # light sample strip (max 400 jittered points) to show the spread
        show = v if len(v) <= 400 else RNG.choice(v, 400, replace=False)
        ax.plot(show, np.full(len(show), yi) + RNG.uniform(-.16, .16, len(show)),
                "o", ms=2, color="0.80", alpha=.45, zorder=1)
        ax.hlines(yi, row["q1"], row["q3"], color="0.3", lw=3.2, zorder=2,
                  label="_" if i else "HC IQR")
        med = row["median_delta (= figure hc_med anchor)"]
        col = median_colors([med])[0]
        ax.plot(med, yi, "|", ms=17, mew=3.2, color=col, zorder=3,
                label="_" if i else "HC median (LOSO)")
        ax.annotate(f"{med:+.2f}", (med, yi - .33), ha="center", fontsize=9,
                    color="0.3", zorder=4)
        ax.annotate(f"n={int(row['n'])}", (58, yi), va="center", ha="left",
                    fontsize=8, color=GRAY)
    ax.axvline(0, color="k", ls="--", lw=.9, zorder=0)
    ax.set_yticks(y, SITE_ORDER)
    ax.set_ylim(-1.0, len(SITE_ORDER) - 0.2)
    ax.set_xlim(-60, 68)
    ax.set_xlabel("HC Δage = predicted age − true age (years)")
    ax.set_title("Leave-studies-out, pooled HC per site", fontsize=10.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="0.92", zorder=0)

    # ---- panel B: LOSO vs random-split medians ---------------------------------
    ax = axes[1]
    w, gap = .3, .06
    med_l = site["median_delta (= figure hc_med anchor)"].values
    med_r = site_rand["median_delta (= figure hc_med anchor)"].values
    xl = y - (w / 2 + gap / 2)
    xr = y + (w / 2 + gap / 2)
    ax.axvline(0, color="k", ls="--", lw=.9, zorder=0)
    b1 = ax.barh(xl, med_l, height=w, color=BLUE, label="leave-studies-out")
    b2 = ax.barh(xr, med_r, height=w, color="0.8", edgecolor="0.4", lw=.5,
                 label="random split")
    for bars, vals in ((b1, med_l), (b2, med_r)):
        for b, v in zip(bars, vals):
            ax.annotate(f"{v:+.1f}", (v + (1.5 if v >= 0 else -1.5), b.get_y() + w / 2),
                        va="center", ha="left" if v >= 0 else "right",
                        fontsize=8, color="0.3")
    ax.set_yticks(y, SITE_ORDER)
    ax.set_ylim(-1.0, len(SITE_ORDER) - 0.2)
    ax.set_xlim(-12, 17)
    ax.set_xlabel("HC Δage median (years)")
    ax.set_title("Median: honest vs optimistic CV", fontsize=10.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="0.92", zorder=0)
    ax.legend(loc="lower left", frameon=False, fontsize=8.5)

    fig.suptitle("HC Δage per site (out-of-fold predictions, HC test samples only)",
                 fontsize=11.5, y=.995)
    fig.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"[done] wrote {OUT}")


if __name__ == "__main__":
    main()
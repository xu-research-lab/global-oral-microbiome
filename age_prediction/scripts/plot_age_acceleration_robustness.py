#!/usr/bin/env python
# =====================================================================
# Supplementary dumbbell plot for the disease excess-age-gap figure:
# for every site x status contrast, show the unadjusted median gap
# difference (as in plot_age_acceleration_1row.py) side by side with the
# age-standardized and within-study estimates from
# robustness_age_gap_adjustment.py.  Rows are ordered like the main
# figure (ascending unadjusted delta); a verdict column flags contrasts
# whose effect size was driven by age-composition or study-level
# confounding.
# Saves: result/figures/age_gap_robustness_dumbbell.{png,pdf}
# =====================================================================
from pathlib import Path

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RESULT = Path(__file__).resolve().parents[1] / "result"
OUT_DIR = RESULT / "figures"

MAIN_STATS = OUT_DIR / "age_disease_microbiome_age_acceleration_1row_stats.csv"
ROB_STATS = OUT_DIR / "age_gap_robustness_stats.csv"

ORAL_SITES = ["saliva", "tooth", "tongue", "oral mucosa"]

# same abbreviations as plot_age_acceleration_1row.py
DISEASE_ABBR = {
    "Aphthous": "AS", "Atopic Dermatitis": "AtD",
    "Axial Spondyloarthritis": "axSpA", "Chagas Disease": "CD",
    "Celiac Disease": "Celiac Disease", "Colorectal Cancer": "CRC",
    "COVID-19": "COVID-19", "Dental Caries": "DC", "Depression": "Depression",
    "Gastroesophageal Reflux Disease": "GERD", "Head and Neck Cancer": "HNC",
    "HIV": "HIV", "IBD": "IBD", "ICU": "ICU", "Lung Cancer": "LC",
    "Obesity": "Obesity", "Oral Lichen Planus": "OLP", "Oral Mucositis": "OM",
    "Oral Squamous Cell Carcinoma": "OSCC", "Oral Submucous Fibrosis": "OSMF",
    "Oral Verrucous Hyperplasia": "OVH", "OSCC": "OSCC",
    "Periodontal Disease": "PD",
}

C_UNADJ, C_AGE, C_STUDY = "#404040", "#2471a3", "#1e8449"
XMIN, XMAX = -42, 42


def classify(r):
    """One-word verdict for a contrast, from the adjusted estimates."""
    was_sig = r.fdr < .05 if np.isfinite(r.get("fdr", np.nan)) else False
    astd_sig = (np.isfinite(r.astd_lo) and np.isfinite(r.astd_hi)
                and r.astd_lo * r.astd_hi > 0)
    flipped = astd_sig and np.sign(r.astd_delta) != np.sign(r.unadj_delta)
    shrunken = (np.isfinite(r.unadj_delta) and abs(r.unadj_delta) >= 8
                and ((np.isfinite(r.instudy_delta)
                      and abs(r.instudy_delta) < abs(r.unadj_delta) / 3)
                     or (astd_sig and abs(r.astd_delta) < abs(r.unadj_delta) / 3)))
    if was_sig and flipped:
        return "sign flip"
    if was_sig and shrunken:
        return "age/study artifact"
    if was_sig and not astd_sig and np.isfinite(r.astd_delta):
        return "ns after age adj."
    if was_sig and astd_sig and np.sign(r.astd_delta) == np.sign(r.unadj_delta) \
            and abs(r.astd_delta) >= 3:
        if np.isfinite(r.instudy_delta) and (
                np.sign(r.instudy_delta) != np.sign(r.unadj_delta)
                or abs(r.instudy_delta) < 3):
            return "study-inconsistent"
        return "robust"
    return ""


def main():
    main_df = pd.read_csv(MAIN_STATS)
    rob = pd.read_csv(ROB_STATS)
    df = rob.merge(main_df[["site", "status", "fdr", "hc_n"]],
                   on=["site", "status"], how="left")
    df["verdict"] = df.apply(classify, axis=1)

    n_rows = [len(df[df.site == j]) for j in ORAL_SITES]
    fig = plt.figure(figsize=(11, 0.34 * sum(n_rows) + 2.6))
    gs = fig.add_gridspec(4, 1, height_ratios=n_rows, hspace=0.08,
                          left=0.17, right=0.80, top=0.925, bottom=0.05)

    axes = [fig.add_subplot(gs[k]) for k in range(4)]
    for ax, j, nr in zip(axes, ORAL_SITES, n_rows):
        sub = df[df.site == j].sort_values("unadj_delta").reset_index(drop=True)
        y = np.arange(len(sub))[::-1]

        # connector line through the available estimates of each row
        for yi, (_, r) in zip(y, sub.iterrows()):
            xs = [r.unadj_delta]
            if np.isfinite(r.instudy_delta):
                xs += [r.astd_delta, r.instudy_delta]
            else:
                xs += [r.astd_delta]
            ax.plot(xs, [yi] * len(xs), color="0.85", lw=1.2, zorder=1)

        hc_ref = sub.hc_med_gap.iloc[0] if len(sub) else np.nan
        ax.axvline(0, color="k", lw=.8, zorder=0)

        # age-standardized CI
        sig = sub.astd_lo * sub.astd_hi > 0
        ax.errorbar(sub.astd_delta, y,
                    xerr=[sub.astd_delta - sub.astd_lo, sub.astd_hi - sub.astd_delta],
                    fmt="none", ecolor=C_AGE, elinewidth=1.4, capsize=2.2, zorder=2)
        ax.scatter(sub.unadj_delta, y, s=42, c=C_UNADJ, zorder=3,
                   label="Unadjusted")
        ax.scatter(sub.astd_delta, y, s=46, marker="s", zorder=3,
                   facecolors=np.where(sig, C_AGE, "white"),
                   edgecolors=C_AGE, linewidths=1.2)
        ins = np.isfinite(sub.instudy_delta).values
        ax.scatter(sub.instudy_delta.values[ins], y[ins], s=46, marker="D",
                   c=C_STUDY, zorder=3)

        ax.set_yticks(y)
        ax.set_yticklabels([DISEASE_ABBR.get(s, s) for s in sub.status],
                           fontsize=10)
        ax.set_xlim(XMIN, XMAX)
        ax.set_ylim(-0.6, len(sub) - 0.4)
        ax.set_title(f"{j}  (HC n = {int(sub.hc_n.iloc[0])}, "
                     f"HC med gap = {hc_ref:+.1f} yr)" if len(sub) else j,
                     fontsize=12, loc="left")
        ax.text(XMAX + 1, len(sub) - 0.4 + .2, "verdict", fontsize=10,
                style="italic", va="bottom", ha="left", fontweight="bold")
        for yi, r in zip(y, sub.itertuples()):
            if r.verdict:
                ax.text(XMAX + 1, yi, r.verdict, fontsize=9, color="0.25",
                        va="center", ha="left")
        if ax is not fig.axes[-1]:
            ax.tick_params(labelbottom=False)

    fig.axes[-1].set_xlabel("Median age-gap difference vs HC (years)",
                            fontsize=11)

    handles = [
        Line2D([], [], marker="o", ls="", color=C_UNADJ,
               label="Unadjusted (main figure)"),
        Line2D([], [], marker="s", ls="", color=C_AGE,
               label="Age-standardized\n(fill = 95% CI excl. 0)"),
        Line2D([], [], marker="D", ls="", color=C_STUDY,
               label="Within-study HC reference\n(batch-corrected; NaN if no in-study HC)"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 0.998),
               ncol=3, fontsize=9, frameon=False)

    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"age_gap_robustness_dumbbell.{ext}", dpi=300)
    print(f"Saved: {OUT_DIR / 'age_gap_robustness_dumbbell.png'}")


if __name__ == "__main__":
    main()
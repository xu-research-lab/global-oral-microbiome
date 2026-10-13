#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_oral_fraction_by_disease.py

Disease association of the gut's oral-derived (paired-oral) fraction
(SourceTracker2): forest of LME coefficients, sourcetracker_oral_fraction ~
status_cleaned (reference = HC) with random intercepts, from the
single_pair model run. Error bars are 95% CIs (beta +/- 1.96 SE); per-
status n is read from the model input table and shown next to each label.

Input : results/single_pair/06_models/lme_all_coefficients.tsv
        results/single_pair/06_models/lme_input_per_pair_table.tsv
Output: results/new/figures/03_oral_fraction/
          oral_fraction_status_lme.tsv
          oral_fraction_by_disease_lme.png / .pdf
"""

from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
# single source of truth for disease abbreviations (user-supplied mapping;
# shared with figs 07/08/09/10)
sys.path.insert(0, str(ROOT / "scripts" / "python"))
from oral_gut_pair.report_figures import DISEASE_LABEL_SHORT as DISEASE_ABBREV  # noqa: E402
IN_COEF = ROOT / "results/single_pair/06_models/lme_all_coefficients.tsv"
IN_PERPAIR = ROOT / "results/single_pair/06_models/lme_input_per_pair_table.tsv"
OUT_DIR = ROOT / "results/new/figures/03_oral_fraction"

TERM_RE = re.compile(r'status_cleaned, Treatment\(reference="HC"\)\)\[T\.(.*)\]')

SIG_COLORS = {"***": "#C44E52", "**": "#DD8452", "*": "#E7C548", "ns": "#9CA3AF"}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    coef = pd.read_csv(IN_COEF, sep="\t", dtype=str, low_memory=False)
    coef = coef.loc[coef["metric_folder"].eq("sourcetracker_oral_fraction")].copy()
    coef["is_random_effect"] = coef["is_random_effect"].astype(str).str.lower()

    sub = coef.loc[
        coef["model"].eq("model3_status")
        & coef["is_random_effect"].eq("false")
        & coef["term_raw"].str.contains("status_cleaned", na=False)
    ].copy()
    sub["disease"] = sub["term_raw"].str.extract(TERM_RE)
    for c in ("estimate", "std_error", "ci95_low", "ci95_high", "p_value"):
        sub[c] = pd.to_numeric(sub[c], errors="coerce")
    sub["q_value_global"] = pd.to_numeric(sub["q_value_global"], errors="coerce")
    sub["sig_label"] = sub["p_significance"].fillna("ns")

    res = sub[["disease", "estimate", "std_error", "ci95_low", "ci95_high",
               "p_value", "q_value_global", "p_significance",
               "q_significance_global"]].sort_values("estimate")
    tsv_fp = OUT_DIR / "oral_fraction_status_lme.tsv"
    res.to_csv(tsv_fp, sep="\t", index=False)

    # ===== forest plot =====
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    res = res.reset_index(drop=True)
    y = np.arange(len(res))[::-1]

    # vertical spacing: per-row height factor tightened to 0.19 in (user
    # request: even more compact vertically) — at this spacing the 9 pt
    # y labels would collide, so y ticks render at 8.5 pt (user-approved
    # trade-off); x axis doubled earlier (user request) so the upper-right
    # region is empty enough to host the legend inside the axes
    fig, ax = plt.subplots(figsize=(4.2, 0.19 * len(res) + 1.1))
    for (_, r), yi in zip(res.iterrows(), y):
        col = SIG_COLORS.get(r["p_significance"], "grey")
        ax.errorbar(r["estimate"], yi,
                    xerr=[[r["estimate"] - r["ci95_low"]],
                          [r["ci95_high"] - r["estimate"]]],
                    fmt="o", ms=5, color=col, ecolor=col, elinewidth=1.2,
                    capsize=2, zorder=3)
    ax.axvline(0, color="grey", lw=0.9, ls="--", zorder=1)

    # per-status n from the model input table (statuses the model actually
    # fitted — rows with a status but non-missing paired_oral_fraction)
    pair = pd.read_csv(IN_PERPAIR, sep="\t", dtype=str, low_memory=False)
    vals = pd.to_numeric(pair["paired_oral_fraction"], errors="coerce")
    pair = pair.loc[vals.notna(), "status_cleaned"]
    n_by_status = pair.dropna().value_counts().to_dict()

    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{DISEASE_ABBREV.get(d, d)}  (n={n_by_status.get(d, 0)})"
         for d in res["disease"]],
        fontsize=8.5,
    )
    ax.set_xlabel("LME β (disease vs HC)", fontsize=8)

    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=k)
               for k, c in SIG_COLORS.items()]
    # legend to the right, inside the axes (top rows' CIs stay left of it)
    ax.legend(handles=handles, title="Significance", frameon=False,
              fontsize=7.5, title_fontsize=7.5, loc="upper right",
              labelspacing=0.4, handletextpad=0.4)

    # fixed ticks at 0 / 0.05 / 0.10 (user request)
    ax.set_xticks([0.0, 0.05, 0.10])
    ax.set_xticklabels(["0.00", "0.05", "0.10"], fontsize=8)
    # title / axis font sizes matched to fig 06 (alpha pair-vs-random):
    # title 10 pt, x-axis (ticks + label) 8 pt, y ticks 8.5 pt.
    # Title centered over the WHOLE figure (suptitle, not ax.title) so the
    # long title text stops driving the tight-bbox width — without this the
    # axes sit right of the long y labels and a title centered on the axes
    # leaves blank canvas on the right
    fig.suptitle("Oral-derived fraction in the gut vs. disease status",
                 fontsize=10)
    fig.tight_layout()

    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"oral_fraction_by_disease_lme.{ext}", dpi=300,
                    bbox_inches="tight")
    print(f"done: {OUT_DIR}")
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
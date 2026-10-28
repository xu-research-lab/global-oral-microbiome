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

# unified plot-area geometry (user request 2026-10-20): identical vertical
# rect to figs 05/06/07/11 — canvas 3.40 in tall, bottom edge 0.70 in /
# top edge 2.96 in from the canvas bottom — saved WITHOUT
# bbox_inches="tight". Y labels moved to the RIGHT (user request
# 2026-10-20), and (follow-up request) the boxed plot area is now the SAME
# BOX as fig 07's square panel: axes at 0.62 in from the canvas left, box
# 2.26 x 2.26 in; only the right-hand label margin makes FIG_W 4.75
# (longest, "COVID-19 discharge  (n=18)", 1.68 in at 8.5 pt)
CANVAS_H = 3.40
AX_BOTTOM_IN = 0.70
AX_TOP_IN = 2.96
FIG_W = 4.75
AX_RECT_F = [0.62 / FIG_W, AX_BOTTOM_IN / CANVAS_H,
             2.26 / FIG_W, (AX_TOP_IN - AX_BOTTOM_IN) / CANVAS_H]


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
    # main figure colours by RAW p (user request 2026-10-22; q-based colours
    # made ICU ns, user reverted to raw p); the BH-q-adjusted version is an
    # ADDITIONAL file (oral_fraction_by_disease_lme_qadjusted.{png,pdf},
    # run the BH-q variant as a SEPARATE file, never overwriting the raw-p one) so the two stay in step
    sub["sig_label"] = sub["p_significance"].fillna("ns")

    res = sub[["disease", "estimate", "std_error", "ci95_low", "ci95_high",
               "p_value", "q_value_global", "p_significance",
               "q_significance_global"]].sort_values("estimate")
    tsv_fp = OUT_DIR / "oral_fraction_status_lme.tsv"
    res.to_csv(tsv_fp, sep="\t", index=False)

    # per-status n from the model input table (statuses the model actually
    # fitted — rows with a status but non-missing paired_oral_fraction)
    pair = pd.read_csv(IN_PERPAIR, sep="\t", dtype=str, low_memory=False)
    vals = pd.to_numeric(pair["paired_oral_fraction"], errors="coerce")
    pair = pair.loc[vals.notna(), "status_cleaned"]
    n_by_status = pair.dropna().value_counts().to_dict()

    # ===== forest plot =====
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def draw_forest(sig_col: str, stem: str) -> None:
        """One forest figure; colours follow sig_col ('p_significance' |
        'q_significance_global'); stem appended to the output filename."""
        res2 = res.reset_index(drop=True)
        y = np.arange(len(res2))[::-1]

        # vertical spacing: the 15 coefficient rows must fit the SHARED axes
        # height (2.26 in across 03/05/06/07/11 — user request) — 15 rows over
        # the ~9 pt font needs the full box; explicit inches-based layout,
        # saved WITHOUT bbox_inches="tight" so the axes rect is exact
        fig = plt.figure(figsize=(FIG_W, CANVAS_H))
        ax = fig.add_axes(AX_RECT_F)
        for (_, r), yi in zip(res2.iterrows(), y):
            col = SIG_COLORS.get(r[sig_col], "grey")
            ax.errorbar(r["estimate"], yi,
                        xerr=[[r["estimate"] - r["ci95_low"]],
                              [r["ci95_high"] - r["estimate"]]],
                        fmt="o", ms=5, color=col, ecolor=col, elinewidth=1.2,
                        capsize=2, zorder=3)
        ax.axvline(0, color="grey", lw=0.9, ls="--", zorder=1)

        ax.set_yticks(y)
        ax.set_yticklabels(
            [f"{DISEASE_ABBREV.get(d, d)}  (n={n_by_status.get(d, 0)})"
             for d in res2["disease"]],
            fontsize=8.5,
        )
        # disease labels on the right edge (user request 2026-10-20) so the
        # forest box matches fig 07's panel for vertical stacking
        ax.yaxis.tick_right()
        ax.tick_params(axis="y", labelsize=8.5)
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
        # no tight_layout / bbox_inches="tight" — the axes rect is fixed (see
        # AX_RECT_F, unified geometry across figure families)

        for ext in ("png", "pdf"):
            fig.savefig(OUT_DIR / f"oral_fraction_by_disease_lme{stem}.{ext}",
                        dpi=300)
        plt.close(fig)

    draw_forest("p_significance", "")          # main: unadjusted p stars
    draw_forest("q_significance_global", "_qadjusted")  # extra: BH-q stars
    print(f"done: {OUT_DIR}")
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
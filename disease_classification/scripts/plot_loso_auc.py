"""Draw prevalence-filtered LOSO performance from per-project AUC tables."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from disease_common import MODULE_DIR, PLOT_SITES, RESULT_DIR

SITE_COLORS = {"saliva": "#2166AC", "tooth": "#B2182B"}


def plot_auc(output_dir=RESULT_DIR):
    labels = pd.read_csv(MODULE_DIR / "data/disease_labels.csv").set_index("disease").abbreviation.to_dict()
    data = pd.concat([
        pd.read_csv(RESULT_DIR / f"rf_loso_auc_all_diseases_genus_{site}.csv").assign(site=site)
        for site in PLOT_SITES
    ], ignore_index=True).dropna(subset=["auc"])
    if data.empty:
        raise ValueError("No valid AUC values to plot")
    order = (data.groupby("disease").agg(n_sites=("site", "nunique"), med=("auc", "median"))
             .sort_values(["n_sites", "med"], ascending=[False, False]).index)
    sites = [site for site in PLOT_SITES if site in set(data.site)]
    present = {d: [s for s in sites if ((data.disease == d) & (data.site == s)).any()]
               for d in order}
    site_slot, group_gap = 0.75, 0.40
    fig, ax = plt.subplots(figsize=(5.5, 2.45))
    ticks, cursor = [], 0.0
    for disease in order:
        width = len(present[disease]) * site_slot
        for j, site in enumerate(present[disease]):
            vals = data.loc[(data.disease == disease) & (data.site == site), "auc"]
            x = cursor + (j + 0.5) * site_slot
            ax.scatter(np.full(len(vals), x), vals.values, marker="o", s=44,
                       facecolors=SITE_COLORS[site], edgecolors=SITE_COLORS[site],
                       linewidths=1.1, alpha=0.5, zorder=3)
            ax.hlines(vals.mean(), x - 0.10, x + 0.10, colors=SITE_COLORS[site],
                      linewidths=2.6, zorder=4)
        ticks.append(cursor + width / 2)
        cursor += width + group_gap
    ax.axhline(0.5, ls="--", c="grey", lw=1, zorder=2)
    ax.set_ylabel("LOSO AUC")
    ax.set_ylim(0, 1.06)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xlim(-group_gap / 2, cursor - group_gap / 2)
    ax.set_xticks(ticks)
    ax.set_xticklabels([labels.get(d, d) for d in order], rotation=35,
                       ha="right", rotation_mode="anchor",
                       fontsize=float(np.clip(140.0 / len(order), 8, 12)))
    ax.spines[["top", "right"]].set_visible(False)
    handles = [Patch(facecolor=SITE_COLORS[s], edgecolor="none", label=s.title()) for s in sites]
    ax.legend(handles=handles, loc="upper center", ncol=len(sites), fontsize=9,
              frameon=True, framealpha=0.9, edgecolor="#DDDDDD", bbox_to_anchor=(0.5, -0.38))
    fig.tight_layout()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / f"loso_rf_1row_genus.{suffix}", dpi=200, bbox_inches="tight")
    plt.close(fig)

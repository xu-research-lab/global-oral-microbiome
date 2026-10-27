#!/usr/bin/env python3
"""
pipeline_v2 Step 14: supported-phylotype catalog (BEFORE HOMD novelty_class).

Definition (user, 2026-08-18): phylotypes passing the support gate
    max_subject_prevalence > 50 & n_bioprojects_total >= 2 & n_continents >= 2
AND clean_phylotype (no contamination-flagged member) — the SUPPORTED CATALOG
(cross-CONTINENT evidence). The HEADLINE layer = supported catalog + GG2>=99% +
HOMD-uncaptured.

This is the "strongly supported oral taxa landscape" result layer; the 12 STRICT
candidates are the clean + GG2>=99% + HOMD-uncaptured subset of this catalog.

Outputs (pipeline_v2/):
  - 14_supported_phylotypes.tsv   : the supported catalog (one row per phylotype)
  - 14_supported_summary.json
  - analysis_report/figs_v2/fig_v2_supported.png : composition figure
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
FIG = PIPE_V2 / "analysis_report" / "figs_v2"

SUPPORT = {"subject": 50, "bioprojects": 2, "continents": 2, "subject_gt": True,
           "bioprojects_ge": True, "continents_ge": True}
HL_PIDS = set(pd.read_csv(PIPE_V2 / "11_headline_candidates.tsv", sep="\t",
                          usecols=["phylotype_id"])["phylotype_id"])

BLUE = "#4E79A7"; RED = "#E15759"; GRAY = "#BAB0AC"; TEAL = "#76B7B2"; INK = "#2B2B2B"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.edgecolor": "#BBBBBB", "axes.linewidth": 0.8,
                     "axes.labelcolor": INK, "text.color": INK, "figure.dpi": 150})
log = lambda s: print(f"[14_supported] {s}", flush=True)


def main():
    nov = pd.read_csv(PIPE_V2 / "11_phylotypes_homd_novelty.tsv", sep="\t", low_memory=False)
    for c in ["rep_identity", "max_subject_prevalence", "n_bioprojects_total", "n_regions"]:
        nov[c] = pd.to_numeric(nov[c], errors="coerce")
    nov["clean_phylotype"] = nov["clean_phylotype"].astype(str).str.lower().isin(["true", "1"])
    nov["oral_clade"] = nov["oral_clade"].astype(str).str.lower().isin(["true", "1"])
    nov["novelty_class"] = nov["novelty_class"].fillna("homd_uncaptured")

    m = ((nov["max_subject_prevalence"] > 50)
         & (nov["n_bioprojects_total"] >= 2)
         & (nov["n_continents"] >= 2)
         & nov["clean_phylotype"])
    sup = nov[m].copy()
    sup["headline"] = sup["phylotype_id"].isin(HL_PIDS)
    sup["gg2_tier"] = np.where(sup["rep_identity"] >= 99, ">=99% (species)",
                      np.where(sup["rep_identity"] >= 97, "97-99% (SGB)", "<97% (unresolved)"))
    sup = sup.sort_values(["max_subject_prevalence", "n_bioprojects_total", "n_regions"],
                          ascending=False)
    sup.to_csv(PIPE_V2 / "14_supported_phylotypes.tsv", sep="\t", index=False)
    log(f"saved 14_supported_phylotypes.tsv ({len(sup):,} rows)")

    # composition figure: GG2 tier x HOMD novelty, annotated with clean/strict
    ctab = pd.crosstab(sup["gg2_tier"], sup["novelty_class"])
    order = ["<97% (unresolved)", "97-99% (SGB)", ">=99% (species)"]
    ctab = ctab.reindex([c for c in order if c in ctab.index])

    fig, ax = plt.subplots(figsize=(8, 4.4))
    y = np.arange(len(ctab))
    left = np.zeros(len(ctab))
    for cls, color in [("homd_captured", GRAY), ("homd_uncaptured", RED)]:
        if cls not in ctab.columns:
            continue
        vals = ctab[cls].to_numpy()
        ax.barh(y, vals, left=left, color=color, height=0.58,
                label="HOMD-captured" if cls == "homd_captured" else "HOMD-uncaptured")
        left += vals
    for yi, (idx, row) in enumerate(ctab.iterrows()):
        total = int(row.sum())
        ax.text(total + 2, yi, f"{total}", va="center", fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels(ctab.index, fontsize=9)
    ax.set_xlim(0, ctab.sum(axis=1).max() * 1.15)
    ax.set_xlabel("Supported phylotypes")
    ax.set_title(f"{len(sup)} clean supported phylotypes (subjects >50, BioProjects ≥2, continents ≥2)\n"
                 "by GG2 identity and HOMD status — headline highlighted",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.7)
    fig.tight_layout()
    fig.savefig(FIG / "fig_v2_supported.png", bbox_inches="tight")
    plt.close(fig)
    log("saved fig_v2_supported.png")

    summary = {
        "support_cutoffs": SUPPORT,
        "n_supported_clean": int(len(sup)),
        "by_gg2_tier": sup["gg2_tier"].value_counts().to_dict(),
        "by_novelty": sup["novelty_class"].value_counts().to_dict(),
        "clean": int(sup["clean_phylotype"].sum()),
        "oral_clade": int(sup["oral_clade"].sum()),
        "n_headline_among_supported": int(sup["headline"].sum()),
    }
    with open(PIPE_V2 / "14_supported_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log(f"summary: {summary}")
    log("Done!")


if __name__ == "__main__":
    main()

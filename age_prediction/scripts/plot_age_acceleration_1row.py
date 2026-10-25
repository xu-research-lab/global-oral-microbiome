#!/usr/bin/env python
# =====================================================================
# Replot "Microbiome age acceleration relative to healthy controls"
# (notebook age_predict.ipynb cells 39-40) from the per-fold prediction
# CSVs produced by run_age_acceleration_cv.py, with the four oral sites
# laid out in ONE ROW (1x4) instead of the original 2x2 grid.
# Saves PNG + PDF to
#   age_prediction/result/figures/
#     age_disease_microbiome_age_acceleration.{png,pdf}
# =====================================================================
from pathlib import Path

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

BASE = Path(os.environ.get("GLOBAL_ORAL_ROOT") or Path(__file__).resolve().parents[2])
RESULT = Path(__file__).resolve().parents[1] / "result"
DATA_DIR = RESULT / "leave_studies_out_rf"
METADATA_PATH = (BASE / "data" / "metadata-and-biom-v2-qc"
                 / "metadata_oral_qc_dedup.tsv")
OUT_DIR = RESULT / "figures"

ORAL_SITES = ["saliva", "tooth", "tongue", "oral mucosa"]
MIN_N = 10
N_BOOT = 2000

# Disease status -> abbreviation shown on the y axis (fallback: full name)
DISEASE_ABBR = {
    "Aphthous":                          "AS",
    "Atopic Dermatitis":                 "AtD",
    "Axial Spondyloarthritis":           "axSpA",
    "Chagas Disease":                    "CD",
    "Celiac Disease":                    "Celiac Disease",
    "Colorectal Cancer":                 "CRC",
    "COVID-19":                          "COVID-19",
    "Dental Caries":                     "DC",
    "Depression":                        "Depression",
    "Gastroesophageal Reflux Disease":   "GERD",
    "Head and Neck Cancer":              "HNC",
    "HIV":                               "HIV",
    "IBD":                               "IBD",
    "ICU":                               "ICU",
    "Lung Cancer":                       "LC",
    "Obesity":                           "Obesity",
    "Oral Lichen Planus":                "OLP",
    "Oral Mucositis":                    "OM",
    "Oral Squamous Cell Carcinoma":      "OSCC",
    "Oral Submucous Fibrosis":           "OSMF",
    "Oral Verrucous Hyperplasia":        "OVH",
    "OSCC":                              "OSCC",
    "Periodontal Disease":               "PD",
}


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = sum(np.sum(x > b) for x in a)
    lt = sum(np.sum(x < b) for x in a)
    return (gt - lt) / (len(a) * len(b))


def boot_ci_median_diff(dis, hc, rng, n_boot=N_BOOT):
    d = np.array([np.median(rng.choice(dis, len(dis), True)) -
                  np.median(rng.choice(hc, len(hc), True)) for _ in range(n_boot)])
    return np.percentile(d, [2.5, 97.5])


def stars(q):
    return "***" if q < .001 else "**" if q < .01 else "*" if q < .05 else "ns"


def main():
    metadata = pd.read_csv(METADATA_PATH, sep="\t", low_memory=False)
    metadata["sample_id"] = metadata["sample_id"].astype(str)

    long_rows, stat_rows = [], []
    rng = np.random.default_rng(0)

    for j in ORAL_SITES:
        site_tag = j.replace(" ", "_")
        # CV script writes the raw site name in the filename (may contain
        # a space, e.g. "oral mucosa"); accept both spaced and underscored.
        files = sorted(DATA_DIR.glob(f"age_predict_test_*_{site_tag}_genus_range_1_80.csv"))
        if not files:
            files = sorted(DATA_DIR.glob(
                f"age_predict_test_*_{j}_genus_range_1_80.csv"))
        if not files:
            print(f"[{j}] no prediction CSVs found, skipping")
            continue
        preds = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        preds["sample_id"] = preds["sample_id"].astype(str)

        meta_pick = metadata.loc[metadata.site_cleaned == j].set_index("sample_id")
        sid = np.intersect1d(meta_pick.index.values, preds.sample_id.values)
        if len(sid) == 0:
            print(f"[{j}] no sample overlap with metadata, skipping")
            continue

        df = preds.set_index("sample_id").loc[sid].copy()
        df["status"] = meta_pick.loc[sid, "status_cleaned"].values
        df["delta_age"] = df.pred_age.values - df.true_age.values
        df["site"] = j

        counts = df["status"].value_counts()
        keep = [s for s in counts.index if s == "HC" or counts[s] >= MIN_N]
        dropped = {s: int(counts[s]) for s in counts.index if s not in keep}
        if dropped:
            print(f"[{j}] delect n<{MIN_N} : {dropped}")
        df = df[df["status"].isin(keep)]
        long_rows.append(df)

        hc = df.loc[df.status == "HC", "delta_age"].values
        recs = []
        for s in [x for x in keep if x != "HC"]:
            v = df.loc[df.status == s, "delta_age"].values
            p = mannwhitneyu(v, hc, alternative="two-sided")[1]
            lo, hi = boot_ci_median_diff(v, hc, rng)
            recs.append(dict(site=j, status=s, n=len(v),
                             med=np.median(v),
                             delta=np.median(v) - np.median(hc),
                             ci_lo=lo, ci_hi=hi,
                             cliff=cliffs_delta(v, hc), p=p))
        if recs:
            r = pd.DataFrame(recs)
            r["fdr"] = bh_fdr(r["p"].values)
            r["hc_med"] = np.median(hc)
            r["hc_n"] = len(hc)
            stat_rows.append(r)

    long_df = pd.concat(long_rows, ignore_index=True)
    stats = pd.concat(stat_rows, ignore_index=True)
    stats["sig"] = stats["fdr"].map(stars)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stats.to_csv(OUT_DIR / "age_disease_microbiome_age_acceleration_1row_stats.csv", index=False)

    # ---------------- figure: saliva full left column; right column = 3 rows
    # tooth / tongue / oral mucosa. GridSpec with unequal width/height. -----
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(12, 4.3))
    gs = GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.75,
                  left=0.16, right=0.82, top=0.96, bottom=0.08,
                  width_ratios=[1.0, 1.0])

    def draw_panel(ax, j):
        s = stats[stats.site == j].sort_values("delta")
        if len(s) == 0:
            ax.set_axis_off()
            ax.set_title(f"{j} — no data")
            return
        ax.axvline(0, color="k", ls="--", lw=0.9, zorder=1)

        y = np.arange(len(s))
        col = np.where(s.fdr < .05, np.where(s.delta > 0, "#c0392b", "#2471a3"), "0.55")
        ax.hlines(y, s.ci_lo, s.ci_hi, color=col, lw=1.6, zorder=2)
        ax.scatter(s.delta, y, s=40 + 3 * np.sqrt(s.n), c=col,
                   edgecolor="k", lw=.5, zorder=3)

        for yi, row in zip(y, s.itertuples()):
            # Both labels sit OUTSIDE the axes (x>1.0 in axes coords) so they
            # never overlap the data points. n= is left-aligned at a fixed x,
            # and the sig stars are left-aligned at a fixed x just to the right
            # — so within a figure column every n= starts at the same place
            # and every star column starts at the same place, with a small gap.
            ax.text(1.02, yi, f"n={row.n}",
                    transform=ax.get_yaxis_transform(),
                    va="center", ha="left", fontsize=11, color="0.3")
            ax.text(1.25, yi, f"{row.sig}",
                    transform=ax.get_yaxis_transform(),
                    va="center", ha="left", fontsize=11, color="0.3")

        ax.set_yticks(y)
        ax.set_yticklabels([DISEASE_ABBR.get(v, v) for v in s.status], fontsize=12)
        ax.set_ylim(-0.7, len(s) - 0.3)
        # compress y axis vertically: less whitespace above/below points
        ax.margins(y=0.02)
        ax.set_xlabel("Age deviation relative to HC (years)", fontsize=12)
        ax.set_title(j,
                     fontsize=13, fontweight="normal")
        ax.spines[["top", "right"]].set_visible(False)

    # saliva spans the full left column (rows 0:3, col 0)
    ax_saliva = fig.add_subplot(gs[:, 0])
    draw_panel(ax_saliva, "saliva")

    # right column: tooth / oral mucosa / tongue stacked top-to-bottom
    # only the last (bottom) panel shows the x label to avoid overlap
    right_sites = ["tooth", "oral mucosa", "tongue"]
    for i, j in enumerate(right_sites):
        ax = fig.add_subplot(gs[i, 1])
        draw_panel(ax, j)
        if i < len(right_sites) - 1:
            ax.set_xlabel("")

    for ext in ("png", "pdf"):
        fig.savefig(OUT_DIR / f"age_disease_microbiome_age_acceleration_2col.{ext}",
                    dpi=300, bbox_inches="tight")
    print(f"Saved: {OUT_DIR / 'age_disease_microbiome_age_acceleration_2col.png'}")
    print(f"Stats: {OUT_DIR / 'age_disease_microbiome_age_acceleration_1row_stats.csv'}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import sys

BASE = Path("/beegfs/huangxiaochang/global_oral/new_analysis")
OUTDIR = BASE / "overview"
OUTDIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(BASE))

from analysis_inputs import load_metadata

meta = load_metadata(dedup=False)

meta["Geography"] = meta["Geography"].fillna("Unknown")
meta["Site_Cleaned"] = meta["Site_Cleaned"].fillna("Unknown")
meta["Status_Cleaned"] = meta["Status_Cleaned"].fillna("Unknown")

meta.groupby("Geography").size().reset_index(name="n_samples").sort_values("n_samples", ascending=False).to_csv(
    OUTDIR / "samples_by_geography.csv", index=False
)
meta.groupby("Site_Cleaned").size().reset_index(name="n_samples").sort_values("n_samples", ascending=False).to_csv(
    OUTDIR / "samples_by_site.csv", index=False
)
# --- samples_by_status.csv: disease study design ---
# 1) Pure HC BioProjects (all samples have Status_Cleaned == "Healthy Control")
pure_hc_mask = meta.groupby("BioProject")["Status_Cleaned"].transform(lambda x: (x == "Healthy Control").all())
pure_hc = meta[pure_hc_mask]
rows = [{
    "Status_Cleaned": "Healthy Control",
    "n_bioprojects": pure_hc["BioProject"].nunique(),
    "n_hc_samples": len(pure_hc),
    "n_disease_samples": 0,
    "n_total_samples": len(pure_hc),
}]

# 2) Disease studies: for each disease type, find BioProjects that contain it,
#    then count HC and disease samples within those BioProjects.
#    n_disease_samples = samples with Status_Cleaned == d (the disease itself)
#    n_other_samples   = other non-HC samples in the same BioProjects
disease_types = meta.loc[meta["Status_Cleaned"] != "Healthy Control", "Status_Cleaned"].unique()
for d in disease_types:
    bp_set = meta.loc[meta["Status_Cleaned"] == d, "BioProject"].unique()
    subset = meta[meta["BioProject"].isin(bp_set)]
    n_hc = int((subset["Status_Cleaned"] == "Healthy Control").sum())
    n_dis = int((subset["Status_Cleaned"] == d).sum())
    n_other = int(((subset["Status_Cleaned"] != "Healthy Control") & (subset["Status_Cleaned"] != d)).sum())
    rows.append({
        "Status_Cleaned": d,
        "n_bioprojects": len(bp_set),
        "n_hc_samples": n_hc,
        "n_disease_samples": n_dis,
        "n_other_samples": n_other,
        "n_total_samples": n_hc + n_dis + n_other,
    })

status_df = pd.DataFrame(rows).sort_values("n_total_samples", ascending=False)
status_df.to_csv(OUTDIR / "samples_by_status.csv", index=False)

# --- Stacked bar plot: All disease types + Healthy Control ---
# Read directly from the just-built status_df; exclude Unknown for the plot.
plot_df = status_df[status_df["Status_Cleaned"] != "Unknown"].copy()
plot_df = plot_df.sort_values("n_total_samples", ascending=False)

# Y-axis label: Disease type (n_bioprojects)
labels = [f"{row['Status_Cleaned']} ({row['n_bioprojects']})" for _, row in plot_df.iterrows()]

fig, ax = plt.subplots(figsize=(10, len(plot_df) * 0.28 + 1.5))
colors = ["#4C78A8", "#E45756"]  # blue=HC, red=disease
y_pos = np.arange(len(plot_df))

bars1 = ax.barh(y_pos, plot_df["n_hc_samples"], color=colors[0], label="Healthy Control")
bars2 = ax.barh(y_pos, plot_df["n_disease_samples"], left=plot_df["n_hc_samples"],
                color=colors[1], label="Disease")

# annotate HC / Disease counts inside bars (only if bar is wide enough)
thresh = max(plot_df["n_total_samples"]) * 0.03
for i, (_, row) in enumerate(plot_df.iterrows()):
    hc = row["n_hc_samples"]
    dis = row["n_disease_samples"]
    if hc > 0 and hc >= thresh:
        ax.text(hc / 2, i, f"{int(hc):,}", ha="center", va="center",
                color="white", fontsize=6)
    if dis > 0 and dis >= thresh:
        ax.text(hc + dis / 2, i, f"{int(dis):,}", ha="center", va="center",
                color="white", fontsize=6)

ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=7)
ax.set_xlabel("Number of samples")
ax.invert_yaxis()  # largest total at top
ax.set_title("Disease Studies + Healthy Control (sample composition)")
ax.legend(loc="lower right")
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(OUTDIR / "status_stacked_barplot.png", dpi=300)
plt.close()

# --- Stacked bar plot: Top 15 disease types (same style as full plot) ---
top15_plot_df = status_df[
    (~status_df["Status_Cleaned"].isin(["Healthy Control", "Unknown"]))
].sort_values("n_total_samples", ascending=False).head(15).copy()
top15_plot_df = top15_plot_df.sort_values("n_total_samples", ascending=True)

labels15 = [f"{row['Status_Cleaned']} ({row['n_bioprojects']})" for _, row in top15_plot_df.iterrows()]

fig, ax = plt.subplots(figsize=(10, len(top15_plot_df) * 0.45 + 1.5))
colors = ["#4C78A8", "#E45756"]  # blue=HC, red=disease
y_pos15 = np.arange(len(top15_plot_df))

bars1 = ax.barh(y_pos15, top15_plot_df["n_hc_samples"], color=colors[0], label="Healthy Control")
bars2 = ax.barh(y_pos15, top15_plot_df["n_disease_samples"], left=top15_plot_df["n_hc_samples"],
                color=colors[1], label="Disease")

# annotate HC / Disease counts inside bars
thresh15 = max(top15_plot_df["n_total_samples"]) * 0.03
for i, (_, row) in enumerate(top15_plot_df.iterrows()):
    hc = row["n_hc_samples"]
    dis = row["n_disease_samples"]
    if hc > 0 and hc >= thresh15:
        ax.text(hc / 2, i, f"{int(hc):,}", ha="center", va="center",
                color="white", fontsize=7)
    if dis > 0 and dis >= thresh15:
        ax.text(hc + dis / 2, i, f"{int(dis):,}", ha="center", va="center",
                color="white", fontsize=7)

ax.set_yticks(y_pos15)
ax.set_yticklabels(labels15, fontsize=8)
ax.set_xlabel("Number of samples")
ax.invert_yaxis()  # largest total at top
ax.set_title("Top 15 Disease Studies (sample composition)")
ax.legend(loc="lower right")
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(OUTDIR / "top15_disease_counts.png", dpi=300)
plt.close()

meta.groupby("BioProject").size().reset_index(name="n_samples").sort_values("n_samples", ascending=False).to_csv(
    OUTDIR / "samples_by_bioproject.csv", index=False
)

site_geo = meta.groupby(["Site_Cleaned", "Geography"]).size().reset_index(name="n_samples")
site_geo.to_csv(OUTDIR / "samples_by_site_geography.csv", index=False)

top20 = site_geo.sort_values("n_samples", ascending=False).head(20)
plt.figure(figsize=(7, 4.5))
ax = sns.barplot(data=top20, x="n_samples", y="Site_Cleaned", hue="Geography")
ax.set_xscale("log")
ax.set_xlabel("n_samples (log scale)")
ax.set_title("Top 20 Sites by Sample Count (colored by Geography)")
# Annotate each bar with the actual count
for p in ax.patches:
    w = p.get_width()
    if w > 0:
        ax.text(w * 1.02, p.get_y() + p.get_height() / 2,
                f"{int(w):,}", ha="left", va="center", fontsize=7)
plt.tight_layout()
plt.savefig(OUTDIR / "top_site_geography_counts.png", dpi=300)
plt.close()

heat = site_geo.pivot(index="Site_Cleaned", columns="Geography", values="n_samples").fillna(0)
heat_log = np.log1p(heat)  # log(1 + x) to handle zeros gracefully
plt.figure(figsize=(7, 4.5))
sns.heatmap(heat_log, cmap="Blues", annot=heat.astype(int), fmt="d",
            annot_kws={"fontsize": 7}, linewidths=0.5,
            cbar_kws={"label": "log1p(n_samples)"})
plt.title("Site × Geography Sample Count (annotated = actual, color = log scale)")
plt.xlabel("")
plt.ylabel("")
plt.tight_layout()
plt.savefig(OUTDIR / "site_geography_heatmap.png", dpi=300)
plt.close()

# --- DEDUP: Stacked bar plot: Top 15 disease types ---
meta_dedup = load_metadata(dedup=True)
meta_dedup["Status_Cleaned"] = meta_dedup["Status_Cleaned"].fillna("Unknown")

# Pure HC
dedup_pure_hc_mask = meta_dedup.groupby("BioProject")["Status_Cleaned"].transform(lambda x: (x == "Healthy Control").all())
dedup_pure_hc = meta_dedup[dedup_pure_hc_mask]
dedup_rows = [{
    "Status_Cleaned": "Healthy Control",
    "n_bioprojects": dedup_pure_hc["BioProject"].nunique(),
    "n_hc_samples": len(dedup_pure_hc),
    "n_disease_samples": 0,
    "n_total_samples": len(dedup_pure_hc),
}]

dedup_disease_types = meta_dedup.loc[meta_dedup["Status_Cleaned"] != "Healthy Control", "Status_Cleaned"].unique()
for d in dedup_disease_types:
    bp_set = meta_dedup.loc[meta_dedup["Status_Cleaned"] == d, "BioProject"].unique()
    subset = meta_dedup[meta_dedup["BioProject"].isin(bp_set)]
    n_hc = int((subset["Status_Cleaned"] == "Healthy Control").sum())
    n_dis = int((subset["Status_Cleaned"] == d).sum())
    dedup_rows.append({
        "Status_Cleaned": d,
        "n_bioprojects": len(bp_set),
        "n_hc_samples": n_hc,
        "n_disease_samples": n_dis,
        "n_total_samples": n_hc + n_dis,
    })

dedup_status_df = pd.DataFrame(dedup_rows).sort_values("n_total_samples", ascending=False)

dedup_top15 = dedup_status_df[
    (~dedup_status_df["Status_Cleaned"].isin(["Healthy Control", "Unknown"]))
].sort_values("n_total_samples", ascending=False).head(15).copy()
dedup_top15 = dedup_top15.sort_values("n_total_samples", ascending=True)

dedup_labels = [f"{row['Status_Cleaned']} ({row['n_bioprojects']})" for _, row in dedup_top15.iterrows()]

fig, ax = plt.subplots(figsize=(10, len(dedup_top15) * 0.45 + 1.5))
colors = ["#4C78A8", "#E45756"]
y_pos = np.arange(len(dedup_top15))

ax.barh(y_pos, dedup_top15["n_hc_samples"], color=colors[0], label="Healthy Control")
ax.barh(y_pos, dedup_top15["n_disease_samples"], left=dedup_top15["n_hc_samples"],
        color=colors[1], label="Disease")

thresh = max(dedup_top15["n_total_samples"]) * 0.03
for i, (_, row) in enumerate(dedup_top15.iterrows()):
    hc = row["n_hc_samples"]
    dis = row["n_disease_samples"]
    if hc > 0 and hc >= thresh:
        ax.text(hc / 2, i, f"{int(hc):,}", ha="center", va="center", color="white", fontsize=7)
    if dis > 0 and dis >= thresh:
        ax.text(hc + dis / 2, i, f"{int(dis):,}", ha="center", va="center", color="white", fontsize=7)

ax.set_yticks(y_pos)
ax.set_yticklabels(dedup_labels, fontsize=8)
ax.set_xlabel("Number of samples")
ax.invert_yaxis()
ax.set_title("Top 15 Disease Studies (dedup, sample composition)")
ax.legend(loc="lower right")
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(OUTDIR / "top15_disease_counts_dedup.png", dpi=300)
plt.close()

print("Done.")


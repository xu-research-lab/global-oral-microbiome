#!/usr/bin/env python3
"""
Generate publication figures for HOMD uncapture analysis.

Figures:
  1. ASV-weighted + read-weighted HOMD capture across identity thresholds
  2. Subject prevalence distribution: identity >=97% vs <97%
  3. Analysis pipeline funnel (horizontal bar)
  4. Cross-study robustness scatter (uncaptured ASVs)
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

# ============================================================
# Style
# ============================================================
sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

COLORS = {
    "ASV-weighted": "#6C5B7B",
    "Read-weighted": "#C06C84",
    "captured": "#355C7D",
    "uncaptured": "#C06C84",
    "low_identity": "#e74c3c",
    "no_hit": "#f39c12",
    "low_qcov": "#95a5a6",
}

THRESHOLDS = [80.0, 85.0, 90.0, 95.0, 97.0, 99.0, 99.5, 100.0]


def log(msg: str) -> None:
    print(f"[plot] {msg}", flush=True)


# ============================================================
# Load data
# ============================================================
def load_data(tsv_path: str) -> pd.DataFrame:
    log(f"Loading {tsv_path} ...")
    df = pd.read_csv(
        tsv_path,
        sep="\t",
        usecols=["n_asvs", "total_reads", "prevalence_samples",
                  "n_bioprojects", "identity", "uncaptured_reason",
                  "subject_prevalence"],
        dtype={
            "n_asvs": "float64",
            "total_reads": "float64",
            "prevalence_samples": "float64",
            "n_bioprojects": "int32",
            "identity": "float32",
            "uncaptured_reason": "str",
            "subject_prevalence": "float64",
        },
    )
    log(f"  {len(df):,} rows loaded")
    return df


# ============================================================
# Figure 1: ASV-weighted + read-weighted identity capture curve
# ============================================================
def fig1_identity_capture(df: pd.DataFrame, outdir: Path) -> None:
    log("Fig 1: Identity threshold capture curve ...")

    total_asvs = df["n_asvs"].sum()
    total_reads = df["total_reads"].sum()

    rows = []
    labels = {
        80.0: "≥80%", 85.0: "≥85%", 90.0: "≥90%", 95.0: "≥95%",
        97.0: "≥97%", 99.0: "≥99%", 99.5: "≥99.5%", 100.0: "100%",
    }

    for t in THRESHOLDS:
        mask = df["identity"] >= t
        asv_cap = df.loc[mask, "n_asvs"].sum()
        read_cap = df.loc[mask, "total_reads"].sum()
        rows.append({
            "identity_threshold": t,
            "identity_label": labels[t],
            "asv_captured": asv_cap,
            "asv_fraction": asv_cap / total_asvs,
            "read_captured": read_cap,
            "read_fraction": read_cap / total_reads,
            "total_asvs": total_asvs,
            "total_reads": total_reads,
        })

    summary = pd.DataFrame(rows)

    # Write TSV
    tsv_path = outdir / "fig1_asv_weighted_identity_capture.tsv"
    summary.to_csv(tsv_path, sep="\t", index=False)
    log(f"  -> {tsv_path}")

    # Plot
    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    x_labels = summary["identity_label"].values
    y_asv = summary["asv_fraction"].values
    y_read = summary["read_fraction"].values

    ax.plot(x_labels, y_asv,
            color=COLORS["ASV-weighted"], linewidth=2.5, marker="o",
            markersize=7, label="ASV-weighted", zorder=5)
    ax.plot(x_labels, y_read,
            color=COLORS["Read-weighted"], linewidth=2.0, marker="s",
            markersize=6, linestyle="--", alpha=0.8, label="Read-weighted", zorder=4)

    # Annotate key thresholds
    for _, row in summary.iterrows():
        if row["identity_threshold"] in (97.0, 100.0):
            ax.annotate(
                f'{row["asv_fraction"]:.1%}',
                xy=(row["identity_label"], row["asv_fraction"]),
                xytext=(0, 12), textcoords="offset points",
                ha="center", fontsize=9, fontweight="bold",
                color=COLORS["ASV-weighted"],
            )

    ax.set_ylim(0, 1.05)
    ax.set_xlabel("HOMD identity threshold", fontsize=12)
    ax.set_ylabel("Captured fraction of total ASVs / reads", fontsize=12)
    ax.set_title("Global HOMD capture across identity thresholds", fontsize=13, fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    ax.legend(title="Weighting", frameon=True, fontsize=10, title_fontsize=10)
    fig.tight_layout()

    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"fig1_asv_weighted_identity_capture.{ext}",
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"  -> fig1_asv_weighted_identity_capture.pdf/.png")


# ============================================================
# Figure 2: Subject prevalence distribution — identity >=97% vs <97%
# ============================================================
def fig2_prevalence_distribution(df: pd.DataFrame, outdir: Path) -> None:
    log("Fig 2: Prevalence distribution by identity group ...")

    # Assign groups
    df = df.copy()
    df["identity_group"] = np.where(
        df["identity"] >= 97.0, "≥97% (HOMD-captured)", "<97% (uncaptured)"
    )

    # Aggregate per row (avoid iterrows for speed)
    prev = df["subject_prevalence"].fillna(0.5).clip(lower=0.3, upper=20000).values
    grp = df["identity_group"].values
    counts = df["n_asvs"].astype(int).values

    rows = []
    for p, g, c in zip(prev, grp, counts):
        if c > 0:
            pts = [(p, g)] * c
            rows.extend(pts)
    plot_df = pd.DataFrame(rows, columns=["subject_prevalence", "identity_group"])
    log(f"  Built {len(plot_df):,} ASV-weighted rows")

    # Stats
    stats = {}
    for g in plot_df["identity_group"].unique():
        sub = plot_df[plot_df["identity_group"] == g]
        vals = sub["subject_prevalence"].values
        stats[g] = {"median": float(np.median(vals)), "mean": float(np.mean(vals)), "n": len(vals)}
        log(f"  {g}: {len(vals):,} ASVs, median prevalence={stats[g]['median']:.0f}")

    # Plot: single-panel ECDF
    fig, ax = plt.subplots(figsize=(8, 5))

    groups = ["<97% (uncaptured)", "≥97% (HOMD-captured)"]
    colors = [COLORS["uncaptured"], COLORS["captured"]]

    for g, color in zip(groups, colors):
        sub = plot_df[plot_df["identity_group"] == g]
        vals = np.sort(sub["subject_prevalence"].values)
        ecdf = np.arange(1, len(vals) + 1) / len(vals)
        ax.step(vals, ecdf, color=color, linewidth=2.2, where="post", label=g)

        # Median marker
        med = stats[g]["median"]
        ax.axvline(x=med, color=color, linestyle=":", linewidth=1.5, alpha=0.7)
        ax.text(med, 0.02, f"median={med:.0f}", color=color, fontsize=9,
                ha="center", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7, edgecolor="none"))

    ax.set_xscale("log")
    ax.set_xlabel("Subject prevalence (number of samples)", fontsize=12)
    ax.set_ylabel("Cumulative fraction of ASVs", fontsize=12)
    ax.set_title("Subject prevalence distribution by HOMD identity", fontsize=13, fontweight="bold")
    ax.set_xlim(0.3, 20000)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=10, frameon=True, loc="center right")

    # KS-test annotation
    from scipy.stats import ks_2samp
    grp_a = plot_df[plot_df["identity_group"] == groups[0]]["subject_prevalence"].values
    grp_b = plot_df[plot_df["identity_group"] == groups[1]]["subject_prevalence"].values
    ks_stat, ks_p = ks_2samp(grp_a, grp_b)
    p_str = f"{ks_p:.2e}" if ks_p >= 1e-300 else "< 2.2e-16"
    ax.text(0.98, 0.98,
            f"KS stat = {ks_stat:.3f}\np = {p_str}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"fig2_prevalence_by_identity.{ext}",
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"  -> fig2_prevalence_by_identity.pdf/.png")


# ============================================================
# Figure 3: Analysis pipeline funnel (horizontal bar)
# ============================================================
def fig3_pipeline_funnel(outdir: Path, base_dir: Path) -> None:
    log("Fig 3: Pipeline funnel ...")

    # Funnel steps from pipeline_report_comprehensive.md (Funnel Overview)
    # Values verified from data files and the report
    total_derep = 2_844_174
    steps = [
        # (label, value, color, phase_label)
        ("Exact-derep\nsequences", 2_844_174, "#2c3e50", "Phase 1 input"),
        ("HOMD-uncaptured\n(id<97% or no hit, no qcov)", 1_843_301, "#2980b9", "HOMD mapping"),
        ("Support floor\n(≥2 BP, ≥50 reads)", 69_820, "#16a085", "prevalence filter"),
        ("Chimera removal\n(vsearch uchime_ref)", 54_254, "#27ae60", "chimera removal"),
        ("Contamination flags\n(2-stage, flags only)", 54_254, "#8e44ad", "Phase 2 firewall"),
        ("SEPP placement\n(GG2 backbone)", 54_132, "#d35400", "Phase 3 SEPP"),
        ("Phylotype clustering\n(tree-edge)", 23_088, "#e67e22", "phylotype cluster"),
        ("Phylotype gates\n(gates 1-5)", 56, "#c0392b", "Phase 4 gates"),
        ("novel_all\n(<97% to all 3 DBs)", 14, "#7f8c8d", "novelty"),
        ("known_not_in_homd\n(HOMD gap)", 42, "#95a5a6", "novelty"),
    ]

    fig, ax = plt.subplots(figsize=(11, 8))

    n = len(steps)
    y_positions = np.arange(n)[::-1]  # top to bottom

    # Display pass-through / filter type
    is_filter = [True] * n  # all are reductions except contamination flags
    is_filter[4] = False  # contamination flags - no removal
    is_filter[8] = False  # novel_all - split from 56
    is_filter[9] = False  # known_not_in_homd - split from 56

    max_val = float(steps[0][1])
    for i, (label, val, color, phase) in enumerate(steps):
        # Bar width proportional to value
        width = val / max_val
        # For the split 56 -> 14 + 42, show them as a branch
        if i in (8, 9):
            # These are the two sub-branches of 56; start from the 56 level width
            parent_width = steps[7][1] / max_val
            bar = ax.barh(y_positions[i], width, height=0.6,
                          color=color, alpha=0.85, edgecolor="white", linewidth=0.5,
                          left=parent_width - width if i == 9 else 0)
        else:
            ax.barh(y_positions[i], width, height=0.6,
                    color=color, alpha=0.85, edgecolor="white", linewidth=0.5)

        ax.text(width + 0.01, y_positions[i],
                f" {val:,}",
                va="center", fontsize=10, fontweight="bold", color="#333333")

    # Y-axis labels
    labels = [s[0] for s in steps]
    phases = [s[3] for s in steps]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Relative count (bar width proportional)", fontsize=11)
    ax.set_title("Analysis pipeline: from 2.84M derep sequences to 56 candidate novel phylotypes",
                 fontsize=13, fontweight="bold")

    # Add phase annotations on the right side
    for i, (label, val, color, phase) in enumerate(steps):
        ax.text(1.02, y_positions[i], phase, va="center", fontsize=8,
                color=color, fontweight="bold", transform=ax.get_yaxis_transform())

    ax.set_xlim(0, 1.18)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Add a connecting line note for the 14+42 split
    ax.annotate("", xy=(0, y_positions[7]), xytext=(0, y_positions[8] + 0.3),
                arrowprops=dict(arrowstyle="-", color="#7f8c8d", lw=0.5, alpha=0.5))
    ax.annotate("", xy=(0, y_positions[7]), xytext=(0, y_positions[9] + 0.3),
                arrowprops=dict(arrowstyle="-", color="#7f8c8d", lw=0.5, alpha=0.5))

    fig.subplots_adjust(left=0.28, right=0.72, top=0.93, bottom=0.08)

    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"fig3_pipeline_funnel.{ext}",
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"  -> fig3_pipeline_funnel.pdf/.png")


# ============================================================
# Figure 4: Cross-study robustness scatter
# ============================================================
def fig4_cross_study_robustness(df: pd.DataFrame, outdir: Path) -> None:
    log("Fig 4: Cross-study robustness scatter ...")

    # Only uncaptured sequences
    unc = df[df["uncaptured_reason"] != "captured"].copy()
    log(f"  Uncaptured sequences: {len(unc):,}")

    fig, ax = plt.subplots(figsize=(9, 6))

    category_order = ["low_identity", "no_hit", "low_qcov"]
    palette = {k: COLORS[k] for k in category_order if k in unc["uncaptured_reason"].values}
    markers = {"low_identity": "o", "no_hit": "D", "low_qcov": "v"}
    labels_cn = {
        "low_identity": "Low identity",
        "no_hit": "No HOMD hit",
        "low_qcov": "Low query cov",
    }

    for cat in category_order:
        sub = unc[unc["uncaptured_reason"] == cat]
        if len(sub) == 0:
            continue
        sizes = np.clip(sub["prevalence_samples"].values, 1, None)
        sizes_log = np.log10(sizes + 1) * 3
        ax.scatter(
            sub["n_bioprojects"], sub["total_reads"],
            s=sizes_log, c=palette[cat], alpha=0.35,
            marker=markers.get(cat, "o"), label=labels_cn.get(cat, cat),
            edgecolors="none", rasterized=True,
        )

    ax.set_yscale("log")
    ax.set_xlabel("Number of BioProjects", fontsize=12)
    ax.set_ylabel("Total reads (log scale)", fontsize=12)
    ax.set_title("Cross-study robustness of HOMD-uncaptured ASVs", fontsize=13, fontweight="bold")

    # Reference lines
    ax.axvline(x=3, color="#2c3e50", linestyle="--", linewidth=1, alpha=0.5, label="_nolegend_")
    ax.text(3.2, ax.get_ylim()[1] * 0.5 if ax.get_ylim()[1] > 0 else 1e6,
            "≥3 BioProjects", fontsize=8, color="#2c3e50", alpha=0.7, rotation=90, va="center")

    legend_handles = []
    for cat in category_order:
        if cat in palette:
            legend_handles.append(
                plt.Line2D([0], [0], marker=markers.get(cat, "o"), color="w",
                           markerfacecolor=palette[cat], markersize=8,
                           label=labels_cn.get(cat, cat))
            )
    legend_handles.append(
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
                   markersize=8, label="Bubble size = Prevalence (log)")
    )
    ax.legend(handles=legend_handles, fontsize=9, loc="upper left", frameon=True,
              title="Category", title_fontsize=10)

    # Annotate top sequences
    top = unc.nlargest(8, "total_reads")
    for _, row in top.iterrows():
        ax.annotate(
            f'{row["n_bioprojects"]}BPs',
            xy=(row["n_bioprojects"], row["total_reads"]),
            xytext=(8, 0), textcoords="offset points",
            fontsize=7, color="#333333", alpha=0.8,
        )
    fig.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.10)

    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"fig4_cross_study_robustness.{ext}",
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"  -> fig4_cross_study_robustness.pdf/.png")


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basedir",
                        default="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
    parser.add_argument("--outdir", default=None,
                        help="Output directory for figures (default: basedir/figures)")
    args = parser.parse_args()

    base_dir = Path(args.basedir)
    outdir = Path(args.outdir) if args.outdir else base_dir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)

    tsv_path = base_dir / "homd_mapping" / "exact_derep_homd_annotated.tsv"
    df = load_data(str(tsv_path))

    fig1_identity_capture(df, outdir)
    fig2_prevalence_distribution(df, outdir)
    fig3_pipeline_funnel(outdir, base_dir)
    fig4_cross_study_robustness(df, outdir)

    log("All figures done.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Plot global HOMD capture over full 80-100% identity thresholds at qcov>=99%."""

from __future__ import annotations

import argparse
import os

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
THRESHOLDS = [80.0, 85.0, 90.0, 95.0, 97.0, 99.0, 99.5, 100.0]

sns.set_theme(style="whitegrid", context="paper", font_scale=1.25)
plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

COLORS = {
    "Exact sequence": "#355C7D",
    "ASV-weighted": "#6C5B7B",
    "Read-weighted": "#C06C84",
    "Prevalence-weighted": "#F67280",
}


def log(msg):
    print(f"[full80] {msg}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--chunksize", type=int, default=200000)
    args = parser.parse_args()

    in_path = os.path.join(args.outdir, "homd_coverage_by_derep_sequence.tsv")
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    out_base = os.path.join(figdir, "fig1_global_identity_thresholds_qcov99_full80to100")
    summary_path = os.path.join(args.outdir, "global_identity_thresholds_qcov99_full80to100.tsv")

    totals = {"derep": 0, "asvs": 0.0, "reads": 0.0, "prevalence": 0.0}
    captured = {t: {"derep": 0, "asvs": 0.0, "reads": 0.0, "prevalence": 0.0} for t in THRESHOLDS}

    usecols = ["identity", "query_cov", "n_asvs", "total_reads", "prevalence_samples"]
    for chunk in pd.read_csv(in_path, sep="\t", usecols=usecols, chunksize=args.chunksize):
        for col in usecols:
            chunk[col] = pd.to_numeric(chunk[col], errors="coerce").fillna(0)
        totals["derep"] += len(chunk)
        totals["asvs"] += float(chunk["n_asvs"].sum())
        totals["reads"] += float(chunk["total_reads"].sum())
        totals["prevalence"] += float(chunk["prevalence_samples"].sum())
        q = chunk["query_cov"] >= 99
        for t in THRESHOLDS:
            m = q & (chunk["identity"] >= t)
            captured[t]["derep"] += int(m.sum())
            captured[t]["asvs"] += float(chunk.loc[m, "n_asvs"].sum())
            captured[t]["reads"] += float(chunk.loc[m, "total_reads"].sum())
            captured[t]["prevalence"] += float(chunk.loc[m, "prevalence_samples"].sum())
        log(f"processed {totals['derep']:,} rows")

    rows = []
    labels = {80.0: ">=80%", 85.0: ">=85%", 90.0: ">=90%", 95.0: ">=95%", 97.0: ">=97%", 99.0: ">=99%", 99.5: ">=99.5%", 100.0: "100%"}
    metric_labels = {"derep": "Exact sequence", "asvs": "ASV-weighted", "reads": "Read-weighted", "prevalence": "Prevalence-weighted"}
    for t in THRESHOLDS:
        for metric, label in metric_labels.items():
            rows.append({
                "identity_threshold": t,
                "identity_label": labels[t],
                "metric": metric,
                "metric_label": label,
                "captured_count": captured[t][metric],
                "denominator": totals[metric],
                "captured_fraction": captured[t][metric] / totals[metric] if totals[metric] else 0,
            })
    df = pd.DataFrame(rows)
    df.to_csv(summary_path, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    sns.lineplot(
        data=df,
        x="identity_label",
        y="captured_fraction",
        hue="metric_label",
        style="metric_label",
        markers=True,
        dashes=False,
        linewidth=2.2,
        markersize=7,
        palette=COLORS,
        ax=ax,
    )
    ax.set_ylim(0, 1)
    ax.set_xlabel("HOMD identity threshold (query coverage >=99%)")
    ax.set_ylabel("Captured fraction")
    ax.set_title("Global HOMD capture across 80-100% identity thresholds")
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.legend(title="Weighting", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    fig.savefig(out_base + ".png", bbox_inches="tight")
    log(f"wrote {out_base}.pdf/.png")
    log(f"wrote {summary_path}")


if __name__ == "__main__":
    main()

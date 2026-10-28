#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s04_beta_diversity_pcoa_permanova_paired_adonis2.py

Paired-sample PCoA + PERMANOVA for oral-gut beta diversity.

Design:
  1. Use only complete paired samples:
       exactly one Oral + one Gut within each BioProject x pair_id by default.

  2. PCoA:
       computed from precomputed sample x sample distance matrices.

  3. PERMANOVA:
       calls R vegan::adonis2().

       Per-BioProject:
         D ~ site
         permutation restricted within pair_id blocks.

       Pooled:
         D ~ bioproject + site
         site is tested after BioProject when by = "terms".
         permutation restricted within BioProject x pair_id blocks.

Default input:
  03.beta_diversity/distance_matrix/
    braycurtis_filter_ra_distance_matrix.tsv
    weighted_unifrac_filter_distance_matrix.tsv
    unweighted_unifrac_filter_distance_matrix.tsv

  01.filter_data/filter_metadata.csv

Default output:
  03.beta_diversity/pcoa_paired_adonis2/

Required R packages:
  vegan
  permute

Before running on the server:
  module load R/4.1.3
  export R_LIBS_USER=~/Rlibs

Notes:
  - This script does not recalculate Bray-Curtis or UniFrac distances.
  - The distance matrices must be full square sample x sample matrices.
  - For pooled adonis2, formula order is bioproject + site.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_ROOT = Path("/beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp")
DEFAULT_DISTANCE_DIR = DEFAULT_ROOT / "03.beta_diversity" / "distance_matrix"
DEFAULT_METADATA = DEFAULT_ROOT / "01.filter_data" / "filter_metadata.csv"
DEFAULT_OUT_DIR = DEFAULT_ROOT / "03.beta_diversity" / "pcoa_paired_adonis2"


# -----------------------------
# Helpers
# -----------------------------

def parse_list(x):
    if x is None:
        return []
    x = str(x).strip()
    if x == "":
        return []
    return [i.strip() for i in x.split(",") if i.strip()]


def parse_bool(x):
    if isinstance(x, bool):
        return x
    s = str(x).strip().lower()
    if s in {"true", "t", "1", "yes", "y"}:
        return True
    if s in {"false", "f", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Cannot parse boolean value: {x}")


def clean_for_filename(x):
    x = str(x)
    x = re.sub(r"[^\w.\-]+", "_", x)
    x = x.strip("_")
    return x if x else "NA"


def read_table(path):
    path = str(path)
    if path.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def standardize_site(x):
    x0 = str(x).strip()
    xl = x0.lower()

    site_map = {
        "oral": "Oral",
        "mouth": "Oral",
        "saliva": "Oral",
        "tongue": "Oral",
        "buccal": "Oral",
        "dental": "Oral",
        "plaque": "Oral",
        "gut": "Gut",
        "stool": "Gut",
        "fecal": "Gut",
        "faecal": "Gut",
        "feces": "Gut",
        "faeces": "Gut",
    }

    return site_map.get(xl, x0)


def normalize_status_value(x):
    x = str(x).strip()
    if x == "" or x.lower() in {"nan", "none", "na", "null"}:
        return "NA"
    return x


def is_hc_status(x):
    s = str(x).strip().lower()
    return s in {
        "hc", "h", "healthy", "health", "control", "controls", "normal",
        "healthy_control", "healthy control", "non-disease", "non_disease",
        "non disease", "non-diseased", "non_diseased", "non diseased",
        "noncase", "non-case", "non_case",
    }


def status_sort_key(x):
    x = normalize_status_value(x)
    if x == "NA":
        return (2, x)
    if is_hc_status(x):
        return (0, x)
    return (1, x)


def classify_bp_status_group(status_values):
    values = [normalize_status_value(x) for x in status_values]
    values = [x for x in values if x != "NA"]

    if not values:
        return 2, "disease_study"

    has_hc = any(is_hc_status(x) for x in values)
    has_disease = any(not is_hc_status(x) for x in values)

    if has_hc and not has_disease:
        return 0, "HC"
    if has_hc and has_disease:
        return 1, "HC+disease"
    return 2, "disease_study"


def metric_display_name(metric):
    mapping = {
        "braycurtis": "Bray-Curtis",
        "weighted_unifrac": "Weighted UniFrac",
        "unweighted_unifrac": "Unweighted UniFrac",
    }
    return mapping.get(metric, metric)


def metric_to_dataset_label(metric):
    if metric == "braycurtis":
        return "filter_ra"
    if metric in {"weighted_unifrac", "unweighted_unifrac"}:
        return "filter"
    return "unknown"


def format_p(p):
    try:
        p = float(p)
    except Exception:
        return "NA"
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def p_to_stars(p):
    try:
        p = float(p)
    except Exception:
        return "NA"
    if not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def json_safe_records(df):
    tmp = df.copy()
    tmp = tmp.replace({np.nan: None})
    records = tmp.to_dict(orient="records")

    def convert(v):
        if isinstance(v, np.integer):
            return int(v)
        if isinstance(v, np.floating):
            if np.isnan(v):
                return None
            return float(v)
        if isinstance(v, np.bool_):
            return bool(v)
        return v

    return [{k: convert(v) for k, v in rec.items()} for rec in records]


def validate_columns(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} missing required columns:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns:\n  "
            + "\n  ".join(df.columns)
        )


# -----------------------------
# Metadata and paired sample filter
# -----------------------------

def build_bp_status_labels_and_order(md, bioproject_col, status_col):
    rows = []
    labels = {}

    for bp, sub in md.groupby(bioproject_col, dropna=False):
        if status_col in sub.columns:
            statuses = (
                sub[status_col]
                .dropna()
                .astype(str)
                .map(normalize_status_value)
                .unique()
                .tolist()
            )
        else:
            statuses = ["NA"]

        statuses_clean = [s for s in statuses if s != "NA"]
        status_text = "NA" if not statuses_clean else "+".join(sorted(statuses_clean, key=status_sort_key))

        status_rank, status_group = classify_bp_status_group(statuses)

        labels[bp] = f"{bp}({status_text})"

        rows.append({
            "bioproject": bp,
            "bp_status_label": labels[bp],
            "status_rank": status_rank,
            "status_group_for_order": status_group,
            "status_values": status_text,
        })

    order_df = pd.DataFrame(rows).sort_values(
        ["status_rank", "bioproject"],
        ascending=[True, True],
    )
    order_df["bp_plot_order"] = np.arange(order_df.shape[0], dtype=int)

    return labels, order_df


def prepare_metadata(metadata, sample_col, bioproject_col, site_col, pair_col, status_col):
    md = metadata.copy()

    required = [sample_col, bioproject_col, site_col, pair_col]
    validate_columns(md, required, "metadata")

    if status_col not in md.columns:
        md[status_col] = "NA"

    for c in [sample_col, bioproject_col, site_col, pair_col, status_col]:
        md[c] = md[c].astype(str).str.strip()

    md = md.drop_duplicates(subset=[sample_col]).copy()

    md["sample_id_for_analysis"] = md[sample_col]
    md["bioproject_for_analysis"] = md[bioproject_col]
    md["site_for_analysis"] = md[site_col].apply(standardize_site)
    md["pair_id_for_analysis"] = md[pair_col]
    md["status_for_analysis"] = md[status_col]

    md = md.loc[md["site_for_analysis"].isin(["Oral", "Gut"]), :].copy()

    bp_status_labels, bp_order_df = build_bp_status_labels_and_order(
        md,
        "bioproject_for_analysis",
        "status_for_analysis",
    )

    md = md.merge(
        bp_order_df[["bioproject", "bp_status_label", "status_group_for_order", "bp_plot_order"]],
        left_on="bioproject_for_analysis",
        right_on="bioproject",
        how="left",
    )
    md = md.drop(columns=["bioproject"], errors="ignore")

    md["pair_block_for_permutation"] = (
        md["bioproject_for_analysis"].astype(str)
        + "__"
        + md["pair_id_for_analysis"].astype(str)
    )

    return md, bp_status_labels, bp_order_df


def filter_complete_pairs(md, strict_one_to_one=True):
    rows = []
    keep_blocks = []

    for (bp, pair_id), sub in md.groupby(["bioproject_for_analysis", "pair_id_for_analysis"], dropna=False):
        n_oral = int((sub["site_for_analysis"] == "Oral").sum())
        n_gut = int((sub["site_for_analysis"] == "Gut").sum())
        n_samples = int(sub["sample_id_for_analysis"].nunique())

        if strict_one_to_one:
            keep = (n_oral == 1 and n_gut == 1 and n_samples == 2)
            rule = "exactly_one_oral_one_gut"
        else:
            keep = (n_oral >= 1 and n_gut >= 1)
            rule = "at_least_one_oral_one_gut"

        pair_block = str(bp) + "__" + str(pair_id)

        rows.append({
            "bioproject": bp,
            "pair_id": pair_id,
            "pair_block_for_permutation": pair_block,
            "n_samples": n_samples,
            "n_oral": n_oral,
            "n_gut": n_gut,
            "keep_paired": bool(keep),
            "pair_filter_rule": rule,
        })

        if keep:
            keep_blocks.append(pair_block)

    summary_df = pd.DataFrame(rows)
    out = md.loc[md["pair_block_for_permutation"].isin(set(keep_blocks)), :].copy()

    return out, summary_df


# -----------------------------
# Distance matrix and PCoA
# -----------------------------

def read_distance_matrix(path):
    """
    Read square distance matrix with sample IDs in both row index and columns.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", index_col=0)

    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    common = [sid for sid in df.index if sid in set(df.columns)]
    df = df.loc[common, common].copy()

    df = df.apply(pd.to_numeric, errors="coerce")

    if df.isna().values.any():
        n_na = int(df.isna().sum().sum())
        raise RuntimeError(f"Distance matrix contains {n_na} NA values: {path}")

    arr = df.values.astype(float)

    if not np.allclose(arr, arr.T, atol=1e-8, rtol=1e-8):
        max_diff = float(np.max(np.abs(arr - arr.T)))
        print(f"WARNING: distance matrix is not perfectly symmetric; max diff={max_diff}. Symmetrizing.", flush=True)
        arr = (arr + arr.T) / 2.0

    np.fill_diagonal(arr, 0.0)
    df = pd.DataFrame(arr, index=df.index, columns=df.columns)

    return df


def subset_dm_and_md(dm_df, md, bp=None):
    sample_set = set(dm_df.index)
    sub_md = md.loc[md["sample_id_for_analysis"].isin(sample_set), :].copy()

    if bp is not None:
        sub_md = sub_md.loc[sub_md["bioproject_for_analysis"].astype(str).eq(str(bp)), :].copy()

    sample_ids = [sid for sid in dm_df.index if sid in set(sub_md["sample_id_for_analysis"])]

    if len(sample_ids) == 0:
        return np.zeros((0, 0), dtype=float), [], sub_md.iloc[0:0].copy()

    sub_dm = dm_df.loc[sample_ids, sample_ids].values.astype(float)
    sub_md = sub_md.set_index("sample_id_for_analysis").loc[sample_ids].reset_index()

    return sub_dm, sample_ids, sub_md


def run_pcoa(D, sample_ids):
    try:
        from skbio import DistanceMatrix
        from skbio.stats.ordination import pcoa
    except ImportError as e:
        raise ImportError(
            "Cannot import scikit-bio PCoA functions. Install scikit-bio first."
        ) from e

    dm = DistanceMatrix(D, ids=sample_ids)
    ord_res = pcoa(dm)

    coords = ord_res.samples.copy()
    coords["sample_id_for_analysis"] = coords.index.astype(str)

    if "PC1" not in coords.columns:
        coords["PC1"] = np.nan
    if "PC2" not in coords.columns:
        coords["PC2"] = np.nan

    var = ord_res.proportion_explained.copy()
    pc1_pct = float(var.iloc[0]) * 100 if len(var) > 0 else np.nan
    pc2_pct = float(var.iloc[1]) * 100 if len(var) > 1 else np.nan

    return coords, pc1_pct, pc2_pct


# -----------------------------
# R adonis2 backend
# -----------------------------

R_ADONIS2_SCRIPT = r"""
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 9) {
  stop("Usage: Rscript run_adonis2.R <dm_file> <meta_file> <out_file> <formula_rhs> <nperm> <seed> <block_col> <by_arg> <log_file>")
}

dm_file <- args[[1]]
meta_file <- args[[2]]
out_file <- args[[3]]
formula_rhs <- args[[4]]
nperm <- as.integer(args[[5]])
seed <- as.integer(args[[6]])
block_col <- args[[7]]
by_arg <- args[[8]]
log_file <- args[[9]]

sink(log_file, split = FALSE)
cat("R adonis2 backend\n")
cat("dm_file: ", dm_file, "\n", sep = "")
cat("meta_file: ", meta_file, "\n", sep = "")
cat("out_file: ", out_file, "\n", sep = "")
cat("formula_rhs: ", formula_rhs, "\n", sep = "")
cat("nperm: ", nperm, "\n", sep = "")
cat("seed: ", seed, "\n", sep = "")
cat("block_col: ", block_col, "\n", sep = "")
cat("by_arg: ", by_arg, "\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")
cat(".libPaths():\n")
print(.libPaths())

suppressPackageStartupMessages({
  library(vegan)
  library(permute)
})

cat("vegan version: ", as.character(packageVersion("vegan")), "\n", sep = "")
cat("permute version: ", as.character(packageVersion("permute")), "\n", sep = "")

dm <- read.table(
  dm_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

md <- read.table(
  meta_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

dm_mat <- as.matrix(dm)
storage.mode(dm_mat) <- "double"

if (!all(rownames(dm_mat) == colnames(dm_mat))) {
  stop("Distance matrix row names and column names are not identical.")
}

D <- as.dist(dm_mat)
ids <- attr(D, "Labels")

md <- md[match(ids, md$sample_id_for_analysis), , drop = FALSE]

if (any(is.na(md$sample_id_for_analysis))) {
  stop("Some distance matrix samples were not found in metadata.")
}

if (!all(md$sample_id_for_analysis == ids)) {
  stop("Metadata order does not match distance labels after matching.")
}

md$site <- factor(md$site_for_analysis)
md$bioproject <- factor(md$bioproject_for_analysis)

if (!(block_col %in% colnames(md))) {
  stop(paste("Block column not found in metadata:", block_col))
}

md$perm_block <- factor(md[[block_col]])

cat("n samples: ", nrow(md), "\n", sep = "")
cat("site table:\n")
print(table(md$site))
cat("n bioprojects: ", length(unique(md$bioproject)), "\n", sep = "")
cat("n permutation blocks: ", length(unique(md$perm_block)), "\n", sep = "")
cat("permutation block size table:\n")
print(table(table(md$perm_block)))

set.seed(seed)

perm <- how(
  nperm = nperm,
  blocks = md$perm_block
)

formula_txt <- paste("D ~", formula_rhs)
cat("formula: ", formula_txt, "\n", sep = "")

fit <- adonis2(
  as.formula(formula_txt),
  data = md,
  permutations = perm,
  by = by_arg
)

out <- as.data.frame(fit)
out$term <- rownames(out)
out <- out[, c("term", setdiff(colnames(out), "term")), drop = FALSE]

write.table(
  out,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done.\n")
sink()
"""


def find_rscript(rscript_arg):
    if rscript_arg is None or str(rscript_arg).strip() == "":
        rscript_arg = "Rscript"

    rscript_arg = str(rscript_arg)

    if os.path.isabs(rscript_arg):
        if os.path.exists(rscript_arg) and os.access(rscript_arg, os.X_OK):
            return rscript_arg
        raise FileNotFoundError(f"Rscript path was provided but is not executable: {rscript_arg}")

    found = shutil.which(rscript_arg)
    if found:
        return found

    raise FileNotFoundError(
        "Cannot find Rscript in PATH.\n"
        "Load R first, for example:\n"
        "  module load R/4.1.3\n"
        "  export R_LIBS_USER=~/Rlibs\n"
        "Or pass --rscript /path/to/Rscript.\n"
        f"Current PATH:\n{os.environ.get('PATH', '')}"
    )


def write_distance_and_metadata_for_r(D, sample_ids, md, dm_file, meta_file):
    dm_df = pd.DataFrame(D, index=sample_ids, columns=sample_ids)
    dm_df.to_csv(dm_file, sep="\t")

    r_md = md.set_index("sample_id_for_analysis").loc[sample_ids].reset_index().copy()

    keep_cols = [
        "sample_id_for_analysis",
        "bioproject_for_analysis",
        "site_for_analysis",
        "pair_id_for_analysis",
        "pair_block_for_permutation",
        "status_for_analysis",
        "bp_status_label",
        "status_group_for_order",
        "bp_plot_order",
    ]
    keep_cols = [c for c in keep_cols if c in r_md.columns]

    r_md[keep_cols].to_csv(meta_file, sep="\t", index=False)


def run_adonis2_r(
    rscript,
    r_backend_script,
    D,
    sample_ids,
    md,
    formula_rhs,
    permutations,
    seed,
    block_col,
    by_arg,
    work_dir,
    result_stem,
):
    dm_file = work_dir / f"{result_stem}.distance_matrix.tsv"
    meta_file = work_dir / f"{result_stem}.metadata.tsv"
    out_file = work_dir / f"{result_stem}.adonis2.tsv"
    log_file = work_dir / f"{result_stem}.R.log"

    write_distance_and_metadata_for_r(D, sample_ids, md, dm_file, meta_file)

    cmd = [
        rscript,
        str(r_backend_script),
        str(dm_file),
        str(meta_file),
        str(out_file),
        str(formula_rhs),
        str(int(permutations)),
        str(int(seed)),
        str(block_col),
        str(by_arg),
        str(log_file),
    ]

    completed = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "R adonis2 failed.\n"
            f"Command:\n  {' '.join(cmd)}\n"
            f"STDOUT:\n{completed.stdout}\n"
            f"STDERR:\n{completed.stderr}\n"
            f"R log file:\n{log_file}"
        )

    if not out_file.exists():
        raise RuntimeError(f"R adonis2 did not create output file: {out_file}")

    ad = pd.read_csv(out_file, sep="\t")
    ad["r_adonis2_file"] = str(out_file)
    ad["r_log_file"] = str(log_file)

    return ad


def extract_term_from_adonis(ad, term):
    if ad.empty:
        return {}

    sub = ad.loc[ad["term"].astype(str).eq(str(term)), :].copy()

    if sub.empty:
        return {
            "adonis2_term_found": False,
            "adonis2_note": f"term_not_found:{term}",
        }

    row = sub.iloc[0].to_dict()

    def get_num(*names):
        for n in names:
            if n in row:
                return pd.to_numeric(pd.Series([row[n]]), errors="coerce").iloc[0]
        return np.nan

    p = get_num("Pr(>F)", "Pr..F.", "p.value", "p_value")
    Fv = get_num("F", "F.Model")
    R2 = get_num("R2")
    ss = get_num("SumOfSqs", "SumsOfSqs", "SumsOfSqs.")

    return {
        "adonis2_term_found": True,
        "adonis2_tested_term": term,
        "df_term": get_num("Df"),
        "sum_of_squares": ss,
        "r2_term": R2,
        "pseudo_F": Fv,
        "p_value": p,
        "p_significance_label": p_to_stars(p),
        "adonis2_note": "ok",
        "r_adonis2_file": row.get("r_adonis2_file", ""),
        "r_log_file": row.get("r_log_file", ""),
    }


# -----------------------------
# Plotting
# -----------------------------

def plot_pcoa_single(coords, pc1_pct, pc2_pct, stats_row, title, out_png, pooled=False):
    site_colors = {
        "Oral": "#D65F4B",
        "Gut": "#63B6C7",
    }

    fig, ax = plt.subplots(figsize=(6.2, 5.4))

    for site in ["Oral", "Gut"]:
        sub = coords.loc[coords["site_for_analysis"].eq(site), :].copy()
        if sub.empty:
            continue

        ax.scatter(
            sub["PC1"],
            sub["PC2"],
            s=28 if not pooled else 14,
            alpha=0.75 if not pooled else 0.45,
            c=site_colors.get(site, "gray"),
            label=site,
            edgecolors="none",
        )

    ax.axhline(0, color="lightgray", linewidth=0.7)
    ax.axvline(0, color="lightgray", linewidth=0.7)

    pc1_label = "NA" if not np.isfinite(pc1_pct) else f"{pc1_pct:.2f}"
    pc2_label = "NA" if not np.isfinite(pc2_pct) else f"{pc2_pct:.2f}"

    ax.set_xlabel(f"PC1 ({pc1_label}%)")
    ax.set_ylabel(f"PC2 ({pc2_label}%)")
    ax.set_title(title, fontsize=11)

    p_text = format_p(stats_row.get("p_value", np.nan))
    r2 = stats_row.get("r2_term", np.nan)
    Fv = stats_row.get("pseudo_F", np.nan)

    if pooled:
        formula = "adonis2: D ~ BioProject + site"
        detail = "pair-blocked permutation"
    else:
        formula = "adonis2: D ~ site"
        detail = "pair-blocked permutation"

    if np.isfinite(r2) and np.isfinite(Fv):
        annot = (
            f"{formula}\n"
            f"{detail}\n"
            f"R²(site)={r2:.4f}, F={Fv:.3f}, P={p_text}"
        )
    else:
        note = stats_row.get("adonis2_note", "")
        annot = f"{formula}\n{detail}\nP={p_text}\n{note}"

    ax.text(
        0.02,
        0.02,
        annot,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.85),
    )

    ax.legend(frameon=False, loc="best")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_pcoa_grid(coord_df, stat_df, metric, out_png, n_cols=4):
    d = coord_df.loc[coord_df["metric"] == metric, :].copy()
    if d.empty:
        return

    bp_info = (
        d[["bioproject", "bp_plot_order", "bp_status_label"]]
        .drop_duplicates()
        .sort_values(["bp_plot_order", "bioproject"])
        .copy()
    )

    bps = bp_info["bioproject"].tolist()

    n_bp = len(bps)
    n_rows = int(np.ceil(n_bp / n_cols))

    site_colors = {
        "Oral": "#D65F4B",
        "Gut": "#63B6C7",
    }

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * 4.0, n_rows * 3.6),
        squeeze=False,
    )

    for ax in axes.ravel():
        ax.axis("off")

    for idx, bp in enumerate(bps):
        ax = axes.ravel()[idx]
        ax.axis("on")

        sub = d.loc[d["bioproject"].eq(bp), :].copy()

        for site in ["Oral", "Gut"]:
            s = sub.loc[sub["site_for_analysis"].eq(site), :].copy()
            if s.empty:
                continue

            ax.scatter(
                s["PC1"],
                s["PC2"],
                s=16,
                alpha=0.70,
                c=site_colors.get(site, "gray"),
                label=site,
                edgecolors="none",
            )

        ax.axhline(0, color="lightgray", linewidth=0.6)
        ax.axvline(0, color="lightgray", linewidth=0.6)

        stat_sub = stat_df.loc[
            (stat_df["metric"] == metric)
            & (stat_df["bioproject"] == bp),
            :
        ].copy()

        if not stat_sub.empty:
            row = stat_sub.iloc[0]
            p_text = format_p(row.get("p_value", np.nan))
            r2 = float(row["r2_term"]) if pd.notna(row.get("r2_term", np.nan)) else np.nan
            title_extra = f"R²={r2:.3f}, P={p_text}" if np.isfinite(r2) else f"P={p_text}"
        else:
            title_extra = "adonis2 NA"

        label = sub["bp_status_label"].iloc[0] if "bp_status_label" in sub.columns and not sub.empty else bp

        ax.set_title(f"{label}\n{title_extra}", fontsize=8)
        ax.set_xlabel("PC1", fontsize=8)
        ax.set_ylabel("PC2", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.45)

    handles = [
        plt.Line2D([], [], marker="o", linestyle="None", color=site_colors["Oral"], label="Oral"),
        plt.Line2D([], [], marker="o", linestyle="None", color=site_colors["Gut"], label="Gut"),
    ]

    fig.legend(handles=handles, loc="upper right", frameon=False)

    fig.suptitle(
        f"{metric_display_name(metric)} paired-sample PCoA by BioProject",
        fontsize=13,
        y=0.995,
    )

    fig.tight_layout(rect=[0, 0, 0.98, 0.975])
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


# -----------------------------
# One analysis
# -----------------------------

def run_one_pcoa_adonis2(
    D,
    sample_ids,
    md,
    metric,
    scope,
    bioproject,
    permutations,
    seed,
    rscript,
    r_backend_script,
    adonis_dir,
    by_arg,
):
    md2 = md.set_index("sample_id_for_analysis").loc[sample_ids].reset_index()

    coords, pc1_pct, pc2_pct = run_pcoa(D, sample_ids)
    coords = coords.merge(md2, on="sample_id_for_analysis", how="left")

    coords["dataset_label"] = metric_to_dataset_label(metric)
    coords["metric"] = metric
    coords["analysis_scope"] = scope
    coords["bioproject"] = bioproject

    if scope == "per_bioproject":
        formula_rhs = "site"
        tested_term = "site"
        formula_label = "D ~ site"
        permutation_scheme = "site labels permuted within pair_id blocks"
        result_stem = f"per_bp__{clean_for_filename(metric)}__{clean_for_filename(bioproject)}"
    else:
        formula_rhs = "bioproject + site"
        tested_term = "site"
        formula_label = "D ~ bioproject + site"
        permutation_scheme = "site labels permuted within BioProject x pair_id blocks"
        result_stem = f"pooled__{clean_for_filename(metric)}"

    try:
        ad_full = run_adonis2_r(
            rscript=rscript,
            r_backend_script=r_backend_script,
            D=D,
            sample_ids=sample_ids,
            md=md2,
            formula_rhs=formula_rhs,
            permutations=permutations,
            seed=seed,
            block_col="pair_block_for_permutation",
            by_arg=by_arg,
            work_dir=adonis_dir,
            result_stem=result_stem,
        )
        term_stats = extract_term_from_adonis(ad_full, tested_term)
        term_stats["test_note"] = "ok"
    except Exception as e:
        term_stats = {
            "adonis2_term_found": False,
            "adonis2_tested_term": tested_term,
            "df_term": np.nan,
            "sum_of_squares": np.nan,
            "r2_term": np.nan,
            "pseudo_F": np.nan,
            "p_value": np.nan,
            "p_significance_label": "NA",
            "adonis2_note": f"adonis2_failed:{type(e).__name__}:{str(e)}",
            "test_note": "adonis2_failed",
            "r_adonis2_file": "",
            "r_log_file": "",
        }

    n_oral = int((md2["site_for_analysis"] == "Oral").sum())
    n_gut = int((md2["site_for_analysis"] == "Gut").sum())
    n_bp = int(md2["bioproject_for_analysis"].nunique())
    n_pairs = int(md2["pair_block_for_permutation"].nunique())

    stat_row = {
        "analysis_scope": scope,
        "dataset_label": metric_to_dataset_label(metric),
        "metric": metric,
        "bioproject": bioproject,
        "formula": formula_label,
        "tested_term": "site",
        "adonis2_by": by_arg,
        "permutation_scheme": permutation_scheme,
        "n_samples": int(len(sample_ids)),
        "n_bioprojects": n_bp,
        "n_pair_blocks": n_pairs,
        "n_oral_samples": n_oral,
        "n_gut_samples": n_gut,
        "pc1_percent": pc1_pct,
        "pc2_percent": pc2_pct,
        "permutations": int(permutations),
        "seed": int(seed),
    }
    stat_row.update(term_stats)

    var_row = {
        "analysis_scope": scope,
        "dataset_label": metric_to_dataset_label(metric),
        "metric": metric,
        "bioproject": bioproject,
        "n_samples": int(len(sample_ids)),
        "n_pair_blocks": n_pairs,
        "pc1_percent": pc1_pct,
        "pc2_percent": pc2_pct,
    }

    return coords, stat_row, var_row


# -----------------------------
# Main
# -----------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Paired-sample PCoA + R vegan::adonis2 PERMANOVA with pair_id-blocked permutation."
    )

    parser.add_argument("--distance-dir", default=str(DEFAULT_DISTANCE_DIR))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))

    parser.add_argument("--braycurtis-matrix", default=None)
    parser.add_argument("--weighted-unifrac-matrix", default=None)
    parser.add_argument("--unweighted-unifrac-matrix", default=None)

    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--bioproject-col", default="bioproject")
    parser.add_argument("--site-col", default="site_group_cleaned")
    parser.add_argument("--pair-col", default="pair_id")
    parser.add_argument("--status-col", default="status_cleaned")

    parser.add_argument("--metrics", default="braycurtis,weighted_unifrac,unweighted_unifrac")
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--adonis2-by", default="terms", choices=["terms", "margin", "onedf", "none"])

    parser.add_argument(
        "--paired",
        type=parse_bool,
        default=True,
        help="Use only complete paired samples. Default: TRUE. This script is designed for TRUE.",
    )
    parser.add_argument(
        "--strict-one-to-one",
        type=parse_bool,
        default=True,
        help="Require exactly one Oral and one Gut in each BioProject x pair_id. Default: TRUE.",
    )

    parser.add_argument("--rscript", default="Rscript", help="Rscript executable or full path. Default: Rscript")
    parser.add_argument("--grid-n-cols", type=int, default=4)
    parser.add_argument("--save-individual-bp-plots", action="store_true", default=True)
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    if not args.paired:
        raise ValueError("This paired adonis2 script is intended to run with --paired TRUE.")

    distance_dir = Path(args.distance_dir).resolve()
    metadata_fp = Path(args.metadata).resolve()
    out_dir = Path(args.out_dir).resolve()

    plot_dir = out_dir / "plots"
    pooled_plot_dir = plot_dir / "pooled"
    per_bp_plot_dir = plot_dir / "per_bioproject"
    grid_plot_dir = plot_dir / "grids"
    adonis_dir = out_dir / "adonis2_tables"
    r_tmp_dir = out_dir / "_r_backend"

    for d in [out_dir, plot_dir, pooled_plot_dir, per_bp_plot_dir, grid_plot_dir, adonis_dir, r_tmp_dir]:
        d.mkdir(parents=True, exist_ok=True)

    metrics = parse_list(args.metrics)
    supported_metrics = {"braycurtis", "weighted_unifrac", "unweighted_unifrac"}
    bad = sorted(set(metrics) - supported_metrics)
    if bad:
        raise ValueError(f"Unsupported metrics: {bad}. Supported: {sorted(supported_metrics)}")

    matrix_paths = {
        "braycurtis": Path(args.braycurtis_matrix).resolve() if args.braycurtis_matrix else distance_dir / "braycurtis_filter_ra_distance_matrix.tsv",
        "weighted_unifrac": Path(args.weighted_unifrac_matrix).resolve() if args.weighted_unifrac_matrix else distance_dir / "weighted_unifrac_filter_distance_matrix.tsv",
        "unweighted_unifrac": Path(args.unweighted_unifrac_matrix).resolve() if args.unweighted_unifrac_matrix else distance_dir / "unweighted_unifrac_filter_distance_matrix.tsv",
    }

    outputs = {
        "coords_per_bp": out_dir / "pcoa_coordinates_per_bp.tsv",
        "coords_pooled": out_dir / "pcoa_coordinates_pooled.tsv",
        "permanova_by_bp": out_dir / "pcoa_permanova_by_bp.tsv",
        "permanova_pooled": out_dir / "pcoa_permanova_pooled.tsv",
        "variance": out_dir / "pcoa_variance_explained.tsv",
        "bp_order": out_dir / "pcoa_bioproject_order.tsv",
        "paired_filter_summary": out_dir / "pcoa_paired_sample_filter_summary.tsv",
        "paired_metadata_used": out_dir / "pcoa_paired_metadata_used.tsv",
        "summary_json": out_dir / "pcoa_summary.json",
    }

    for p in outputs.values():
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  distance_dir: {distance_dir}", flush=True)
    print(f"  metadata:     {metadata_fp}", flush=True)
    print(f"  out_dir:      {out_dir}", flush=True)

    if not metadata_fp.exists():
        raise FileNotFoundError(f"Metadata not found: {metadata_fp}")

    for metric in metrics:
        if not matrix_paths[metric].exists():
            raise FileNotFoundError(f"Distance matrix for {metric} not found: {matrix_paths[metric]}")

    print("[0/7] Checking Rscript ...", flush=True)
    rscript = find_rscript(args.rscript)
    print(f"  Rscript: {rscript}", flush=True)
    print(f"  R_LIBS_USER: {os.environ.get('R_LIBS_USER', '')}", flush=True)

    r_backend_script = r_tmp_dir / "run_adonis2_pair_blocked.R"
    r_backend_script.write_text(R_ADONIS2_SCRIPT)

    print("[1/7] Reading metadata and filtering complete paired samples ...", flush=True)
    raw_md = read_table(metadata_fp)

    md, bp_status_labels, bp_order_df = prepare_metadata(
        metadata=raw_md,
        sample_col=args.sample_col,
        bioproject_col=args.bioproject_col,
        site_col=args.site_col,
        pair_col=args.pair_col,
        status_col=args.status_col,
    )

    paired_md, pair_filter_summary = filter_complete_pairs(
        md,
        strict_one_to_one=args.strict_one_to_one,
    )

    if paired_md.empty:
        raise RuntimeError("No complete paired samples remained after paired filtering.")

    bp_order_df = bp_order_df.loc[
        bp_order_df["bioproject"].astype(str).isin(set(paired_md["bioproject_for_analysis"].astype(str))),
        :
    ].copy()
    bp_order_df = bp_order_df.sort_values(["status_rank", "bioproject"]).copy()
    bp_order_df["bp_plot_order"] = np.arange(bp_order_df.shape[0], dtype=int)

    paired_md = paired_md.drop(columns=["bp_plot_order"], errors="ignore").merge(
        bp_order_df[["bioproject", "bp_plot_order"]],
        left_on="bioproject_for_analysis",
        right_on="bioproject",
        how="left",
    ).drop(columns=["bioproject"], errors="ignore")

    bp_order_df.to_csv(outputs["bp_order"], sep="\t", index=False)
    pair_filter_summary.to_csv(outputs["paired_filter_summary"], sep="\t", index=False)
    paired_md.to_csv(outputs["paired_metadata_used"], sep="\t", index=False)

    print(f"  paired samples kept: {paired_md['sample_id_for_analysis'].nunique()}", flush=True)
    print(f"  paired blocks kept:  {paired_md['pair_block_for_permutation'].nunique()}", flush=True)
    print(f"  BioProjects kept:    {paired_md['bioproject_for_analysis'].nunique()}", flush=True)

    coords_per_bp_rows = []
    coords_pooled_rows = []
    stats_per_bp_rows = []
    stats_pooled_rows = []
    var_rows = []
    plot_files = []
    metric_input_summaries = {}

    print("[2/7] Reading distance matrices and running paired PCoA/adonis2 ...", flush=True)

    for metric in metrics:
        print(f"  Metric: {metric}", flush=True)

        dm_df = read_distance_matrix(matrix_paths[metric])

        sample_ids_in_metadata = set(paired_md["sample_id_for_analysis"])
        sample_ids_used = [sid for sid in dm_df.index if sid in sample_ids_in_metadata]
        dm_df = dm_df.loc[sample_ids_used, sample_ids_used].copy()

        metric_input_summaries[metric] = {
            "distance_matrix": str(matrix_paths[metric]),
            "n_samples_in_distance_matrix_after_paired_filter": int(dm_df.shape[0]),
            "dataset_label": metric_to_dataset_label(metric),
        }

        # Pooled.
        print(f"    Pooled paired PCoA/adonis2: {metric}", flush=True)
        D, sample_ids, pooled_md = subset_dm_and_md(dm_df, paired_md, bp=None)

        if D.shape[0] >= 4 and pooled_md["site_for_analysis"].nunique() >= 2 and pooled_md["pair_block_for_permutation"].nunique() >= 2:
            coords, stat_row, var_row = run_one_pcoa_adonis2(
                D=D,
                sample_ids=sample_ids,
                md=pooled_md,
                metric=metric,
                scope="pooled_all_bioprojects",
                bioproject="ALL",
                permutations=args.permutations,
                seed=args.seed,
                rscript=rscript,
                r_backend_script=r_backend_script,
                adonis_dir=adonis_dir,
                by_arg=args.adonis2_by,
            )

            coords_pooled_rows.append(coords)
            stats_pooled_rows.append(stat_row)
            var_rows.append(var_row)

            out_png = pooled_plot_dir / f"pcoa_pooled_paired_adonis2_{clean_for_filename(metric)}.png"
            plot_pcoa_single(
                coords=coords,
                pc1_pct=var_row["pc1_percent"],
                pc2_pct=var_row["pc2_percent"],
                stats_row=stat_row,
                title=f"Pooled paired samples\n{metric_display_name(metric)} PCoA",
                out_png=out_png,
                pooled=True,
            )
            plot_files.append(str(out_png))
        else:
            print(f"    skipped pooled {metric}: too few samples, sites, or pair blocks.", flush=True)

        # Per BP.
        print(f"    Per-BioProject paired PCoA/adonis2: {metric}", flush=True)
        per_metric_coords_for_grid = []
        per_metric_stats_for_grid = []

        for bp in bp_order_df["bioproject"].astype(str).tolist():
            D_bp, sample_ids_bp, md_bp = subset_dm_and_md(dm_df, paired_md, bp=bp)

            if D_bp.shape[0] < 4:
                continue
            if md_bp["site_for_analysis"].nunique() < 2:
                continue
            if md_bp["pair_block_for_permutation"].nunique() < 2:
                continue

            try:
                coords_bp, stat_bp, var_bp = run_one_pcoa_adonis2(
                    D=D_bp,
                    sample_ids=sample_ids_bp,
                    md=md_bp,
                    metric=metric,
                    scope="per_bioproject",
                    bioproject=bp,
                    permutations=args.permutations,
                    seed=args.seed,
                    rscript=rscript,
                    r_backend_script=r_backend_script,
                    adonis_dir=adonis_dir,
                    by_arg=args.adonis2_by,
                )

                coords_per_bp_rows.append(coords_bp)
                stats_per_bp_rows.append(stat_bp)
                var_rows.append(var_bp)

                per_metric_coords_for_grid.append(coords_bp)
                per_metric_stats_for_grid.append(stat_bp)

                if args.save_individual_bp_plots:
                    metric_bp_plot_dir = per_bp_plot_dir / clean_for_filename(metric)
                    metric_bp_plot_dir.mkdir(parents=True, exist_ok=True)

                    out_png = metric_bp_plot_dir / f"pcoa_paired_adonis2_{clean_for_filename(bp)}.png"
                    plot_pcoa_single(
                        coords=coords_bp,
                        pc1_pct=var_bp["pc1_percent"],
                        pc2_pct=var_bp["pc2_percent"],
                        stats_row=stat_bp,
                        title=f"{bp_status_labels.get(bp, bp)}\n{metric_display_name(metric)} paired PCoA",
                        out_png=out_png,
                        pooled=False,
                    )
                    plot_files.append(str(out_png))

            except Exception as e:
                print(f"      skipped {bp} / {metric}: {type(e).__name__}: {e}", flush=True)
                continue

        if per_metric_coords_for_grid:
            grid_coords = pd.concat(per_metric_coords_for_grid, axis=0, ignore_index=True)
            grid_stats = pd.DataFrame(per_metric_stats_for_grid)

            out_grid_png = grid_plot_dir / f"pcoa_grid_per_bp_paired_adonis2_{clean_for_filename(metric)}.png"
            plot_pcoa_grid(
                coord_df=grid_coords,
                stat_df=grid_stats,
                metric=metric,
                out_png=out_grid_png,
                n_cols=args.grid_n_cols,
            )
            plot_files.append(str(out_grid_png))

    print("[3/7] Combining output tables ...", flush=True)

    coords_per_bp_df = pd.concat(coords_per_bp_rows, axis=0, ignore_index=True) if coords_per_bp_rows else pd.DataFrame()
    coords_pooled_df = pd.concat(coords_pooled_rows, axis=0, ignore_index=True) if coords_pooled_rows else pd.DataFrame()
    stats_per_bp_df = pd.DataFrame(stats_per_bp_rows)
    stats_pooled_df = pd.DataFrame(stats_pooled_rows)
    var_df = pd.DataFrame(var_rows)

    print("[4/7] Writing output tables ...", flush=True)

    coords_per_bp_df.to_csv(outputs["coords_per_bp"], sep="\t", index=False)
    coords_pooled_df.to_csv(outputs["coords_pooled"], sep="\t", index=False)
    stats_per_bp_df.to_csv(outputs["permanova_by_bp"], sep="\t", index=False)
    stats_pooled_df.to_csv(outputs["permanova_pooled"], sep="\t", index=False)
    var_df.to_csv(outputs["variance"], sep="\t", index=False)

    print("[5/7] Writing summary JSON ...", flush=True)

    summary = {
        "inputs": {
            "distance_dir": os.path.abspath(str(distance_dir)),
            "metadata": os.path.abspath(str(metadata_fp)),
            "distance_matrices": {m: str(matrix_paths[m]) for m in metrics},
            "rscript": rscript,
            "r_libs_user": os.environ.get("R_LIBS_USER", ""),
            "r_backend_script": str(r_backend_script),
        },
        "out_dir": os.path.abspath(str(out_dir)),
        "paired_filter": {
            "paired": True,
            "strict_one_to_one": bool(args.strict_one_to_one),
            "sample_col": args.sample_col,
            "bioproject_col": args.bioproject_col,
            "site_col": args.site_col,
            "pair_col": args.pair_col,
            "n_paired_samples_kept": int(paired_md["sample_id_for_analysis"].nunique()),
            "n_pair_blocks_kept": int(paired_md["pair_block_for_permutation"].nunique()),
            "n_bioprojects_kept": int(paired_md["bioproject_for_analysis"].nunique()),
        },
        "metrics": metrics,
        "metric_input_summaries": metric_input_summaries,
        "permutations": int(args.permutations),
        "seed": int(args.seed),
        "permanova_backend": "R vegan::adonis2",
        "adonis2_by": args.adonis2_by,
        "permanova_design": {
            "per_bioproject": "D ~ site",
            "pooled": "D ~ bioproject + site",
            "pooled_note": "site term is tested after BioProject when adonis2 by='terms'",
            "permutation": "restricted within BioProject x pair_id blocks; with one Oral and one Gut per block this is equivalent to independent swap/no-swap within each pair",
        },
        "outputs": {k: str(v) for k, v in outputs.items()},
        "adonis2_tables_dir": str(adonis_dir),
        "plots": plot_files,
        "n_per_bp_coordinate_rows": int(coords_per_bp_df.shape[0]),
        "n_pooled_coordinate_rows": int(coords_pooled_df.shape[0]),
        "n_per_bp_permanova_rows": int(stats_per_bp_df.shape[0]),
        "n_pooled_permanova_rows": int(stats_pooled_df.shape[0]),
        "pooled_permanova": json_safe_records(stats_pooled_df),
    }

    with open(outputs["summary_json"], "w") as f:
        json.dump(summary, f, indent=2)

    print("[6/7] Done.", flush=True)
    print("\nMain outputs:", flush=True)
    for k, v in outputs.items():
        print(f"  {k}: {v}", flush=True)

    print("\nPlot files:", flush=True)
    for p in plot_files:
        print(f"  {p}", flush=True)

    print("\n[7/7] Pooled PERMANOVA:", flush=True)
    if not stats_pooled_df.empty:
        cols = [
            "dataset_label",
            "metric",
            "formula",
            "tested_term",
            "n_samples",
            "n_bioprojects",
            "n_pair_blocks",
            "n_oral_samples",
            "n_gut_samples",
            "r2_term",
            "pseudo_F",
            "p_value",
            "p_significance_label",
            "test_note",
        ]
        cols = [c for c in cols if c in stats_pooled_df.columns]
        print(stats_pooled_df[cols].to_string(index=False), flush=True)
    else:
        print("  No pooled PERMANOVA rows.", flush=True)


if __name__ == "__main__":
    main()

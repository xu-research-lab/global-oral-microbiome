# Auto-extracted from disease_predict.ipynb (genus version) -- plotting cells
# Run AFTER run_rf_genus.py. Each cell is self-contained (reads csv files).
import matplotlib
matplotlib.use('Agg')
# ===== notebook cell 19 =====
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Patch
from sklearn.metrics import roc_auc_score

# =============================================================================
# PART 1 -- build the RF AUC table (kept in memory as rf_long)
# =============================================================================

# ------------------------------ CONFIG ---------------------------------------
RUN_INFO   = "run_info_genus.csv"          # columns: site, disease, study
RESULT_DIR = "result"                # where the prediction csv files live
RF_FILE    = os.path.join(RESULT_DIR, "RF_results_lso.csv")
OUTPUT_DIR = "result"                # where the (optional) summary csv files go
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Site tokens to drop from the whole analysis (as they appear in run_info).
EXCLUDE_SITES = {"oral_mucosa"}

# Map plot display name -> the site token used in filenames / run_info / groups.
SITE_TOKENS = {
    "Saliva":      "saliva",
    "Tooth":       "tooth",
    "Oral mucosa": "oral_mucosa",
}
TOKEN_TO_DISPLAY = {v: k for k, v in SITE_TOKENS.items()}

# ---- Ground-truth labels ----------------------------------------------------
META_PATH      = "meta_genus.tsv"
META_ID_COL    = "sample_id"   # sample id column; must match SampleID / sampel_id
META_LABEL_COL = "group"       # binary target column (case=1 / control=0)


def load_labels() -> pd.Series:
    """Return a Series: sample_id -> binary label (0/1)."""
    sep = "\t" if META_PATH.lower().endswith((".tsv", ".txt")) else ","
    meta = pd.read_csv(META_PATH, sep=sep)
    missing = [c for c in (META_ID_COL, META_LABEL_COL) if c not in meta.columns]
    if missing:
        raise KeyError(
            f"{missing} not found in {META_PATH}. "
            f"Available columns: {list(meta.columns)}"
        )
    labels = meta[META_LABEL_COL]
    labels.index = meta[META_ID_COL].astype(str)
    return labels[~labels.index.duplicated(keep="first")]


def safe_auc(y_true: pd.Series, y_score: pd.Series) -> float:
    """AUC that returns NaN instead of raising on degenerate inputs."""
    df = pd.DataFrame({"y": y_true.values, "s": y_score.values}).dropna()
    if len(df) < 2 or df["y"].nunique() < 2:
        return np.nan
    return roc_auc_score(df["y"], df["s"])


# ---- Load inputs ------------------------------------------------------------
run_info = pd.read_csv(RUN_INFO)
run_info.columns = [c.strip().lower() for c in run_info.columns]   # site/disease/study
run_info = run_info[~run_info["site"].isin(EXCLUDE_SITES)].copy()  # drop excluded sites
labels = load_labels()

# ---- Random Forest: a single file, split by the `group` column --------------
# group looks like "<site>_<disease>_<study>". We rebuild the expected group
# from run_info instead of splitting the string (disease names have underscores).
rf_raw = pd.read_csv(RF_FILE)                       # sampel_id, group, pred_prob
rf_raw["group"] = rf_raw["group"].astype(str)

rf_records = []
for _, r in run_info.iterrows():
    site, disease, study = str(r["site"]), str(r["disease"]), str(r["study"])
    group = f"{site}_{disease}_{study}"
    sub = rf_raw[rf_raw["group"] == group]
    if sub.empty:
        print(f"[RF] no rows for group, skipped: {group}")
        continue
    y_score = sub.set_index(sub["sampel_id"].astype(str))["pred_prob"]
    y_true = labels.reindex(y_score.index)
    rf_records.append(
        {"site": site, "disease": disease, "test_study": study,
         "auc": safe_auc(y_true, y_score)}
    )
rf_long = pd.DataFrame(rf_records)

# ---- (Optional) also write per-site summary tables to disk ------------------
cols = ["disease", "test_study", "auc"]
for display, token in SITE_TOKENS.items():
    if token in EXCLUDE_SITES:
        continue
    (rf_long[rf_long["site"] == token][cols]
        .to_csv(os.path.join(OUTPUT_DIR, f"rf_loso_auc_all_diseases_genus_{token}.csv"),
                index=False))
    print(f"[write] {token}: RF={len(rf_long[rf_long['site'] == token])} rows")


# =============================================================================
# PART 2 -- plot (single row), RF only
# =============================================================================

SITE_COLORS = {'Saliva': '#2166AC',       # blue
               'Tooth': '#B2182B',        # red
               'Oral mucosa': '#1B7837'}  # green

# ---- Model styling: only RF is drawn ----------------------------------------
MODELS        = ['RF']
MODEL_MARKERS = {'RF': 'o'}              # filled circle
MODEL_FILL    = {'RF': True}             # filled marker in the site colour
ANCHOR_MODELS = ['RF']
POINT_ALPHA = 0.5
# Layout (in data units)
SITE_SLOT = 0.75      # width reserved for one site inside a disease group
GROUP_GAP = 0.40      # empty space between two disease groups
MEAN_HALFWIDTH = 0.10

# Single model -> centred inside the site slot
MODEL_OFFSETS = {'RF': 0.0}

# Single-row layout
N_ROWS = 1               # number of stacked rows of subplots
FIG_ROW_HEIGHT = 3.4     # figure inches per row
INCHES_PER_UNIT = 0.95   # figure inches per data unit of x-axis
MIN_WIDTH, MAX_WIDTH = 8.0, 22.0

MARK_BEST_SITE = False     # star above the best site per disease (RF mean)

# Disease display name -> abbreviation shown on the x axis (fallback: full name)
DISEASE_ABBR = {
    "Periodontal_Disease":               "PD",
    "Dental_Caries":                     "DC",
    "Gastroesophageal_Reflux_Disease":   "GERD",
    "HIV":                               "HIV",
    "Head_and_Neck_Cancer":              "HNC",
    "OSCC":                              "OSCC",
    "Lung_Cancer":                       "LC",
    "Obesity":                           "Obesity",
    "Colorectal_Cancer":                 "CRC",
    "Hypertension":                      "HTN",
    "Thyroid_Cancer":                    "TC",
    "Breast_Cancer":                     "BC",
    "IBD":                               "IBD",
}

# ---------------------------------- BUILD LONG TABLE FROM MEMORY -------------
rf = rf_long.dropna(subset=['auc']).copy()
rf['model'] = 'RF'

data = rf.copy()
data['site'] = data['site'].map(TOKEN_TO_DISPLAY)   # token -> display name
data = data.dropna(subset=['site']).copy()          # drop unmapped / excluded sites

# Keep only diseases present for the anchor model (RF) within each site
kept = []
for site_display, g in data.groupby('site'):
    shared = set(g[g.model == 'RF']['disease'])
    kept.append(g[g['disease'].isin(shared)])
data = pd.concat(kept, ignore_index=True)

# Display sites actually present, in the configured order
SITE_ORDER = [d for d in SITE_TOKENS if d in set(data['site'])]

# ----------------------------------------------------------- DISEASE ORDER ---
rf_all = data[data.model == 'RF']
n_sites = rf_all.groupby('disease')['site'].nunique()
med_auc = rf_all.groupby('disease')['auc'].median()
order = (pd.DataFrame({'n_sites': n_sites, 'med': med_auc})
         .sort_values(['n_sites', 'med'], ascending=[False, False])
         .index.tolist())

# Which sites are actually present for each disease
present_sites = {
    d: [s for s in SITE_ORDER
        if not data[(data.disease == d) & (data.site == s)].empty]
    for d in order
}

# ---------------------------------- SPLIT DISEASES INTO N_ROWS BALANCED ROWS -
def disease_width(d):
    """Data-unit width occupied by one disease group (incl. trailing gap)."""
    return len(present_sites[d]) * SITE_SLOT + GROUP_GAP


def split_contiguous(items, n):
    """Split an ordered list into n contiguous groups with balanced total width."""
    widths = [disease_width(d) for d in items]
    target = sum(widths) / n if n else sum(widths)
    rows, cur, cur_w = [], [], 0.0
    for i, (d, w) in enumerate(zip(items, widths)):
        cur.append(d)
        cur_w += w
        rows_left_after = n - len(rows) - 1
        items_left = len(items) - i - 1
        # close the row once we reach the target, keeping >=1 item per remaining row
        if len(rows) < n - 1 and cur_w >= target and items_left > rows_left_after:
            rows.append(cur)
            cur, cur_w = [], 0.0
    if cur:
        rows.append(cur)
    while len(rows) < n:          # pad with empty rows if needed
        rows.append([])
    return rows


rows = split_contiguous(order, N_ROWS)

# Row width in data units (drop the trailing GROUP_GAP of the last group)
def row_units(row_order):
    return sum(disease_width(d) for d in row_order) - GROUP_GAP if row_order else 0.0

max_units = max((row_units(r) for r in rows), default=1.0)

# Tick font sized to the busiest row
max_per_row = max((len(r) for r in rows), default=1)
tick_fs = float(np.clip(140.0 / max(max_per_row, 1), 8, 12))


# ------------------------------------------------------- RENDER ONE ROW ------
def render_row(ax, row_order):
    site_x, tick_pos, band = {}, [], []
    cursor = 0.0
    for d in row_order:
        k = len(present_sites[d])
        width = k * SITE_SLOT
        for j, s in enumerate(present_sites[d]):
            site_x[(d, s)] = cursor + (j + 0.5) * SITE_SLOT
        tick_pos.append(cursor + width / 2)
        band.append((cursor - GROUP_GAP / 2, cursor + width + GROUP_GAP / 2))
        cursor += width + GROUP_GAP

    for d in row_order:
        for s in present_sites[d]:
            d_site = data[(data.disease == d) & (data.site == s)]
            xc = site_x[(d, s)]

            # Scatter + mean bar for RF
            vals = d_site[d_site.model == 'RF'].set_index('test_study')['auc']
            if len(vals) == 0:
                continue
            x = xc + MODEL_OFFSETS['RF']
            face = SITE_COLORS[s] if MODEL_FILL['RF'] else 'white'
            ax.scatter(np.full(len(vals), x), vals.values,
                       marker=MODEL_MARKERS['RF'], s=44,
                       facecolors=face, edgecolors=SITE_COLORS[s],
                       linewidths=1.1, alpha=POINT_ALPHA, zorder=3)
            ax.hlines(vals.mean(), x - MEAN_HALFWIDTH, x + MEAN_HALFWIDTH,
                      colors=SITE_COLORS[s], linewidths=2.6, zorder=4)

        # Star the best-performing site for this disease (RF mean AUC)
        if MARK_BEST_SITE and len(present_sites[d]) >= 2:
            means = (data[(data.disease == d) & (data.model == 'RF')]
                     .groupby('site')['auc'].mean())
            if not means.empty:
                best = means.idxmax()
                ax.text(site_x[(d, best)], 1.005, '*',
                        ha='center', va='bottom', fontsize=16,
                        color=SITE_COLORS[best], zorder=5)

    ax.axhline(0.5, ls='--', c='grey', lw=1, zorder=2)
    ax.set_ylabel('LOSO AUC')
    ax.set_ylim(0, 1.06)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xlim(-GROUP_GAP / 2, max_units + GROUP_GAP / 2)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels([DISEASE_ABBR.get(d, d) for d in row_order],
                       rotation=35, ha='right', rotation_mode='anchor',
                       fontsize=tick_fs)
    ax.spines[['top', 'right']].set_visible(False)


# ------------------------------------------------------------------ FIGURE ---
fig, ax = plt.subplots(1, 1, figsize=(5.5, 2.45))

render_row(ax, rows[0])

# ---------------------------------------------- LEGEND -----
handles = [Patch(facecolor=SITE_COLORS[s], edgecolor='none', label=s) for s in SITE_ORDER]
if MARK_BEST_SITE:
    handles.append(mlines.Line2D([], [], color='black', marker='*', linestyle='None',
                                 markersize=9, label='Best site (RF mean)'))

ax.legend(handles=handles, loc='upper center', ncol=len(SITE_ORDER), fontsize=9,
          frameon=True, framealpha=0.9, edgecolor='#DDDDDD',
          bbox_to_anchor=(0.5, -0.38))

fig.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "loso_rf_1row_genus.png"),
             dpi=200, bbox_inches="tight")
plt.savefig(os.path.join(OUTPUT_DIR, "loso_rf_1row_genus.pdf"),
             bbox_inches="tight")
plt.show()

# ===== notebook cell 34 =====
"""
Five-fold CV, paired Oral vs Gut samples: RF only.

PART 1 builds df_all_rf in memory from the raw prediction files:
    run_info_oral_gut_genus.csv                      -> site / fold / study
    result/RF_results_five_fold.csv            -> RF (pred_prob)
    meta_oral_gut_genus.tsv                          -> TRUE case/control labels

One AUC is computed per (site, study, fold); the 5 folds form the box.
PART 2 plots, per disease/study, the Gut vs Oral boxes and connects each
fold's Gut point to its matching Oral point with a line.
"""

import math
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score

# =============================================================================
# PART 1 -- build df_all_rf  (RF only)
# =============================================================================

# ------------------------------------------------------------------ CONFIG ---
RUN_INFO   = "run_info_oral_gut_genus.csv"        # columns: site, fold, study
RESULT_DIR = "result"

# ---- Ground-truth labels ----------------------------------------------------
META_PATH      = "meta_oral_gut_genus.tsv"
META_ID_COL     = "sample_id"        # sample id column; must match SampleID
META_LABEL_COL  = "group"            # binary target column (case=1 / control=0)
META_STATUS_COL = "status_cleaned"   # disease / HC label column, used for titles


def load_labels() -> pd.Series:
    """Return a Series: sample_id -> binary label (0/1)."""
    sep = "\t" if META_PATH.lower().endswith((".tsv", ".txt")) else ","
    meta = pd.read_csv(META_PATH, sep=sep)
    missing = [c for c in (META_ID_COL, META_LABEL_COL) if c not in meta.columns]
    if missing:
        raise KeyError(
            f"{missing} not found in {META_PATH}. "
            f"Available columns: {list(meta.columns)}"
        )
    labels = meta[META_LABEL_COL]
    labels.index = meta[META_ID_COL].astype(str)
    return labels[~labels.index.duplicated(keep="first")]


def load_status() -> pd.Series:
    """Return a Series: sample_id -> status_cleaned (disease name or 'HC')."""
    sep = "\t" if META_PATH.lower().endswith((".tsv", ".txt")) else ","
    meta = pd.read_csv(META_PATH, sep=sep)
    missing = [c for c in (META_ID_COL, META_STATUS_COL) if c not in meta.columns]
    if missing:
        raise KeyError(
            f"{missing} not found in {META_PATH}. "
            f"Available columns: {list(meta.columns)}"
        )
    status = meta[META_STATUS_COL].astype(str)
    status.index = meta[META_ID_COL].astype(str)
    return status[~status.index.duplicated(keep="first")]


def safe_auc(y_true: pd.Series, y_score: pd.Series) -> float:
    """AUC that returns NaN instead of raising on degenerate inputs."""
    df = pd.DataFrame({"y": y_true.values, "s": y_score.values}).dropna()
    if len(df) < 2 or df["y"].nunique() < 2:
        return np.nan
    return roc_auc_score(df["y"], df["s"])


def _fold_str(v) -> str:
    """Fold as a bare integer string ('1'), tolerant of int/float/str input."""
    try:
        return str(int(v))
    except (ValueError, TypeError):
        return str(v).strip()


# ---- Load inputs ------------------------------------------------------------
run_info = pd.read_csv(RUN_INFO)
run_info.columns = [c.strip().lower() for c in run_info.columns]   # site / fold / study
labels = load_labels()
status = load_status()

# ---- RF: one row-group per (site, study, fold) ------------------------------
rf_raw = pd.read_csv("result/RF_results_five_fold.csv")            # sampel_id, group, pred_prob
rf_raw["group"] = rf_raw["group"].astype(str)

rf_records = []
for _, r in run_info.iterrows():
    site, study, fold = str(r["site"]).strip(), str(r["study"]).strip(), _fold_str(r["fold"])
    group = f"{site}_{study}_{fold}"
    sub = rf_raw[rf_raw["group"] == group]
    if sub.empty:
        print(f"[RF] no rows for group, skipped: {group}")
        continue
    sids = sub["sampel_id"].astype(str)
    y_score = sub.set_index(sids)["pred_prob"]
    y_true = labels.reindex(y_score.index)

    # Disease name for the title: look up this group's samples in the metadata
    # by sample id, then take the non-HC value from status_cleaned.
    grp_status = status.reindex(y_score.index).dropna()
    non_hc = grp_status[grp_status.str.upper() != "HC"]
    disease_name = non_hc.mode().iloc[0] if not non_hc.empty else study

    rf_records.append(
        {"site": site, "study": study, "fold": fold,
         "disease": disease_name, "auc": safe_auc(y_true, y_score)}
    )
df_all_rf = pd.DataFrame(rf_records)

# Drop rows where AUC could not be computed (single-class / no-label folds).
df_all_rf = df_all_rf.dropna(subset=["auc"]).reset_index(drop=True)

assert not df_all_rf.empty, (
    "No AUC values were computed. Check META_ID_COL / META_LABEL_COL and that "
    "sample ids in the prediction files match the metadata."
)
print(f"[built] RF rows={len(df_all_rf)}")


# =============================================================================
# PART 2 -- plot: per-disease Gut vs Oral boxes, folds connected by lines
# =============================================================================

# x = study, boxes grouped by site (Gut / Oral), RF only.
# For each fold, the Gut point is connected to its matching Oral point.
# Points sit exactly at each box center (no horizontal jitter/spread).

# -- INPUT --
combined = df_all_rf.copy()
combined["model"] = "RF"

# -- Global style --
sns.set_theme(style="whitegrid", context="notebook")
mpl.rcParams.update({
    "font.family":     "DejaVu Sans",
    "axes.titlesize":  13, "axes.titleweight": "bold",
    "axes.labelsize":  11, "axes.edgecolor":   "#666666",
    "axes.linewidth":  0.8,
    "grid.color":      "#E6E6E6", "grid.linewidth": 0.8,
    "xtick.labelsize": 10, "ytick.labelsize":  9,
})

site_order = ["Gut", "Oral"]
palette = {
    "Gut":  "#0B2D6B",
    "Oral": "#9C3D00",
}

# -- Spacing parameters --
BOX_WIDTH = 0.55
BOX_GAP   = 0.15
N_HUES    = len(site_order)             # 2
SLOT_W    = BOX_WIDTH / N_HUES          # width allocated per hue group


def box_center_offset(h_idx: int) -> float:
    """Offset of hue box center from the category tick (matches seaborn dodge)."""
    return -BOX_WIDTH / 2 + (h_idx + 0.5) * SLOT_W


gut_off  = box_center_offset(site_order.index("Gut"))
oral_off = box_center_offset(site_order.index("Oral"))

diseases = combined["disease"].unique()
n = len(diseases); ncols = 2; nrows = math.ceil(n / ncols)
fig, axes = plt.subplots(nrows, ncols,
                         figsize=(ncols * 3, nrows * 2),
                         squeeze=False)
axes = axes.flatten()

handles, labels = None, None

for i, disease in enumerate(diseases):
    ax = axes[i]
    sub = combined[combined["disease"] == disease]

    # Explicit study order so boxplot and scatter share the same x mapping
    study_order = list(sub["study"].unique())
    study_to_x  = {s: j for j, s in enumerate(study_order)}

    # -- Box plot --
    sns.boxplot(
        data=sub, x="study", y="auc",
        hue="site", hue_order=site_order,
        order=study_order,          # anchor x-axis order
        ax=ax,
        width=BOX_WIDTH,
        gap=BOX_GAP,
        showfliers=False,
        palette=palette,
        linewidth=1.0,
        boxprops=dict(alpha=0.8, edgecolor="#3a3a3a"),
        medianprops=dict(color="#222222", linewidth=1.5),
        whiskerprops=dict(color="#666666", linewidth=1.0),
        capprops=dict(color="#666666", linewidth=1.0),
    )

    # -- Per-fold paired points + connecting lines (no jitter) --
    for study in study_order:
        x_tick = study_to_x[study]
        s_sub  = sub[sub["study"] == study]

        gut  = s_sub[s_sub["site"] == "Gut"].groupby("fold")["auc"].mean()
        oral = s_sub[s_sub["site"] == "Oral"].groupby("fold")["auc"].mean()

        all_folds = sorted(set(gut.index) | set(oral.index), key=lambda v: int(v))

        # Points sit at the box center; same x for every fold (no spread).
        xg = x_tick + gut_off
        xo = x_tick + oral_off

        for fold in all_folds:
            has_g = fold in gut.index
            has_o = fold in oral.index

            # connecting line: Gut fold -> matching Oral fold
            if has_g and has_o:
                ax.plot([xg, xo], [gut[fold], oral[fold]],
                        color="#9A9A9A", lw=0.9, alpha=0.7, zorder=2)

            # scatter points
            if has_g:
                ax.scatter([xg], [gut[fold]], color=palette["Gut"],
                           s=3.5 ** 2, alpha=0.95,
                           edgecolors="white", linewidths=0.4, zorder=3)
            if has_o:
                ax.scatter([xo], [oral[fold]], color=palette["Oral"],
                           s=3.5 ** 2, alpha=0.95,
                           edgecolors="white", linewidths=0.4, zorder=3)

    ax.set_title(str(disease), pad=10)
    ax.set_xlabel(""); ax.set_ylabel("AUC")
    ax.grid(axis="x", visible=False)
    ax.axhline(0.5, color="#BBBBBB", ls="--", lw=0.8, zorder=0)
    sns.despine(ax=ax)

    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    if ax.get_legend() is not None:
        ax.get_legend().remove()

# -- Shared legend --
for j in range(n, len(axes)):
    axes[j].axis("off")

if n < len(axes):
    axes[n].legend(handles, labels, title="Site", loc="center",
                   frameon=True, fontsize=11, title_fontsize=12, markerscale=1.3)
else:
    fig.legend(handles, labels, title="Site", loc="upper right",
               frameon=True, fontsize=10, title_fontsize=11)

fig.tight_layout() 
plt.savefig(os.path.join(RESULT_DIR, "fivefold_oral_gut_rf_paired_genus.png"),
             dpi=200, bbox_inches="tight")
plt.show()

# ===== notebook cell 40 =====
"""
RF / DeepPhylo -- within-continent (wc, analysis A) vs leave-one-continent-out
(loco, analysis B), drawn in ONE panel for easy comparison.

*** Changes in this version ***
1) Each loco fold now holds out one continent and tests on a SINGLE study from
   that continent, so wc and loco share one naming scheme:
       <split>_<site>_<disease>_<continent>_<test_study>
2) All mean markers are removed (the former continent-mean triangles are gone).
   Every point on the panel is the AUC of one test study.
3) In the "Both (wc vs loco)" section, each connector links the SAME test study
   under the two training setups (within-continent vs cross-continent), so the
   slope of a line is that study's change in generalization performance.

Layout: diseases are grouped left-to-right into three sections
    [ 1. BOTH wc + loco ]   |   [ 2. wc only ]   |   [ 3. loco only ]
Each disease occupies one x slot; when both splits exist, wc sits on the left
half of the slot and loco on the right half.
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Patch
import seaborn as sns
from sklearn.metrics import roc_auc_score

# =============================================================================
# PART 1 -- build the AUC table
# =============================================================================

# ------------------------------------------------------------------ CONFIG ---
RUN_INFO   = "run_info_continent_genus.csv"    # site, disease, continent, test_study, split_type
RESULT_DIR = "result"
RF_RESULTS = "result/RF_results_continent.csv"

# ---- Ground-truth labels ----------------------------------------------------
META_PATH      = "meta_genus.tsv"
META_ID_COL    = "sample_id"   # sample id column; must match SampleID (SRR...)
META_LABEL_COL = "group"       # binary target column (case=1 / control=0)

SPLIT_LABEL = {"wc": "Within-continent", "loco": "Cross-continent"}

DROP_SITES = ["oral_mucosa"]   # sites to exclude from the panel


def _to_binary(s: pd.Series) -> pd.Series:
    """Coerce a label column to 0/1 floats, accepting common string encodings."""
    if s.dtype != object:
        return pd.to_numeric(s, errors="coerce")
    mapping = {
        "case": 1, "control": 0, "disease": 1, "healthy": 0, "patient": 1,
        "positive": 1, "negative": 0, "yes": 1, "no": 0,
        "true": 1, "false": 0, "1": 1, "0": 0,
    }
    return s.astype(str).str.strip().str.lower().map(mapping)


def load_labels() -> pd.Series:
    """Return a Series: sample_id -> binary label (0/1)."""
    sep = "\t" if META_PATH.lower().endswith((".tsv", ".txt")) else ","
    meta = pd.read_csv(META_PATH, sep=sep)
    missing = [c for c in (META_ID_COL, META_LABEL_COL) if c not in meta.columns]
    if missing:
        raise KeyError(
            f"{missing} not found in {META_PATH}. "
            f"Available columns: {list(meta.columns)}"
        )
    labels = _to_binary(meta[META_LABEL_COL])
    labels.index = meta[META_ID_COL].astype(str)
    return labels[~labels.index.duplicated(keep="first")]


def safe_auc(y_true: pd.Series, y_score: pd.Series) -> float:
    """AUC that returns NaN instead of raising on degenerate inputs."""
    d = pd.DataFrame({"y": y_true.values, "s": y_score.values}).dropna()
    if len(d) < 2 or d["y"].nunique() < 2:
        return np.nan
    return roc_auc_score(d["y"], d["s"])


# ---- Load inputs ------------------------------------------------------------
run_info = pd.read_csv(RUN_INFO, dtype={"test_study": str})
run_info.columns = [c.strip().lower() for c in run_info.columns]
labels = load_labels()

rf_raw = pd.read_csv(RF_RESULTS)
rf_raw["group"] = rf_raw["group"].astype(str)

records = []
for _, r in run_info.iterrows():
    split     = str(r["split_type"]).strip()
    site      = str(r["site"]).strip()
    disease   = str(r["disease"]).strip()
    continent = str(r["continent"]).strip()
    study     = str(r["test_study"]).strip()

    # wc and loco now share one scheme: <split>_<site>_<disease>_<continent>_<test_study>
    group = f"{split}_{site}_{disease}_{continent}_{study}"

    sub = rf_raw[rf_raw["group"] == group]
    if sub.empty:
        print(f"[RF] no rows for group, skipped: {group}")
        continue

    sids    = sub["sampel_id"].astype(str)
    y_score = sub.set_index(sids)["pred_prob"]
    y_true  = labels.reindex(y_score.index)

    records.append({
        "site": site, "disease": disease, "continent": continent,
        "test_study": study, "split_type": split,
        "split": SPLIT_LABEL.get(split, split),
        "auc": safe_auc(y_true, y_score),
    })

df = pd.DataFrame(records).dropna(subset=["auc"]).reset_index(drop=True)
df = df.loc[~df.site.isin(DROP_SITES)].reset_index(drop=True)

assert not df.empty, (
    "No AUC values were computed. Check META_ID_COL / META_LABEL_COL, the group "
    "naming, and that sample ids in the prediction file match the metadata."
)
print(f"[built] {len(df)} AUC rows "
      f"(wc={(df.split_type=='wc').sum()}, loco={(df.split_type=='loco').sum()})")


# =============================================================================
# PART 2 -- single panel: one point per test study, no means.
#           [both wc+loco]  |  [wc only]  |  [loco only]
#           connectors link the SAME test study across wc / loco
# =============================================================================

sns.set_theme(style="whitegrid", context="notebook")
mpl.rcParams.update({
    "font.family":     "DejaVu Sans",
    "axes.titlesize":  13, "axes.titleweight": "bold",
    "axes.labelsize":  11, "axes.edgecolor":   "#666666",
    "axes.linewidth":  0.8,
    "grid.color":      "#E6E6E6", "grid.linewidth": 0.8,
    "xtick.labelsize": 10, "ytick.labelsize":  9,
})

# x categories: disease (prefixed with site only if more than one site present)
multi_site = df["site"].nunique() > 1
df["xcat"] = (df["site"] + " | " + df["disease"]) if multi_site else df["disease"]
df["test_study"] = df["test_study"].astype(str)

# which splits each disease actually has (may be wc-only or loco-only)
present = {x: set(df[df.xcat == x]["split_type"]) for x in df["xcat"].unique()}

# --- grouping / ordering -----------------------------------------------------
GROUP_BOTH, GROUP_WC, GROUP_LOCO = 0, 1, 2
GROUP_LABEL = {GROUP_BOTH: "Both  (wc vs loco)",
               GROUP_WC:   "wc only",
               GROUP_LOCO: "loco only"}
GROUP_TINT  = {GROUP_BOTH: "#EAF1FA", GROUP_WC: "#FBF3E7", GROUP_LOCO: "#F3F3F3"}


def group_rank(x):
    p = present[x]
    if p == {"wc", "loco"}:
        return GROUP_BOTH          # paired -> most informative, goes first
    if p == {"wc"}:
        return GROUP_WC
    return GROUP_LOCO              # loco only


def wc_level(x):
    """Ordering scalar for the paired section: mean AUC over all wc folds of
    this disease. Used only to sort the x axis; nothing is drawn for it."""
    v = df[(df.xcat == x) & (df.split_type == "wc")]["auc"]
    return v.mean() if len(v) else np.nan


def sort_key(x):
    g = group_rank(x)
    if g == GROUP_BOTH:
        return (g, -wc_level(x), x)
    return (g, -df[df.xcat == x]["auc"].mean(), x)


all_xcats = sorted(df["xcat"].unique(), key=sort_key)
pos = {x: i for i, x in enumerate(all_xcats)}
group_of = {x: group_rank(x) for x in all_xcats}

# contiguous section spans (group_rank, start_idx, end_idx) in slot coordinates
sections = []
if all_xcats:
    start = 0
    for i in range(1, len(all_xcats) + 1):
        if i == len(all_xcats) or group_of[all_xcats[i]] != group_of[all_xcats[start]]:
            sections.append((group_of[all_xcats[start]], start, i - 1))
            start = i

# consistent continent colours + within-cluster horizontal dodge
continents = sorted(df["continent"].unique())
K = len(continents)
pal = sns.color_palette("tab10", n_colors=max(K, 3))
cont_color = {c: pal[i] for i, c in enumerate(continents)}

HALF       = 0.25                            # when both: wc=pos-HALF, loco=pos+HALF
SUBSTEP    = min(0.40 / max(K, 1), 0.10)     # horizontal dodge between continents
STUDY_SPAN = SUBSTEP * 0.55                  # total spread of studies within a continent
cont_dodge = {c: (i - (K - 1) / 2) * SUBSTEP for i, c in enumerate(continents)}


def study_offsets(studies):
    """Small horizontal offset per study within a continent. wc and loco reuse
    the same offsets so that the connectors line up."""
    n = len(studies)
    if n <= 1:
        return {s: 0.0 for s in studies}
    step = STUDY_SPAN / (n - 1)
    return {s: (i - (n - 1) / 2) * step for i, s in enumerate(studies)}


def sub_centers(x):
    """Centre x of the wc / loco sub-clusters (a lone side is centred on the slot)."""
    i = pos[x]
    both = present[x] == {"wc", "loco"}
    return ((i - HALF) if both else i, (i + HALF) if both else i)


ncols = len(all_xcats)
fig, ax = plt.subplots(figsize=(15, 7))

# ---- section backgrounds + separators ---------------------------------------
for g, s0, s1 in sections:
    ax.axvspan(s0 - 0.5, s1 + 0.5, color=GROUP_TINT[g], zorder=0)
for _, _, s1 in sections[:-1]:
    ax.axvline(s1 + 0.5, color="#AAAAAA", ls="-", lw=1.4, zorder=1)
for x, i in pos.items():
    if i > 0:
        ax.axvline(i - 0.5, color="#DDDDDD", ls=":", lw=0.6, zorder=1)
    if present[x] == {"wc", "loco"}:
        ax.axvline(i, color="#C9C9C9", ls=":", lw=0.8, zorder=1)
for g, s0, s1 in sections:
    if g != GROUP_BOTH:
        continue
    for i in range(s0, s1 + 1):
        ax.axvline(i - 0.5, color="#8C8C8C", ls="--", lw=1.0, zorder=1)
    ax.axvline(s1 + 0.5, color="#8C8C8C", ls="--", lw=1.0, zorder=1)

# ---- data: one point per test study + same-study wc->loco connectors --------
for x in all_xcats:
    sub          = df[df["xcat"] == x]
    wc_c, loco_c = sub_centers(x)
    paired       = present[x] == {"wc", "loco"}

    for cont in continents:
        csub = sub[sub.continent == cont]
        if csub.empty:
            continue
        col     = cont_color[cont]
        base_d  = cont_dodge[cont]
        studies = sorted(csub["test_study"].unique())
        off     = study_offsets(studies)

        wc_map = (csub[csub.split_type == "wc"]
                  .set_index("test_study")["auc"].to_dict())
        lc_map = (csub[csub.split_type == "loco"]
                  .set_index("test_study")["auc"].to_dict())

        for st in studies:
            dx     = base_d + off[st]
            xw, xl = wc_c + dx, loco_c + dx
            a_wc   = wc_map.get(st)
            a_lc   = lc_map.get(st)

            # same test study: wc -> loco connector
            if paired and a_wc is not None and a_lc is not None:
                ax.plot([xw, xl], [a_wc, a_lc], color=col, lw=1.2, alpha=0.55,
                        zorder=2, solid_capstyle="round")

            if a_wc is not None:
                ax.scatter([xw], [a_wc], marker="o", s=42, facecolors=col,
                           edgecolors="white", linewidths=0.6, alpha=0.90, zorder=3)
            if a_lc is not None:
                ax.scatter([xl], [a_lc], marker="D", s=52, facecolors=col,
                           edgecolors="white", linewidths=0.8, alpha=0.95, zorder=4)

ax.axhline(0.5, color="#BBBBBB", ls="--", lw=0.8, zorder=1)
ax.set_ylim(0, 1.03)
ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.set_ylabel("AUC")
ax.set_xlim(-0.7, ncols - 0.3)

# ---- two-level x labels: wc/loco sub-labels, then 45-deg disease name -------
sub_ticks, sub_labels = [], []
for x in all_xcats:
    wc_c, loco_c = sub_centers(x)
    if present[x] == {"wc", "loco"}:
        sub_ticks  += [wc_c, loco_c]
        sub_labels += ["wc", "loco"]
    elif "wc" in present[x]:
        sub_ticks.append(wc_c)
        sub_labels.append("wc")
    else:
        sub_ticks.append(loco_c)
        sub_labels.append("loco")
ax.set_xticks(sub_ticks)
ax.set_xticklabels(sub_labels, fontsize=8, color="#555555")

trans = ax.get_xaxis_transform()            # x=data, y=axes-fraction
for x, i in pos.items():
    ax.text(i, -0.10, x.replace("_", " "), transform=trans,
            ha="right", va="top", rotation=45, fontsize=10, fontweight="bold")

# ---- section headers above the panel ---------------------------------------
for g, s0, s1 in sections:
    ax.text((s0 + s1) / 2, 1.045, GROUP_LABEL[g], transform=trans,
            ha="center", va="bottom", fontsize=11, fontweight="bold",
            color="#333333",
            bbox=dict(boxstyle="round,pad=0.3", fc=GROUP_TINT[g],
                      ec="#BBBBBB", lw=0.8))

ax.grid(axis="x", visible=False)
sns.despine(ax=ax)

# ---- legend: continents (colour) + marker meaning --------------------------
cont_handles = [Patch(facecolor=cont_color[c], edgecolor="none",
                      label=c.replace("_", " ")) for c in continents]
marker_handles = [
    mlines.Line2D([], [], color="#555555", marker="o", linestyle="None",
                  markersize=8, markeredgecolor="white",
                  label="wc: test study (trained within continent)"),
    mlines.Line2D([], [], color="#555555", marker="D", linestyle="None",
                  markersize=9, markeredgecolor="white",
                  label="loco: test study (its continent held out)"),
    mlines.Line2D([], [], color="#999999", linestyle="-", linewidth=1.4,
                  alpha=0.6, label="same test study (wc -> loco)"),
]
leg1 = ax.legend(handles=cont_handles, title="Continent",
                 loc="upper left", bbox_to_anchor=(1.015, 1.0),
                 bbox_transform=ax.transAxes,
                 frameon=True, fontsize=10, title_fontsize=11,
                 borderaxespad=0.0)
ax.add_artist(leg1)
ax.legend(handles=marker_handles, title="Marker",
          loc="lower left", bbox_to_anchor=(1.015, 0.0),
          bbox_transform=ax.transAxes,
          frameon=True, fontsize=10, title_fontsize=11,
          borderaxespad=0.0)

ax.set_title("Each line links one test study: evaluated within its own continent "
             "vs with its continent held out of training", pad=26)
fig.suptitle("RF - within- vs cross-continent generalization (per test study)",
             fontsize=14, fontweight="bold", y=1.06)

fig.tight_layout()
plt.savefig(os.path.join(RESULT_DIR, "continent_A_vs_B_combined_genus.png"),
             dpi=200, bbox_inches="tight")
plt.show()

# ===== notebook cell 45 =====
"""
DeepPhylo cross-disease transfer: train on disease A, test on disease B
(all ordered pairs, within each site).

PART 1 -- one AUC per run_info row from the prediction files:
    run_info_cross_disease_genus.csv
        columns include: site, train_disease, test_disease, split_type
    result/predictions_deepphylo_cd_<site>_train_<train_disease>_test_<test_disease>.csv
    meta_genus.tsv  -> TRUE case/control labels (SampleID = ERR...)

PART 2 -- one transfer-matrix heatmap per site:
    rows    = train disease (A)   -- only diseases that actually have results
    columns = test disease  (B)   -- only diseases that actually have results
    colour  = AUC (RdBu_r, 0.5 = white = chance; red better, blue worse)
"""

import math
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import seaborn as sns
from sklearn.metrics import roc_auc_score

# =============================================================================
# PART 1 -- build the AUC table
# =============================================================================

# ------------------------------------------------------------------ CONFIG ---
RUN_INFO   = "run_info_cross_disease_genus.csv"
RESULT_DIR = "result"

# ---- Ground-truth labels ----------------------------------------------------
META_PATH      = "meta_genus.tsv"
META_ID_COL    = "sample_id"   # sample id column; must match SampleID (ERR...)
META_LABEL_COL = "group"       # binary target column (case=1 / control=0)


def _to_binary(s: pd.Series) -> pd.Series:
    """Coerce a label column to 0/1 floats, accepting common string encodings."""
    if s.dtype != object:
        return pd.to_numeric(s, errors="coerce")
    mapping = {
        "case": 1, "control": 0, "disease": 1, "healthy": 0, "patient": 1,
        "positive": 1, "negative": 0, "yes": 1, "no": 0,
        "true": 1, "false": 0, "1": 1, "0": 0,
    }
    return s.astype(str).str.strip().str.lower().map(mapping)


def load_labels() -> pd.Series:
    """Return a Series: sample_id -> binary label (0/1)."""
    sep = "\t" if META_PATH.lower().endswith((".tsv", ".txt")) else ","
    meta = pd.read_csv(META_PATH, sep=sep)
    missing = [c for c in (META_ID_COL, META_LABEL_COL) if c not in meta.columns]
    if missing:
        raise KeyError(
            f"{missing} not found in {META_PATH}. "
            f"Available columns: {list(meta.columns)}"
        )
    labels = _to_binary(meta[META_LABEL_COL])
    labels.index = meta[META_ID_COL].astype(str)
    return labels[~labels.index.duplicated(keep="first")]


def safe_auc(y_true: pd.Series, y_score: pd.Series) -> float:
    """AUC that returns NaN instead of raising on degenerate inputs."""
    df = pd.DataFrame({"y": y_true.values, "s": y_score.values}).dropna()
    if len(df) < 2 or df["y"].nunique() < 2:
        return np.nan
    return roc_auc_score(df["y"], df["s"])


# ---- Load inputs ------------------------------------------------------------
run_info = pd.read_csv(RUN_INFO)
run_info.columns = [c.strip().lower() for c in run_info.columns]
labels = load_labels()

records = []
# for _, r in run_info.iterrows():
#     site = str(r["site"]).strip()
#     train = str(r["train_disease"]).strip()
#     test = str(r["test_disease"]).strip()
#     path = os.path.join(
#         RESULT_DIR,
#         f"predictions_deepphylo_cd_{site}_train_{train}_test_{test}.csv",
#     )
#     if not os.path.exists(path):
#         print(f"[DeepPhylo] missing, skipped: {path}")
#         continue
#     pred = pd.read_csv(path)                        # SampleID, Prob_1, Pred
#     y_score = pred.set_index(pred["SampleID"].astype(str))["Prob_1"]
#     y_true = labels.reindex(y_score.index)
#     records.append({
#         "site": site, "train_disease": train, "test_disease": test,
#         "auc": safe_auc(y_true, y_score),
#     })

rf_raw = pd.read_csv("result/RF_results_cross_disease.csv")
rf_raw["group"] = rf_raw["group"].astype(str)
records = []
for _, r in run_info.iterrows():
    site = str(r["site"]).strip()
    train = str(r["train_disease"]).strip()
    test = str(r["test_disease"]).strip()

    group = f"cd_{site}_train_{train}_test_{test}"
    sub = rf_raw[rf_raw["group"] == group]
    if sub.empty:
        print(f"[RF] no rows for group, skipped: {group}")
        continue
    sids = sub["sampel_id"].astype(str)
    y_score = sub.set_index(sids)["pred_prob"]
    y_true = labels.reindex(y_score.index)

    records.append({"site": site, "train_disease": train, "test_disease": test,
        "auc": safe_auc(y_true, y_score),
    })
df = pd.DataFrame(records)
assert not df.empty, (
    "No rows built. Check RESULT_DIR / file naming and META_ID_COL / META_LABEL_COL."
)
n_ok = df["auc"].notna().sum()
print(f"[built] {len(df)} train->test pairs ({n_ok} with a valid AUC)")

df.to_csv(os.path.join(RESULT_DIR, "cross_disease_auc_long_genus.csv"), index=False)


# =============================================================================
# PART 2 -- transfer-matrix heatmap(s), one per site
# =============================================================================

sns.set_theme(style="white", context="notebook")
mpl.rcParams.update({
    "font.family":    "DejaVu Sans",
    "axes.titlesize": 13, "axes.titleweight": "bold",
    "axes.labelsize": 11,
})

CMAP, VMIN, VMAX = "RdBu_r", 0.0, 1.0   # 0.5 -> white (chance)

valid = df[df["auc"].notna()]
sites = sorted(valid["site"].unique())

n = len(sites)
ncols = min(n, 2)
nrows = math.ceil(n / ncols)
fig, axes = plt.subplots(nrows, ncols,
                         figsize=(ncols * 7.2, nrows * 6.6),
                         squeeze=False, constrained_layout=True)
axes = axes.flatten()

for i, site in enumerate(sites):
    ax = axes[i]
    sub = valid[valid["site"] == site]

    # keep only diseases that really appear for this site (drop all-empty rows/cols)
    rows = sorted(sub["train_disease"].unique())
    cols = sorted(sub["test_disease"].unique())
    mat = (sub.pivot_table(index="train_disease", columns="test_disease", values="auc")
              .reindex(index=rows, columns=cols))

    sns.heatmap(
        mat, ax=ax, cmap=CMAP, vmin=VMIN, vmax=VMAX,
        annot=True, fmt=".2f", annot_kws={"size": 9},
        linewidths=0.6, linecolor="white",
        mask=mat.isna(), cbar=False, square=False,
    )
    ax.set_facecolor("#F0F0F0")                     # missing pairs within this site
    ax.set_title(site, pad=8)
    ax.set_xlabel("Test disease (B)")
    ax.set_ylabel("Train disease (A)")
    ax.set_xticklabels([t.get_text().replace("_", " ") for t in ax.get_xticklabels()],
                       rotation=40, ha="right", fontsize=9)
    ax.set_yticklabels([t.get_text().replace("_", " ") for t in ax.get_yticklabels()],
                       rotation=0, fontsize=9)

for j in range(n, len(axes)):
    axes[j].axis("off")

# Shared colour bar
sm = ScalarMappable(norm=Normalize(VMIN, VMAX), cmap=CMAP)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes.tolist(), fraction=0.03, pad=0.02, shrink=0.85)
cbar.set_label("AUC  (train A -> test B)")
cbar.ax.axhline(0.5, color="#333333", lw=1.0)       # mark the chance level

fig.suptitle("RF cross-disease transfer  -  train (rows) -> test (columns)",
             fontsize=14, fontweight="bold")
plt.savefig(os.path.join(RESULT_DIR, "cross_disease_transfer_genus.png"),
            dpi=200, bbox_inches="tight")
plt.show()

# ===== notebook cell 51 =====
"""
DeepPhylo learning curve: effect of dataset size (number of training studies)
on performance.

PART 1 -- one AUC per run_info row from the prediction files:
    run_info_lc_genus.csv
        columns: site, disease, n_studies, repeat, subset, test_study, split_type
    result/predictions_deepphylo_lc_<site>_<disease>_k<n_studies>_r<repeat>_test_<test_study>.csv
    meta_genus.tsv  -> TRUE case/control labels (SampleID = SRR...)

PART 2 -- AUC vs number of training studies:
    line   = mean AUC per disease at each dataset size
    band   = spread across repeats and held-out studies (std dev)
    dots   = individual (repeat, test_study) AUCs (faint)
    one panel per site.
"""

import math
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score

# =============================================================================
# PART 1 -- build the AUC table
# =============================================================================

# ------------------------------------------------------------------ CONFIG ---
RUN_INFO   = "run_info_lc_genus.csv"
RESULT_DIR = "result"

# ---- Ground-truth labels ----------------------------------------------------
META_PATH      = "meta_genus.tsv"
META_ID_COL    = "sample_id"   # sample id column; must match SampleID (SRR...)
META_LABEL_COL = "group"       # binary target column (case=1 / control=0)

ERRORBAR = "sd"   # band: "sd" (std dev), "se" (std error), or ("ci", 95)


def _to_binary(s: pd.Series) -> pd.Series:
    """Coerce a label column to 0/1 floats, accepting common string encodings."""
    if s.dtype != object:
        return pd.to_numeric(s, errors="coerce")
    mapping = {
        "case": 1, "control": 0, "disease": 1, "healthy": 0, "patient": 1,
        "positive": 1, "negative": 0, "yes": 1, "no": 0,
        "true": 1, "false": 0, "1": 1, "0": 0,
    }
    return s.astype(str).str.strip().str.lower().map(mapping)


def load_labels() -> pd.Series:
    """Return a Series: sample_id -> binary label (0/1)."""
    sep = "\t" if META_PATH.lower().endswith((".tsv", ".txt")) else ","
    meta = pd.read_csv(META_PATH, sep=sep)
    missing = [c for c in (META_ID_COL, META_LABEL_COL) if c not in meta.columns]
    if missing:
        raise KeyError(
            f"{missing} not found in {META_PATH}. "
            f"Available columns: {list(meta.columns)}"
        )
    labels = _to_binary(meta[META_LABEL_COL])
    labels.index = meta[META_ID_COL].astype(str)
    return labels[~labels.index.duplicated(keep="first")]


def safe_auc(y_true: pd.Series, y_score: pd.Series) -> float:
    """AUC that returns NaN instead of raising on degenerate inputs."""
    df = pd.DataFrame({"y": y_true.values, "s": y_score.values}).dropna()
    if len(df) < 2 or df["y"].nunique() < 2:
        return np.nan
    return roc_auc_score(df["y"], df["s"])


def _int(v) -> str:
    try:
        return str(int(v))
    except (ValueError, TypeError):
        return str(v).strip()


# ---- Load inputs ------------------------------------------------------------
run_info = pd.read_csv(RUN_INFO)
run_info.columns = [c.strip().lower() for c in run_info.columns]
labels = load_labels()

# records = []
# for _, r in run_info.iterrows():
#     site = str(r["site"]).strip()
#     disease = str(r["disease"]).strip()
#     k = _int(r["n_studies"])
#     rep = _int(r["repeat"])
#     test = str(r["test_study"]).strip()

#     path = os.path.join(
#         RESULT_DIR,
#         f"predictions_deepphylo_lc_{site}_{disease}_k{k}_r{rep}_test_{test}.csv",
#     )
#     if not os.path.exists(path):
#         print(f"[DeepPhylo] missing, skipped: {path}")
#         continue
#     pred = pd.read_csv(path)                        # SampleID, Prob_1, Pred
#     y_score = pred.set_index(pred["SampleID"].astype(str))["Prob_1"]
#     y_true = labels.reindex(y_score.index)
#     records.append({
#         "site": site, "disease": disease,
#         "n_studies": int(k),                 # pool size (incl. held-out test study)
#         "n_train_studies": int(k) - 1,       # studies actually used for training
#         "repeat": int(rep), "test_study": test,
#         "auc": safe_auc(y_true, y_score),
#     })

rf_raw = pd.read_csv("result/RF_results_lc.csv")
rf_raw["group"] = rf_raw["group"].astype(str)
records = []
for _, r in run_info.iterrows():
    site = str(r["site"]).strip()
    disease = str(r["disease"]).strip()
    k = _int(r["n_studies"])
    rep = _int(r["repeat"])
    test = str(r["test_study"]).strip()

    group = f"lc_{site}_{disease}_k{k}_r{rep}_test_{test}"
    sub = rf_raw[rf_raw["group"] == group]
    if sub.empty:
        print(f"[RF] no rows for group, skipped: {group}")
        continue
    sids = sub["sampel_id"].astype(str)
    y_score = sub.set_index(sids)["pred_prob"]
    y_true = labels.reindex(y_score.index)

    records.append({
        "site": site, "disease": disease,
        "n_studies": int(k),                 # pool size (incl. held-out test study)
        "n_train_studies": int(k) - 1,       # studies actually used for training
        "repeat": int(rep), "test_study": test,
        "auc": safe_auc(y_true, y_score),
    })
df = pd.DataFrame(records).dropna(subset=["auc"]).reset_index(drop=True)
assert not df.empty, (
    "No AUC values were computed. Check file naming and META_ID_COL / META_LABEL_COL."
)
print(f"[built] {len(df)} rows; sizes = {sorted(df['n_studies'].unique())}")
df.to_csv(os.path.join(RESULT_DIR, "learning_curve_auc_long_genus.csv"), index=False)


# =============================================================================
# PART 2 -- learning-curve plot
# =============================================================================

sns.set_theme(style="whitegrid", context="notebook")
mpl.rcParams.update({
    "font.family":    "DejaVu Sans",
    "axes.titlesize": 13, "axes.titleweight": "bold",
    "axes.labelsize": 11, "axes.edgecolor": "#666666", "axes.linewidth": 0.8,
    "grid.color": "#E6E6E6", "grid.linewidth": 0.8,
})

diseases = sorted(df["disease"].unique())
pal = dict(zip(diseases, sns.color_palette("tab10", n_colors=max(len(diseases), 3))))
sites = sorted(df["site"].unique())
xticks = sorted(df["n_train_studies"].unique())

n = len(sites)
ncols = min(n, 2)
nrows = math.ceil(n / ncols)
fig, axes = plt.subplots(nrows, ncols,
                         figsize=(ncols * 6.0, nrows * 4.6),
                         squeeze=False, sharey=True)
axes = axes.flatten()

rng = np.random.default_rng(0)
handles_labels = None

for i, site in enumerate(sites):
    ax = axes[i]
    d = df[df["site"] == site]

    # faint raw points (small horizontal jitter so repeats don't stack)
    for dis in diseases:
        s = d[d["disease"] == dis]
        if s.empty:
            continue
        jit = rng.uniform(-0.08, 0.08, len(s))
        ax.scatter(s["n_train_studies"] + jit, s["auc"], color=pal[dis],
                   s=14, alpha=0.25, edgecolors="none", zorder=2)

    # mean line + spread band per disease
    sns.lineplot(
        data=d, x="n_train_studies", y="auc", hue="disease", hue_order=diseases,
        palette=pal, estimator="mean", errorbar=ERRORBAR,
        marker="o", markersize=6, linewidth=2.0, err_style="band",
        ax=ax, legend=(i == 0), zorder=3,
    )

    ax.axhline(0.5, color="#BBBBBB", ls="--", lw=0.9, zorder=1)
    ax.set_title(site, pad=8)
    ax.set_xlabel("Number of training studies")
    ax.set_ylabel("AUC")
    ax.set_ylim(0, 1.03)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xticks(xticks)

    if i == 0 and ax.get_legend() is not None:
        handles_labels = ax.get_legend_handles_labels()
    if ax.get_legend() is not None:
        ax.get_legend().remove()

for j in range(n, len(axes)):
    axes[j].axis("off")

if handles_labels is not None:
    fig.legend(*handles_labels, title="Disease", loc="center left",
               bbox_to_anchor=(1.005, 0.5), frameon=True,
               fontsize=10, title_fontsize=11)

fig.suptitle("RF learning curve  -  AUC vs number of training studies",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 0.86, 0.96])
plt.savefig(os.path.join(RESULT_DIR, "learning_curve_genus.png"),
            dpi=200, bbox_inches="tight")
plt.show()

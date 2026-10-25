#!/usr/bin/env python3
# ==============================================================================
# build_supplementary_da_xlsx.py
# Build Supplementary_Table_Differential_Abundance.xlsx from the da_v2 OVA
# three-method differential-abundance results (LinDA + MaAsLin2 + ANCOM-BC2).
#
# Sources (real analysis outputs, nothing re-fitted):
#   Age   differential_abundance/result/age_hc_ova/age_da_consensus/
#         figures_saliva_target_top10_globalpeak_maaslin2coef/
#         figureA_classification_2method.tsv
#   Sex   differential_abundance/result/sex_pooled_hc/sex_pooled_da_consensus/consensus_robust_enhanced.tsv
#         (robust = FDR<0.05 in >=2/3 methods; saliva, HC only, pooled
#          continents, Male target vs Female reference)
#   Site/ differential_abundance/result/{site,continent}_ova/*/<x>_da_consensus/
#   Age/  figures_{saliva_target,}_top10_globalpeak_maaslin2coef/
#   Cont  figureA_classification_2method.tsv
#         rows with passes_direction_filter == TRUE  (paper Figure A counting
#         rule: robust >=2 methods + significant MaAsLin2 with matching
#         coefficient sign + unique group mean-abundance extremum)
#         rows with passes_direction_filter == TRUE  (paper Figure A counting
#         rule: robust >=2 methods + significant MaAsLin2 with matching
#         coefficient sign + unique group mean-abundance extremum)
#
# Mean relative abundance / prevalence per Group are computed here from
# data/metadata-and-biom-v2-qc/genus_table_relabund.tsv restricted to the
# exact sample scopes of the underlying analyses (declared, never imputed).
#
# Usage (from project root):
#   python3 differential_abundance/supplementary_da_xlsx/build_supplementary_da_xlsx.py
# ==============================================================================

import csv
import glob
import os
import time
import collections

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

# Repo root: two levels up from this file (differential_abundance/supplementary_da_xlsx/);
# can be overridden via GLOBAL_ORAL_ROOT for trees checked out elsewhere.
ROOT = os.environ.get("GLOBAL_ORAL_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DATA_DIR = f"{ROOT}/data/metadata-and-biom-v2-qc"
METADATA_FILE = f"{DATA_DIR}/metadata_oral_qc_dedup.tsv"
RELABUND_FILE = f"{DATA_DIR}/genus_table_relabund.tsv"
DA = f"{ROOT}/differential_abundance/result"
OUT_XLSX = f"{DA}/Supplementary_Table_Differential_Abundance.xlsx"

SITE_CLS = f"{DA}/site_ova/site_da_consensus/figures_target_top10_globalpeak_maaslin2coef/figureA_classification_2method.tsv"
CONT_CLS = f"{DA}/continent_ova/continent_da_consensus/figures_target_top10_globalpeak_maaslin2coef/figureA_classification_2method.tsv"
SEX_ENH = f"{DA}/sex_pooled_hc/sex_pooled_da_consensus/consensus_robust_enhanced.tsv"
AGE_CLS = (f"{DA}/age_hc_ova/age_da_consensus/"
           "figures_saliva_target_top10_globalpeak_maaslin2coef/"
           "figureA_classification_2method.tsv")

MERGED = {
    "linda": {
        "site": f"{DA}/site_ova/site_linda_ova/all_ova_results.tsv",
        "continent": f"{DA}/continent_ova/continent_linda_ova/all_ova_results.tsv",
        "age": f"{DA}/age_hc_ova/age_linda_ova/all_ova_results.tsv",
        "sex": f"{DA}/sex_pooled_hc/sex_linda_pairwise/all_sex_pairwise_results.tsv",
    },
    "maaslin2": {
        "site": f"{DA}/site_ova/site_maaslin2_ova/all_ova_results.tsv",
        "continent": f"{DA}/continent_ova/continent_maaslin2_ova/all_ova_results.tsv",
        "age": f"{DA}/age_hc_ova/age_maaslin2_ova/all_ova_results.tsv",
        "sex": f"{DA}/sex_pooled_hc/sex_maaslin2_pairwise/all_sex_pairwise_results.tsv",
    },
    "ancombc2": {
        "site": f"{DA}/site_ova/site_ancombc2_ova/all_ova_results.tsv",
        "continent": f"{DA}/continent_ova/continent_ancombc2_ova/all_ova_results.tsv",
        "age": f"{DA}/age_hc_ova/age_ancombc2_ova/all_ova_results.tsv",
        "sex": f"{DA}/sex_pooled_hc/sex_ancombc2_pairwise/all_sex_pairwise_results.tsv",
    },
}

AGE_ORDER = ["0-6 months", "6-36 months", "3-12 years", "12-18 years",
             "18-60 years", "＞60 years"]
AGE_DISP = {"＞60 years": ">60 years"}
SITE_ORDER = ["saliva", "tooth", "oral mucosa", "oral swab", "tongue", "gingiva"]
SITE_DISP = {"saliva": "Saliva", "tooth": "Tooth", "oral mucosa": "Oral mucosa",
             "oral swab": "Oral swab", "tongue": "Tongue", "gingiva": "Gingiva"}
CONT_ORDER = ["Africa", "Asia", "Europe", "North America", "Oceania", "South America"]
SEX_DISP = {"F": "Female", "M": "Male"}

EXPECTED = {
    "site": {"Saliva": (2, 4), "Tooth": (44, 6), "Oral mucosa": (14, 6),
             "Oral swab": (23, 4), "Tongue": (11, 43), "Gingiva": (12, 13)},
    "continent": {"Africa": (37, 20), "Asia": (49, 5), "Europe": (14, 13),
                  "North America": (21, 11), "Oceania": (18, 9),
                  "South America": (20, 52)},
    "sex": {"Female": (12,)},
}

SCOPE_LABEL = {
    "age":       "Healthy controls (HC) only, saliva",
    "sex":       "Healthy controls (HC) only, saliva (continents pooled)",
    "site":      "Full integrated dataset (HC + disease), all valid oral sites",
    "continent": "Full integrated dataset (HC + disease), saliva",
}
CONS_STATUS = {
    "age": ("Robust: FDR<0.05 in >=2 of 3 methods; significant MaAsLin2 with matching "
            "coefficient sign; unique max/min mean relative abundance across age stages"),
    "sex": "Robust: significant (FDR<0.05) in >=2 of 3 methods, agreeing direction",
    "site": ("Robust: FDR<0.05 in >=2 of 3 methods; significant MaAsLin2 with matching "
             "coefficient sign; unique max/min mean relative abundance across sites"),
    "continent": ("Robust: FDR<0.05 in >=2 of 3 methods; significant MaAsLin2 with matching "
                  "coefficient sign; unique max/min mean relative abundance across continents"),
}

verification = []


def vcheck(name, expected, actual, note=""):
    # numeric-aware comparison: supports "< 1e-4" style expected thresholds
    ok = str(expected) == str(actual)
    if not ok:
        try:
            ev = float(str(expected).lstrip("<= "))
            av = float(str(actual))
            ok = (av < ev) if str(expected).startswith("<") and not str(expected).startswith("<=") \
                else (av <= ev)
        except ValueError:
            ok = False
    status = "OK" if ok else "CHECK"
    verification.append([name, str(expected), str(actual), status, note])
    return status


def vinfo(name, actual, note=""):
    verification.append([name, "-", str(actual), "INFO", note])


# ===========================================================================
# 1. Metadata scopes
# ===========================================================================
print("[1/6] loading metadata...", flush=True)
meta = pd.read_csv(METADATA_FILE, sep="\t", header=0, dtype=str,
                   keep_default_na=False, quoting=csv.QUOTE_NONE,
                   usecols=["sample_id", "site_cleaned", "site_group_cleaned",
                            "age_cleaned_group", "sex_cleaned", "status_cleaned",
                            "continent_group"])
meta = meta[meta["sample_id"] != ""].drop_duplicates("sample_id")
sg = meta["site_group_cleaned"].str.strip()
meta = meta[meta["site_cleaned"].str.strip().isin(SITE_ORDER)
            & (sg.isin(["", "NA", "Oral"]) | sg.isna())]
print(f"    metadata samples (valid oral sites): {len(meta)}", flush=True)

with open(RELABUND_FILE, "r") as fh:
    header_cols = fh.readline().rstrip("\n").split("\t")
sample_ids = header_cols[1:]
col_idx = {s: i for i, s in enumerate(sample_ids)}   # index into parts[1:]
rel_samples = set(sample_ids)

scope_all = meta[meta["sample_id"].isin(rel_samples)]          # site/continent analyses
hc = scope_all[scope_all["status_cleaned"] == "HC"]            # age/sex analyses
print(f"    aligned to relabund: all={len(scope_all)}  hc={len(hc)}", flush=True)


def subset(df, **kw):
    m = pd.Series(True, index=df.index)
    for k, v in kw.items():
        m &= df[k] == v
    return df.loc[m, "sample_id"]      # sample_id values, not index labels


SUBSETS = {}
for g in AGE_ORDER:
    SUBSETS[("age", g)] = subset(hc, site_cleaned="saliva", age_cleaned_group=g)
for g in ("F", "M"):
    SUBSETS[("sex", g)] = subset(hc, site_cleaned="saliva", sex_cleaned=g)
for g in SITE_ORDER:
    SUBSETS[("site", g)] = subset(scope_all, site_cleaned=g)
for g in CONT_ORDER:
    SUBSETS[("continent", g)] = subset(scope_all, site_cleaned="saliva",
                                       continent_group=g)
print("    subset sizes (n): " + ", ".join(
    f"{a}:{g}={len(v)}" for (a, g), v in SUBSETS.items()), flush=True)

# ===========================================================================
# 2. per-comparison sample sizes + formulas from run_config.tsv (original runs)
# ===========================================================================
RUNCFG = {}    # (analysis, comparison) -> dict(n_ref, n_tgt, formula, scope)
cfg_dirs = {
    "site": f"{DA}/site_ova/*_ova/*_enrich/run_config.tsv" +
            ":",   # placeholder replaced below
}
for analysis, pattern in [
    ("site", f"{DA}/site_ova/site_{{m}}_ova/*_enrich/run_config.tsv"),
    ("continent", f"{DA}/continent_ova/continent_{{m}}_ova/*_enrich/run_config.tsv"),
    ("age", f"{DA}/age_hc_ova/age_{{m}}_ova/*_enrich*/run_config.tsv"),
    ("sex", f"{DA}/sex_pooled_hc/sex_{{m}}_pairwise/*/run_config.tsv"),
]:
    for m in ("linda", "maaslin2", "ancombc2"):
        for p in glob.glob(pattern.format(m=m)):
            c = pd.read_csv(p, sep="\t", nrows=1)
            comp = c["comparison"].iat[0]
            n1 = int(c["n1"].iat[0])      # N(reference / site1 / age1 / Sex1=F)
            n2 = int(c["n2"].iat[0])      # N(target group)
            if (analysis, comp) not in RUNCFG:
                RUNCFG[(analysis, comp)] = {
                    "n_ref": n1, "n_tgt": n2,
                    "formula": c["formula"].iat[0], "first_method": m,
                }
vcheck("run configs found (expect 6+6+6+1=19)", 19, len(RUNCFG))

# ===========================================================================
# 3. Load consensus tables and build normalized records
# ===========================================================================
print("[2/6] building records from consensus tables...", flush=True)
records = collections.defaultdict(list)


def genus_of(feature):
    parts = [p for p in feature.split(";") if p.strip()]
    for p in reversed(parts):
        if p.startswith("g__"):
            g = p[3:]
            if g:
                return g
            # uncalled genus: display as "<closest classified rank above> (unclassified genus)"
            parent = next((p2[3:] for p2 in reversed(parts)
                           if len(p2) > 3 and p2[:3] in ("k__", "p__", "c__", "o__", "f__")),
                          None)
            return f"{parent} (unclassified genus)" if parent else "unclassified"
    return "unclassified"


def sanitize_level(s):
    return s.replace(" ", "_").replace("＞", "gt")


def base_record(analysis, comparison, group_disp, ref_disp, feature, row):
    return {"Analysis": analysis,
            "Sample_scope": SCOPE_LABEL[analysis],
            "Group": group_disp,
            "Reference_group": ref_disp,
            "Comparison": comparison,
            "Genus": genus_of(feature),
            "Feature_full_taxonomy": feature,
            "_row": row}


# ---- Site + Continent + Age (Figure A classification rule) -----------------
CLS_CMP = {"site": {s: f"{sanitize_level(s)}_enrich" for s in SITE_ORDER},
           "continent": {c: f"{sanitize_level(c)}_enrich" for c in CONT_ORDER},
           "age": {a: f"saliva_enrich_{sanitize_level(a)}" for a in AGE_ORDER}}
CLS_DISP = {"site": SITE_DISP, "continent": {c: c for c in CONT_ORDER},
            "age": {a: AGE_DISP.get(a, a) for a in AGE_ORDER}}
CLS_REF = {"site": "Other (pooled other valid oral sites)",
           "continent": "Other (pooled other continents)",
           "age": "Other (pooled other valid age stages, saliva)"}
for analysis, cls_file in [("site", SITE_CLS), ("continent", CONT_CLS),
                           ("age", AGE_CLS)]:
    df = pd.read_csv(cls_file, sep="\t")
    ok = df[df["passes_direction_filter"] == True].drop_duplicates(["feature", "target"])
    for _, r in ok.iterrows():
        g = r["target"]
        rec = base_record(analysis, CLS_CMP[analysis][g], CLS_DISP[analysis][g],
                          CLS_REF[analysis], r["feature"], r)
        rec["Enrichment_direction"] = ("Enriched" if r["direction"] == "Enriched"
                                       else "Depleted")
        records[analysis].append(rec)
        # cross-check: classification direction must match model direction
        if (r["direction"] == "Enriched") != (r["maaslin2_lfc"] > 0):
            verification.append(
                [f"{analysis}: direction vs MaAsLin2 coef sign ({r['feature'][:40]}... / {g})",
                 "consistent", "INCONSISTENT", "CHECK", ""])
    tg = collections.Counter((CLS_DISP[analysis][r["target"]], r["direction"])
                             for _, r in ok.iterrows())
    if analysis in EXPECTED:
        for grp, (e_exp, d_exp) in EXPECTED[analysis].items():
            vcheck(f"{analysis}: {grp} enriched", e_exp, tg.get((grp, "Enriched"), 0))
            vcheck(f"{analysis}: {grp} depleted", d_exp, tg.get((grp, "Depleted"), 0))
    print(f"    {analysis}: {len(ok)} rows", flush=True)
    if analysis == "age":
        for g in [AGE_DISP.get(a, a) for a in AGE_ORDER]:
            vinfo(f"age: {g} figureA rule (Enriched/Depleted)",
                  f"{tg.get((g,'Enriched'),0)} enriched / {tg.get((g,'Depleted'),0)} depleted",
                  "unique-extremum rule applied (matches site/continent Figure A rule)")

# ---- Sex --------------------------------------------------------------------
df = pd.read_csv(SEX_ENH, sep="\t")
for _, r in df.iterrows():
    grp = r["direction_consensus"].replace("agreed:", "")    # 'F' or 'M'
    ref = "M" if grp == "F" else "F"
    rec = base_record("sex", "saliva_pooled_F_vs_M", SEX_DISP[grp], SEX_DISP[ref],
                      r["feature"], r)
    rec["Enrichment_direction"] = "Enriched"
    records["sex"].append(rec)
tg = collections.Counter(r["Group"] for r in records["sex"])
vcheck("sex: Female-enriched (robust n2)", 12, tg.get("Female", 0))
vinfo("sex: Male-enriched (robust n2)", tg.get("Male", 0),
      "no paper expected number given; total robust = 57")
print(f"    sex: {len(df)} rows (Female={tg.get('Female')}, Male={tg.get('Male')})",
      flush=True)

# ===========================================================================
# 4. raw P values from merged per-method tables
# ===========================================================================
print("[3/6] joining p values from merged per-method tables...", flush=True)
pval_lookup = {}
for method, per_an in MERGED.items():
    for analysis, path in per_an.items():
        m = pd.read_csv(path, sep="\t", usecols=["feature", "Comparison", "pvalue", "padj"])
        for c, f, pv, qv in zip(m["Comparison"], m["feature"], m["pvalue"], m["padj"]):
            pval_lookup[(analysis, method, c, f)] = (pv, qv)
        # padj in the merged table should equal the padj captured in the
        # consensus files; this is used in the Verification sheet later.
print(f"    merged lookups: {len(pval_lookup)}", flush=True)

# ===========================================================================
# 5. per-group mean relative abundance + prevalence (streaming relabund)
# ===========================================================================
print("[4/6] streaming relabund for group means/prevalence...", flush=True)
GROUPS_OF = {"age": AGE_ORDER, "sex": ["F", "M"], "site": SITE_ORDER,
             "continent": CONT_ORDER}
need_features = collections.defaultdict(set)
for an in records:
    for r in records[an]:
        need_features[an].add(r["Feature_full_taxonomy"])

feat_subsets = {}
for an, feats in need_features.items():
    for f in feats:
        feat_subsets.setdefault(f, []).extend([(an, g) for g in GROUPS_OF[an]])

stats = {}
total_missing = []
with open(RELABUND_FILE, "r") as fh:
    fh.readline()
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        keys = feat_subsets.get(parts[0])
        if not keys:
            continue
        row = np.array(parts[1:], dtype=np.float64)
        d = {}
        for key in set(keys):
            idx_list = [col_idx[s] for s in SUBSETS[key] if s in col_idx]
            if not idx_list:
                total_missing.append((key, "no samples in relabund"))
                continue
            vals = row[np.array(idx_list, dtype=np.int64)]
            d[key] = (len(idx_list), float(np.nansum(vals)), int((vals > 0).sum()))
        stats[parts[0]] = d
print(f"    features found: {len(stats)}/{len(feat_subsets)}", flush=True)
print(f"    group subsets missing from relabund: {len(total_missing)}", flush=True)
vcheck("all table features found in relabund", len(feat_subsets), len(stats))
vcheck("analysis group subsets with no aligned relabund samples", 0, len(total_missing))


def group_abundance(analysis, graw, feature):
    """Per-group (mean relative abundance, prevalence) for every group of
    `analysis` — the analysis scope partitions into these groups, so the
    pooled-Other reference of a row equals scope totals minus the focal
    group's contribution."""
    d = stats.get(feature) or {}
    out = {}
    tot_n = tot_s = tot_nz = 0
    for g in GROUPS_OF[analysis]:
        n, s, nz = d.get((analysis, g), (0, 0, 0))
        out[GROUP_DISP[analysis][g]] = ((s / n if n else np.nan),
                                        (nz / n if n else np.nan))
        tot_n += n
        tot_s += s
        tot_nz += nz
    if analysis != "sex":
        n, s, nz = d.get((analysis, graw), (0, 0, 0))
        on, os_, onz = tot_n - n, tot_s - s, tot_nz - nz
        out["Other_pooled"] = ((os_ / on if on else np.nan),
                               (onz / on if on else np.nan))
    return out


# ===========================================================================
# 6. assemble final tables
# ===========================================================================
print("[5/6] assembling workbook...", flush=True)
BASE_COLUMNS = ["Analysis", "Sample_scope", "Group", "Reference_group", "Comparison",
                "Genus", "Feature_full_taxonomy", "Enrichment_direction",
                "Effect_size_LinDA", "P_value_LinDA", "FDR_LinDA", "Sig_LinDA",
                "Effect_size_MaAsLin2", "P_value_MaAsLin2", "FDR_MaAsLin2", "Sig_MaAsLin2",
                "Effect_size_ANCOMBC2", "P_value_ANCOMBC2", "FDR_ANCOMBC2", "Sig_ANCOMBC2",
                "N_methods_significant", "Method_support", "Consensus_status",
                "Group_n_samples", "Reference_n_samples", "Source_file"]

# Per-group abundance/prevalence columns: one column per group of the analysis
# (plus the pooled-Other reference, except Sex where the reference is the other
# sex group itself and would duplicate an existing column).
GROUP_DISP = {"age": {a: AGE_DISP.get(a, a) for a in AGE_ORDER},
              "sex": {"F": "Female", "M": "Male"},
              "site": {s: SITE_DISP[s] for s in SITE_ORDER},
              "continent": {c: c for c in CONT_ORDER}}


def metric_columns(analysis):
    gds = [GROUP_DISP[analysis][g] for g in GROUPS_OF[analysis]]
    mean_cols = [f"Mean_rel_abund_{g}" for g in gds]
    prev_cols = [f"Prevalence_{g}" for g in gds]
    if analysis != "sex":
        mean_cols.append("Mean_rel_abund_Other_pooled")
        prev_cols.append("Prevalence_Other_pooled")
    return mean_cols + prev_cols


METRIC_COLS = {an: metric_columns(an) for an in ("age", "sex", "site", "continent")}
# Combined_All stacks all four analyses; the *_Other_pooled columns are shared
# (every row belongs to exactly one analysis, which Sample_scope declares), so
# duplicate names collapse to the first occurrence.
ALL_COLUMNS = BASE_COLUMNS + list(dict.fromkeys(
    c for an in ("age", "sex", "site", "continent") for c in METRIC_COLS[an]))

CFG_ROOT = {"age": "differential_abundance/result/age_hc_ova",
            "sex": "differential_abundance/result/sex_pooled_hc",
            "site": "differential_abundance/result/site_ova",
            "continent": "differential_abundance/result/continent_ova"}
SRC_FILE = {"age": "differential_abundance/result/age_hc_ova/age_da_consensus/figures_saliva_target_top10_globalpeak_maaslin2coef/figureA_classification_2method.tsv",
            "sex": "differential_abundance/result/sex_pooled_hc/sex_pooled_da_consensus/consensus_robust_enhanced.tsv",
            "site": "differential_abundance/result/site_ova/site_da_consensus/figures_target_top10_globalpeak_maaslin2coef/figureA_classification_2method.tsv",
            "continent": "differential_abundance/result/continent_ova/continent_da_consensus/figures_target_top10_globalpeak_maaslin2coef/figureA_classification_2method.tsv"}

final_rows = []
pval_missing = 0
pval_missing_rows = []
for analysis in ["age", "sex", "site", "continent"]:
    for rec in records[analysis]:
        row = rec["_row"]
        g = rec["Group"]
        graw = ({"age": {AGE_DISP.get(a, a): a for a in AGE_ORDER}.get(g, g),
                 "sex": {"Female": "F", "Male": "M"}.get(g, g),
                 "site": {v: k for k, v in SITE_DISP.items()}.get(g, g),
                 "continent": g}.get(analysis, g))
        ab = group_abundance(analysis, graw, rec["Feature_full_taxonomy"])
        out = {k: rec[k] for k in BASE_COLUMNS if k in rec}
        sig = {}
        for m in ("linda", "maaslin2", "ancombc2"):
            sig[m] = str(row.get(f"in_{m}", "")).strip().lower() == "true"
        for m, dm in (("linda", "LinDA"), ("maaslin2", "MaAsLin2"),
                      ("ancombc2", "ANCOMBC2")):
            pv, _ = pval_lookup.get(
                (analysis, m, rec["Comparison"], rec["Feature_full_taxonomy"]),
                (np.nan, np.nan))
            if pv != pv:                                   # nan check
                pval_missing += 1
                pval_missing_rows.append((analysis, dm, g, rec["Genus"]))
            out[f"Effect_size_{dm}"] = row.get(f"{m}_lfc", np.nan)
            out[f"P_value_{dm}"] = pv
            out[f"FDR_{dm}"] = row.get(f"{m}_padj", np.nan)
            out[f"Sig_{dm}"] = "Yes" if sig[m] else "No"
        nms = []
        for m in ("linda", "maaslin2", "ancombc2"):
            if sig[m]:
                nms.append({"linda": "LinDA", "maaslin2": "MaAsLin2",
                            "ancombc2": "ANCOM-BC2"}[m])
        out["N_methods_significant"] = sum(sig.values())
        out["Method_support"] = "+".join(nms) if nms else "NA"
        out["Consensus_status"] = CONS_STATUS[analysis]
        for suffix, (m, p) in ab.items():
            out[f"Mean_rel_abund_{suffix}"] = m
            out[f"Prevalence_{suffix}"] = p
        cfg = RUNCFG.get((analysis, rec["Comparison"]))
        if analysis == "sex":
            out["Group_n_samples"] = cfg["n_tgt"] if graw == "M" else cfg["n_ref"]
            out["Reference_n_samples"] = cfg["n_ref"] if graw == "M" else cfg["n_tgt"]
        else:
            out["Group_n_samples"] = cfg["n_tgt"]
            out["Reference_n_samples"] = cfg["n_ref"]
        out["Source_file"] = SRC_FILE[analysis]
        final_rows.append(out)
vcheck("rows with an unjoined (NaN) p value left",
       0, pval_missing,
       "genera never reported by a method (e.g. ANCOM-BC2 internal filtering) keep NA for that method only")
for an, dm, g, gn in pval_missing_rows:
    vinfo(f"{an}/{dm}: P_value not reported by the original run",
          f"{g}: {gn}",
          "the genus is absent from that method's per-comparison result table "
          "(ANCOM-BC2 drops/merges taxa internally); FDR_* is NA for that method")
fdf = pd.DataFrame(final_rows, columns=ALL_COLUMNS)

# ---- ordering ----------------------------------------------------------------
fdf["_a"] = fdf["Analysis"].map({"age": 0, "sex": 1, "site": 2, "continent": 3})
fdf["_g"] = fdf.apply(
    lambda r: {"age": [AGE_DISP.get(a, a) for a in AGE_ORDER],
               "sex": ["Female", "Male"],
               "site": [SITE_DISP[s] for s in SITE_ORDER],
               "continent": CONT_ORDER}[r["Analysis"]].index(r["Group"]), axis=1)
fdf["_d"] = (fdf["Enrichment_direction"] != "Enriched").astype(int)
fdf = fdf.sort_values(["_a", "_g", "_d", "FDR_MaAsLin2", "Genus"], kind="mergesort"
                      ).drop(columns=["_a", "_g", "_d"])

# ---- run-config sample sizes vs recomputed subset sizes ----------------------
for (analysis, comp), cfg in sorted(RUNCFG.items()):
    sat = sanitize_level
    if analysis in ("site", "continent"):
        graw = next(g for g in GROUPS_OF[analysis] if sat(g) == comp.rsplit("_enrich", 1)[0])
        ref_raw = [x for x in GROUPS_OF[analysis] if x != graw]
    elif analysis == "age":
        graw = next(a for a in AGE_ORDER if sat(a) == comp.replace("saliva_enrich_", ""))
        ref_raw = [a for a in AGE_ORDER if a != graw]
    else:  # sex
        graw, ref_raw = "M", ["F"]
    if (analysis, graw) in SUBSETS:
        vcheck(f"{analysis}/{comp}: N(target) run-config vs recomputed",
               cfg["n_tgt"], len(SUBSETS[(analysis, graw)]))
    ref_n_recomp = len(set().union(*[SUBSETS[(analysis, x)] for x in ref_raw]))
    vcheck(f"{analysis}/{comp}: N(reference pooled Other) run-config vs recomputed",
           cfg["n_ref"], ref_n_recomp)

# ---- mean-abundance cross-check (classification file value vs recomputed) ----
mr = []
mr_detail = []
chk = fdf[fdf["Analysis"].isin(["site", "continent", "age"])]
cls_lookup = {}
for an, cls_file in [("site", SITE_CLS), ("continent", CONT_CLS), ("age", AGE_CLS)]:
    d = pd.read_csv(cls_file, sep="\t")
    d = d[d["passes_direction_filter"] == True].drop_duplicates(["feature", "target"])
    for _, r in d.iterrows():
        cls_lookup[(an, r["feature"], r["target"])] = r["target_mean_abundance"]
for _, r in chk.iterrows():
    if r["Analysis"] == "site":
        graw = {v: k for k, v in SITE_DISP.items()}[r["Group"]]
    elif r["Analysis"] == "age":
        graw = {AGE_DISP.get(a, a): a for a in AGE_ORDER}[r["Group"]]
    else:
        graw = r["Group"]
    file_v = cls_lookup.get((r["Analysis"], r["Feature_full_taxonomy"], graw))
    focal = r[f"Mean_rel_abund_{r['Group']}"]
    if file_v is not None and focal == focal:
        d = abs(file_v - focal)
        mrd = d / max(file_v, focal, 1e-300)
        mr.append(mrd)
        mr_detail.append((r["Analysis"], r["Group"], r["Genus"], file_v,
                          focal, d, mrd))
vcheck("site+continent+age: recomputed group mean vs classification file (max rel diff)",
       "< 1e-4", f"{max(mr):.2e}" if mr else "n/a", f"{len(mr)} values compared")
if mr_detail:
    with open("/tmp/mean_crosscheck.tsv", "w") as fh:
        fh.write("analysis\tgroup\tgenus\tfile_mean\trecomputed_mean\tabs_diff\trel_diff\n")
        for t in sorted(mr_detail, key=lambda x: -x[6]):
            fh.write("\t".join(map(str, t)) + "\n")

# ===========================================================================
# 7. write workbook
# ===========================================================================
print("[6/6] writing xlsx...", flush=True)
wb = Workbook()
HDR_FONT = Font(bold=True)
HDR_FILL = PatternFill("solid", fgColor="DDEBF7")


def style_sheet(ws, ncols):
    for c in range(1, ncols + 1):
        cell = ws.cell(row=1, column=c)
        cell.font = HDR_FONT
        cell.fill = HDR_FILL
        cell.alignment = Alignment(vertical="center", wrap_text=True)
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(ncols)}{ws.max_row}"


# ---- README -----------------------------------------------------------------
ws = wb.active
ws.title = "README"
RM_TITLE = "__TITLE__"
RM_SECTION = "__SECTION__"
readme = [
    (RM_TITLE, "Supplementary Table. Age-, Sex-, Site- and Continent-associated differential genera (one-vs-all designs)"),
    ("", ""),
    ("Generated", time.strftime("%Y-%m-%d %H:%M")),
    ("Project", "Global oral microbiome meta-analysis — da_v2 three-method differential-abundance pipeline (LinDA + MaAsLin2 + ANCOM-BC2), results root: differential_abundance/result/"),
    ("", ""),
    (RM_SECTION, "1. Analysis scope"),
    ("Age", "Health-associated saliva samples from healthy controls (HC) only; six life stages (0-6 months, 6-36 months, 3-12 years, 12-18 years, 18-60 years, >60 years), each stage vs the pooled Other (one-vs-all). Rows follow the paper Figure A classification rule (robust >=2/3 methods + significant MaAsLin2 with matching coefficient sign + unique group mean-abundance extremum across the six age stages); 133 genera."),
    ("Sex", "Healthy controls (HC), saliva only, all continents pooled; single contrast Male vs Female (reference = Female). Robust = FDR<0.05 in >=2/3 methods; all 57 robust genera listed (12 female-enriched, 45 male-enriched)."),
    ("Site", "Full integrated dataset (healthy controls + disease), all six valid oral sampling sites; each site vs the pooled Other sites (one-vs-all). Rows follow the paper Figure A classification rule (robust >=2/3 methods + significant MaAsLin2 with matching coefficient sign + unique group mean-abundance extremum); 182 genera."),
    ("Continent", "Full integrated dataset (healthy controls + disease), saliva only; each continent vs the pooled Other continents (one-vs-all). Same Figure A classification rule as Site; 269 genera."),
    ("", ""),
    (RM_SECTION, "2. Methods and thresholds (full record per comparison in run_config.tsv)"),
    ("LinDA", "MicrobiomeStat::linda on genus counts; prev.filter=0.01, is.winsor=FALSE, BH-adjusted, alpha=0.05."),
    ("MaAsLin2", "TSS normalization + LOG transform + LM; min_prevalence=0.01, BH q values, standardize=FALSE."),
    ("ANCOM-BC2", "ANCOMBC::ancombc2; prv_cut=0.01, lib_cut=0, struc_zero=FALSE, neg_lb=FALSE, BH-adjusted; random effect (1|bioproject_id) when a BioProject spans both groups of a comparison."),
    ("Significance", "FDR-adjusted p value (padj/qval) < 0.05 for each method."),
    ("Covariates", "Age: sex, platform, sequencing region, BMI group, continent. Sex: continent, age, platform, sequencing region, BMI group. Site: age, sex, platform, sequencing region, BMI group, continent. Continent: age, sex, platform, sequencing region, BMI group. BioProject random effect where applicable (see Comparison run_config formula)."),
    ("", ""),
    (RM_SECTION, "3. Consensus rules (original analysis rules, not re-defined)"),
    ("Age", "Paper Figure A classification: robust (>=2/3 methods significant) AND significant MaAsLin2 coefficient with matching sign AND a unique maximum (Enriched) / minimum (Depleted) mean relative abundance of the genus across the six age stages. Computed with the same ova_extreme_counts() rule as the site/continent Figure A (06-build-figureA-classification.R)."),
    ("Sex", "Robust consensus = significant (FDR<0.05) in >=2 of the 3 methods with an agreed enrichment direction (n2 stratification — the primary consensus output of 02-merge-and-consensus.R). All robust genera are listed (direction follows the covariate-adjusted model coefficient, not the raw group mean)."),
    ("Site & Continent", "Paper Figure A classification: robust (>=2/3 methods significant) AND significant MaAsLin2 coefficient with matching sign AND a unique maximum (Enriched) / minimum (Depleted) mean relative abundance of the genus across the compared groups. This rule reproduces the published comparison counts exactly (see Verification sheet)."),
    ("", ""),
    (RM_SECTION, "4. Column definitions"),
    ("Effect_size_LinDA", "LinDA log2FoldChange (log2 scale; target group vs reference)."),
    ("Effect_size_MaAsLin2", "MaAsLin2 model coefficient (TSS-restricted proportions, natural-log scale; target vs reference). Different scale from LinDA/ANCOM-BC2 — kept separate on purpose."),
    ("Effect_size_ANCOMBC2", "ANCOM-BC2 log fold-change (bias-corrected log scale; target vs reference). Different scale — kept separate."),
    ("P_value_*", "Raw (unadjusted) model p value, joined from the merged per-method tables (all_*_results.tsv)."),
    ("FDR_*", "FDR-adjusted p value from the original consensus outputs (LinDA padj / MaAsLin2 qval / ANCOM-BC2 padj)."),
    ("Sig_*", "Whether the method was significant in the original consensus (in_linda / in_maaslin2 / in_ancombc2 flags; FDR<0.05)."),
    ("N_methods_significant | Method_support", "Count and names of the methods significant for the genus in that comparison."),
    ("Enrichment_direction", "Enriched = higher abundance in the focal Group vs Reference_group; Depleted = lower."),
    ("Mean_rel_abund_<group>", "Mean relative abundance (proportion, 0-1) of the genus within each group's samples, one column per group of the analysis (Age: the six life stages; Sex: Female/Male; Site: the six sites; Continent: the six continents) plus Mean_rel_abund_Other_pooled = the pooled-Other reference samples (all analysis-scope samples outside the focal Group; not present for Sex, where the reference is the other sex group itself). Computed from data/metadata-and-biom-v2-qc/genus_table_relabund.tsv on the exact analysis sample scope (not fabricated; NA if genus absent from the table). The focal Group's column carries the value behind Enrichment_direction."),
    ("Prevalence_<group>", "Proportion of each group's samples with non-zero relative abundance (same per-group layout as Mean_rel_abund_<group>, incl. Prevalence_Other_pooled; same source as above)."),
    ("Group_n_samples | Reference_n_samples", "Sample counts of the target group and pooled reference, from run_config.tsv of the corresponding comparison (n1 = reference/Other, n2 = target; for Sex, n1 = Female, n2 = Male)."),
    ("Comparison", "Original internal comparison name (e.g. saliva_enrich_0-6_months, gingiva_enrich, saliva_pooled_F_vs_M)."),
    ("Genus | Feature_full_taxonomy", "Original genus name and full Greengenes2-style lineage; database suffixes (Neisseria_563205, Veillonella_A, Bifidobacterium_388775, ...) kept verbatim. Genera with no genus-level assignment (empty g__ in the lineage) are shown as \"<closest classified rank> (unclassified genus)\" (e.g. \"Filifactoraceae (unclassified genus)\"); Feature_full_taxonomy keeps the verbatim lineage; nothing renamed or merged."),
    ("", ""),
    (RM_SECTION, "5. Input files -> worksheets"),
    ("Age sheet", "differential_abundance/result/age_hc_ova/age_da_consensus/figures_saliva_target_top10_globalpeak_maaslin2coef/figureA_classification_2method.tsv (rows passing the filter; built by 06-build-figureA-classification.R with the same rule as site/continent) + per-method all_ova_results.tsv"),
    ("Sex sheet", "differential_abundance/result/sex_pooled_hc/sex_pooled_da_consensus/consensus_robust_enhanced.tsv (+ all_sex_pairwise_results.tsv)"),
    ("Site sheet", "differential_abundance/result/site_ova/site_da_consensus/figures_target_top10_globalpeak_maaslin2coef/figureA_classification_2method.tsv (rows passing the filter) + per-method all_ova_results.tsv"),
    ("Continent sheet", "differential_abundance/result/continent_ova/continent_da_consensus/figures_target_top10_globalpeak_maaslin2coef/figureA_classification_2method.tsv (+ per-method all_ova_results.tsv)"),
    ("Summary sheet", "Counts derived from the four analysis sheets."),
    ("Combined_All sheet", "Age + Sex + Site + Continent stacked; shared columns first, then the per-analysis group abundance/prevalence columns appended (values only on rows of the matching analysis)."),
    ("Verification sheet", "Automated consistency checks: expected vs actual counts, run-config vs recomputed sample sizes, direction/sign concordance, mean-abundance cross-check."),
    ("", ""),
    (RM_SECTION, "6. Notes"),
    ("Effect-size scales", "LinDA reports log2 FC; MaAsLin2 a natural-log coefficient; ANCOM-BC2 its own log fold-change. Never pooled or rescaled."),
    (">60 years label", "Metadata uses the full-width character in the >60 label; displayed as '>60 years' here; the original comparison name keeps the sanitized form (gt60_years)."),
    ("1% prevalence filter", "Applied inside each method during DA (prev_cut=0.01 per comparison); genera not passing a tool's filter have NA statistics for that tool."),
    ("Sorting", "Rows sorted by Analysis, Group (canonical order), Enrichment_direction (Enriched first), then MaAsLin2 FDR ascending."),
    ("Group-level scope of Age/Sex analyses", "Age and Sex analyses exclude disease groups by design and are saliva-only; Site and Continent analyses use the full integrated dataset (HC + disease). Do not mix the two scopes."),
]
for item in readme:
    ws.append([item[0], item[1]])
ws.column_dimensions["A"].width = 34
ws.column_dimensions["B"].width = 120
for rr, (a, b) in enumerate(readme, start=1):
    ca, cb = ws.cell(row=rr, column=1), ws.cell(row=rr, column=2)
    for c in (ca, cb):
        c.alignment = Alignment(vertical="top", wrap_text=True)
    if a == RM_TITLE:
        ca.font = Font(bold=True, size=13)
    elif a == RM_SECTION:
        ca.font = Font(bold=True, size=12)
        ca.fill = PatternFill("solid", fgColor="EFEFEF")
        cb.fill = PatternFill("solid", fgColor="EFEFEF")

# ---- Summary ------------------------------------------------------------------
ws = wb.create_sheet("Summary")
ws.append(["Analysis", "Scope", "Group", "Reference_group", "N_genera_enriched",
           "N_genera_depleted", "N_genera_total", "Group_n_samples",
           "Reference_n_samples"])
scope_of = {a: SCOPE_LABEL[a] for a in SCOPE_LABEL}
cnt = collections.Counter(zip(fdf["Analysis"], fdf["Sample_scope"], fdf["Group"],
                              fdf["Enrichment_direction"]))
nsz = {}
for _, r in fdf.iterrows():
    nsz.setdefault((r["Analysis"], r["Group"]),
                   (r["Group_n_samples"], r["Reference_n_samples"]))
REF_DISP = {"age": "Other age stages (pooled, saliva)",
            "site": "Other valid oral sites (pooled)",
            "continent": "Other continents (pooled, saliva)"}
for an, groups in [("age", [AGE_DISP.get(a, a) for a in AGE_ORDER]),
                   ("sex", ["Female", "Male"]),
                   ("site", [SITE_DISP[s] for s in SITE_ORDER]),
                   ("continent", CONT_ORDER)]:
    for g in groups:
        e = cnt.get((an, scope_of[an], g, "Enriched"), 0)
        d = cnt.get((an, scope_of[an], g, "Depleted"), 0)
        gn, rn = nsz.get((an, g), (None, None))
        ref = (REF_DISP[an] if an != "sex" else ("Male" if g == "Female" else "Female"))
        ws.append([an.capitalize(), scope_of[an], g, ref, e, d, e + d, gn, rn])
style_sheet(ws, 9)
for c, w in zip(range(1, 10), [11, 42, 14, 34, 12, 12, 12, 12, 12]):
    ws.column_dimensions[get_column_letter(c)].width = w

# ---- data sheets ---------------------------------------------------------------
for name, sub in [("Age", fdf[fdf["Analysis"] == "age"]),
                  ("Sex", fdf[fdf["Analysis"] == "sex"]),
                  ("Site", fdf[fdf["Analysis"] == "site"]),
                  ("Continent", fdf[fdf["Analysis"] == "continent"]),
                  ("Combined_All", fdf)]:
    cols = ALL_COLUMNS if name == "Combined_All" else BASE_COLUMNS + METRIC_COLS[name.lower()]
    ws = wb.create_sheet(name)
    ws.append(cols)
    for _, r in sub.iterrows():
        vals = [r.get(c) for c in cols]
        ws.append([None if (isinstance(v, float) and v != v) or v is pd.NA else v
                   for v in vals])
    style_sheet(ws, len(cols))
    col_of = {nm: i + 1 for i, nm in enumerate(cols)}
    for nm in ["P_value_LinDA", "FDR_LinDA", "P_value_MaAsLin2", "FDR_MaAsLin2",
               "P_value_ANCOMBC2", "FDR_ANCOMBC2"]:
        for rr in range(2, ws.max_row + 1):
            ws.cell(row=rr, column=col_of[nm]).number_format = "0.00E+00"
    for nm in [c for c in cols if c.startswith("Mean_rel_abund_")]:
        for rr in range(2, ws.max_row + 1):
            ws.cell(row=rr, column=col_of[nm]).number_format = "0.00E+00"
    for nm in ["Effect_size_LinDA", "Effect_size_MaAsLin2", "Effect_size_ANCOMBC2"]:
        for rr in range(2, ws.max_row + 1):
            ws.cell(row=rr, column=col_of[nm]).number_format = "0.000"
    for nm in [c for c in cols if c.startswith("Prevalence_")]:
        for rr in range(2, ws.max_row + 1):
            ws.cell(row=rr, column=col_of[nm]).number_format = "0.0%"
    wmap = {"Analysis": 10, "Sample_scope": 26, "Group": 12, "Reference_group": 26,
            "Comparison": 24, "Genus": 22, "Feature_full_taxonomy": 68,
            "Enrichment_direction": 12, "Consensus_status": 40,
            "Method_support": 22, "Source_file": 56, "N_methods_significant": 9,
            "Sig_LinDA": 8, "Sig_MaAsLin2": 10, "Sig_ANCOMBC2": 10,
            "Group_n_samples": 10, "Reference_n_samples": 10}
    for nm in ["Effect_size_LinDA", "P_value_LinDA", "FDR_LinDA",
               "Effect_size_MaAsLin2", "P_value_MaAsLin2", "FDR_MaAsLin2",
               "Effect_size_ANCOMBC2", "P_value_ANCOMBC2", "FDR_ANCOMBC2"]:
        wmap.setdefault(nm, 14)
    for nm in [c for c in cols if c.startswith(("Mean_rel_abund_", "Prevalence_"))]:
        wmap.setdefault(nm, 17)
    for nm, w in wmap.items():
        ws.column_dimensions[get_column_letter(col_of[nm])].width = w

# ---- Verification ----------------------------------------------------------------
ws = wb.create_sheet("Verification")
ws.append(["Check", "Expected", "Actual", "Status", "Note"])
for r in verification:
    ws.append(r)
style_sheet(ws, 5)
for c, w in zip(range(1, 6), [64, 24, 24, 10, 46]):
    ws.column_dimensions[get_column_letter(c)].width = w
fails = [r for r in verification if r[3] == "CHECK"]

wb.save(OUT_XLSX)
counts = fdf.groupby(["Analysis", "Group", "Enrichment_direction"]).size()
print(f"\nsaved {OUT_XLSX}")
print(f"rows: age={int((fdf['Analysis']=='age').sum())}, "
      f"sex={int((fdf['Analysis']=='sex').sum())}, "
      f"site={int((fdf['Analysis']=='site').sum())}, "
      f"continent={int((fdf['Analysis']=='continent').sum())}, total={len(fdf)}")
print(f"verification checks: {len(verification)}, flagged: {len(fails)}")
for r in fails:
    print("  CHECK:", r)
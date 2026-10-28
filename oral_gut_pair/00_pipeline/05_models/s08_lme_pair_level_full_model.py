#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s08_lme_pair_level_full_model.py

Fit full mixed-effects models for pair-level oral-gut outcomes in the
11 time-series BioProjects.

Recommended model:
  outcome ~ BioProject + sampling_time_z_by_bp + selected covariates
            + (1 | BioProject:subject_id)

Implementation:
  Python statsmodels MixedLM with a random intercept for subject_block.
  Fixed effects are assembled as a design matrix. Rank-deficient columns are
  automatically removed by QR decomposition before model fitting.
"""

import argparse
import json
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


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
    if s in {"1", "true", "t", "yes", "y"}:
        return True
    if s in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Cannot parse boolean value: {x}")


def clean_missing(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "na", "none", "null"}:
        return np.nan
    return s


def bh_fdr(pvals):
    p = pd.to_numeric(pd.Series(pvals), errors="coerce").values.astype(float)
    q = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    q[order] = adj
    return q


def choose_existing_columns(df, cols):
    return [c for c in cols if c in df.columns]


def clean_categorical_series(s, missing_level="Missing"):
    out = s.map(clean_missing)
    out = out.astype("object")
    out = out.where(out.notna(), missing_level)
    return out.astype(str)


def has_within_bp_variation(df, col, bp_col):
    if col not in df.columns or bp_col not in df.columns:
        return False
    tmp = df[[bp_col, col]].copy()
    tmp[col] = clean_categorical_series(tmp[col])
    nun = tmp.groupby(bp_col)[col].nunique(dropna=False)
    return bool((nun >= 2).any())


def make_fixed_effect_matrix(
    df,
    outcome,
    time_var,
    bioproject_col,
    categorical_covariates,
    numeric_covariates,
    include_bioproject=True,
    missing_level="Missing",
    drop_no_within_bp_variation=True,
):
    """Return exog matrix, term metadata, and dropped-term records."""
    dropped = []
    term_meta = []

    exog_parts = []

    intercept = pd.DataFrame({"Intercept": np.ones(df.shape[0], dtype=float)}, index=df.index)
    exog_parts.append(intercept)
    term_meta.append({"column": "Intercept", "source_term": "Intercept", "term_type": "intercept"})

    # BioProject fixed effect.
    if include_bioproject:
        if bioproject_col not in df.columns:
            dropped.append({"term": bioproject_col, "term_type": "categorical", "reason": "missing_column"})
        else:
            s = clean_categorical_series(df[bioproject_col], missing_level=missing_level)
            if s.nunique(dropna=False) < 2:
                dropped.append({"term": bioproject_col, "term_type": "categorical", "reason": "single_level"})
            else:
                dummies = pd.get_dummies(s, prefix=f"C({bioproject_col})", drop_first=True, dtype=float)
                exog_parts.append(dummies)
                for c in dummies.columns:
                    term_meta.append({"column": c, "source_term": bioproject_col, "term_type": "bioproject_fixed_effect"})

    # Time fixed effect.
    if time_var:
        if time_var not in df.columns:
            dropped.append({"term": time_var, "term_type": "numeric", "reason": "missing_column"})
        else:
            x = pd.to_numeric(df[time_var], errors="coerce")
            if x.notna().sum() < 3:
                dropped.append({"term": time_var, "term_type": "numeric", "reason": "too_few_nonmissing_values"})
            elif x.nunique(dropna=True) < 2:
                dropped.append({"term": time_var, "term_type": "numeric", "reason": "single_unique_value"})
            else:
                exog_parts.append(pd.DataFrame({time_var: x.astype(float)}, index=df.index))
                term_meta.append({"column": time_var, "source_term": time_var, "term_type": "sampling_time"})

    # Extra numeric covariates.
    for col in numeric_covariates:
        if col == time_var:
            continue
        if col not in df.columns:
            dropped.append({"term": col, "term_type": "numeric", "reason": "missing_column"})
            continue
        x = pd.to_numeric(df[col], errors="coerce")
        if x.notna().sum() < 3:
            dropped.append({"term": col, "term_type": "numeric", "reason": "too_few_nonmissing_values"})
            continue
        if x.nunique(dropna=True) < 2:
            dropped.append({"term": col, "term_type": "numeric", "reason": "single_unique_value"})
            continue
        # Standardize numeric covariates to improve convergence.
        mean = x.mean()
        sd = x.std(ddof=0)
        if sd and np.isfinite(sd) and sd > 0:
            x2 = (x - mean) / sd
            cname = f"z_{col}"
        else:
            x2 = x
            cname = col
        exog_parts.append(pd.DataFrame({cname: x2.astype(float)}, index=df.index))
        term_meta.append({"column": cname, "source_term": col, "term_type": "numeric_covariate"})

    # Categorical covariates.
    for col in categorical_covariates:
        if col == bioproject_col:
            continue
        if col not in df.columns:
            dropped.append({"term": col, "term_type": "categorical", "reason": "missing_column"})
            continue
        s = clean_categorical_series(df[col], missing_level=missing_level)
        nlev = s.nunique(dropna=False)
        if nlev < 2:
            dropped.append({"term": col, "term_type": "categorical", "reason": "single_level"})
            continue
        if drop_no_within_bp_variation and include_bioproject:
            tmp = pd.DataFrame({bioproject_col: df[bioproject_col].astype(str), col: s})
            nun = tmp.groupby(bioproject_col)[col].nunique(dropna=False)
            if not (nun >= 2).any():
                dropped.append({"term": col, "term_type": "categorical", "reason": "no_within_bioproject_variation_with_bioproject_fixed_effect"})
                continue
        dummies = pd.get_dummies(s, prefix=f"C({col})", drop_first=True, dtype=float)
        if dummies.shape[1] == 0:
            dropped.append({"term": col, "term_type": "categorical", "reason": "no_dummy_columns"})
            continue
        exog_parts.append(dummies)
        for c in dummies.columns:
            term_meta.append({"column": c, "source_term": col, "term_type": "categorical_covariate"})

    X = pd.concat(exog_parts, axis=1)
    X = X.apply(pd.to_numeric, errors="coerce")

    # Drop columns with non-finite values after design construction.
    bad_cols = [c for c in X.columns if not np.isfinite(X[c]).all()]
    if bad_cols:
        for c in bad_cols:
            dropped.append({"term": c, "term_type": "design_column", "reason": "nonfinite_design_values"})
        X = X.drop(columns=bad_cols)
        term_meta = [m for m in term_meta if m["column"] not in set(bad_cols)]

    # Drop zero-variance columns except Intercept.
    zero_var_cols = []
    for c in X.columns:
        if c == "Intercept":
            continue
        if float(np.nanstd(X[c].values)) == 0.0:
            zero_var_cols.append(c)
    if zero_var_cols:
        for c in zero_var_cols:
            dropped.append({"term": c, "term_type": "design_column", "reason": "zero_variance_after_subset"})
        X = X.drop(columns=zero_var_cols)
        term_meta = [m for m in term_meta if m["column"] not in set(zero_var_cols)]

    # QR rank filter.
    rank_before = int(np.linalg.matrix_rank(X.values)) if X.shape[1] else 0
    ncol_before = int(X.shape[1])
    kept_cols = list(X.columns)
    rank_note = "full_rank"

    if X.shape[1] > 0 and rank_before < X.shape[1]:
        try:
            from scipy.linalg import qr
            Q, R, piv = qr(X.values, mode="economic", pivoting=True)
            diag = np.abs(np.diag(R))
            tol = np.finfo(float).eps * max(X.shape) * (diag[0] if len(diag) else 1.0)
            keep_n = int((diag > tol).sum())
            keep_idx = sorted(piv[:keep_n])
            kept_cols = [X.columns[i] for i in keep_idx]
            dropped_cols = [c for c in X.columns if c not in set(kept_cols)]
            for c in dropped_cols:
                dropped.append({"term": c, "term_type": "design_column", "reason": "rank_deficient_dropped_by_qr"})
            X = X[kept_cols].copy()
            term_meta = [m for m in term_meta if m["column"] in set(kept_cols)]
            rank_note = f"rank_filtered_by_qr:{ncol_before}->{X.shape[1]}"
        except Exception as e:
            rank_note = f"rank_deficient_but_qr_failed:{type(e).__name__}:{str(e)}"

    return X, pd.DataFrame(term_meta), pd.DataFrame(dropped), rank_before, ncol_before, rank_note


def fit_one_mixedlm(y, X, groups):
    try:
        import statsmodels.api as sm
    except ImportError as e:
        raise ImportError("statsmodels is required for MixedLM. Install statsmodels in this environment.") from e

    # MixedLM requires aligned arrays.
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        model = sm.MixedLM(endog=y.astype(float), exog=X.astype(float), groups=groups.astype(str))
        last_err = None
        for method in ["lbfgs", "bfgs", "cg", "powell", "nm"]:
            try:
                fit = model.fit(reml=False, method=method, maxiter=1000, disp=False)
                warn_msgs = [str(w.message) for w in wlist]
                return fit, method, "; ".join(warn_msgs)
            except Exception as e:
                last_err = e
                continue
        raise RuntimeError(f"MixedLM failed for all optimizers. Last error: {type(last_err).__name__}: {last_err}")


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


def main():
    ap = argparse.ArgumentParser(description="Fit full LME models for pair-level time-series oral-gut metrics.")
    ap.add_argument("--input", required=True, help="Pair-level metric table, e.g. lme_input_10_metrics.tsv")
    ap.add_argument("--sampling-time", default="00.data/pair_sampling_time_from_pair_time_point.tsv")
    ap.add_argument("--out-dir", default="05.lme_time_series_full_model")
    ap.add_argument("--outcomes", default="delta_alpha_observed,delta_alpha_shannon,delta_alpha_invsimpson,beta_braycurtis_distance,beta_weighted_unifrac_distance,beta_unweighted_unifrac_distance,shared_otu_count,shared_gut_abundance_fraction,shared_otu_count_fraction_of_gut,sourcetracker_paired_oral_fraction")
    ap.add_argument("--time-var", default="sampling_time_z_by_bp")
    ap.add_argument("--bioproject-col", default="bioproject")
    ap.add_argument("--pair-col", default="pair_id")
    ap.add_argument("--subject-col", default="subject_id")
    ap.add_argument("--categorical-covariates", default="age_cleaned_group,sex_cleaned,status_HC_vs_disease,status_cleaned,bmi_cleaned_group,continent_group,sequencing_region_cleaned")
    ap.add_argument("--numeric-covariates", default="")
    ap.add_argument("--include-bioproject", type=parse_bool, default=True)
    ap.add_argument("--drop-no-within-bp-variation", type=parse_bool, default=True)
    ap.add_argument("--missing-level", default="Missing")
    ap.add_argument("--min-n", type=int, default=30)
    ap.add_argument("--min-groups", type=int, default=5)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    input_fp = Path(args.input).resolve()
    time_fp = Path(args.sampling_time).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "merged_input": out_dir / "lme_pair_level_input_with_sampling_time.tsv",
        "fixed_effects": out_dir / "lme_full_model_fixed_effects.tsv",
        "time_effects": out_dir / "lme_full_model_sampling_time_effects.tsv",
        "diagnostics": out_dir / "lme_full_model_diagnostics.tsv",
        "dropped_terms": out_dir / "lme_full_model_dropped_terms.tsv",
        "summary_json": out_dir / "lme_full_model_summary.json",
    }
    for p in outputs.values():
        if p.exists() and not args.overwrite:
            raise FileExistsError(f"Output exists: {p}. Use --overwrite.")

    print("Resolved paths:", flush=True)
    print(f"  input:         {input_fp}", flush=True)
    print(f"  sampling_time: {time_fp}", flush=True)
    print(f"  out_dir:       {out_dir}", flush=True)

    df = pd.read_csv(input_fp, sep="\t", dtype=str)
    time = pd.read_csv(time_fp, sep="\t", dtype=str)

    required_input = [args.bioproject_col, args.pair_col]
    missing = [c for c in required_input if c not in df.columns]
    if missing:
        raise ValueError("Input missing required columns: " + ", ".join(missing))

    required_time = ["bioproject", "pair_id", "subject_id", args.time_var, "subject_block", "pair_block"]
    missing_time = [c for c in required_time if c not in time.columns]
    if missing_time:
        raise ValueError("Sampling-time table missing required columns: " + ", ".join(missing_time))

    # Merge pair-level time variables.
    time_cols = [
        "bioproject", "pair_id", "subject_id", "pair_time_point",
        "sampling_time_numeric", "sampling_time_parse_type", "sampling_time_scale_type",
        "sampling_time_is_parsed", "sampling_time_z_by_bp", "sampling_time_rank_by_bp",
        "sampling_time_rank_z_by_bp", "subject_block", "pair_block",
    ]
    time_cols = [c for c in time_cols if c in time.columns]

    merged = df.merge(
        time[time_cols],
        left_on=[args.bioproject_col, args.pair_col],
        right_on=["bioproject", "pair_id"],
        how="left",
        suffixes=("", "_time"),
    )
    for c in ["bioproject_time", "pair_id_time"]:
        if c in merged.columns:
            merged = merged.drop(columns=[c])

    # Prefer subject_block from time table; if missing, build from input subject_id.
    if "subject_block" not in merged.columns or merged["subject_block"].isna().all():
        subj = merged[args.subject_col] if args.subject_col in merged.columns else merged.get("subject_id", pd.Series(np.nan, index=merged.index))
        merged["subject_block"] = merged[args.bioproject_col].astype(str) + "__" + subj.astype(str)

    outcomes = choose_existing_columns(merged, parse_list(args.outcomes))
    missing_outcomes = [c for c in parse_list(args.outcomes) if c not in merged.columns]
    if missing_outcomes:
        print("WARNING: requested outcomes missing and skipped: " + ", ".join(missing_outcomes), flush=True)
    if not outcomes:
        raise ValueError("No requested outcomes were found in input table.")

    categorical_covariates = choose_existing_columns(merged, parse_list(args.categorical_covariates))
    missing_cat = [c for c in parse_list(args.categorical_covariates) if c not in merged.columns]
    if missing_cat:
        print("WARNING: requested categorical covariates missing and skipped: " + ", ".join(missing_cat), flush=True)

    numeric_covariates = choose_existing_columns(merged, parse_list(args.numeric_covariates))
    missing_num = [c for c in parse_list(args.numeric_covariates) if c not in merged.columns]
    if missing_num:
        print("WARNING: requested numeric covariates missing and skipped: " + ", ".join(missing_num), flush=True)

    merged.to_csv(outputs["merged_input"], sep="\t", index=False)

    fixed_rows = []
    diag_rows = []
    dropped_rows = []

    print("\nFitting models:", flush=True)

    for outcome in outcomes:
        print(f"  outcome: {outcome}", flush=True)
        dat = merged.copy()
        dat["__y__"] = pd.to_numeric(dat[outcome], errors="coerce")
        dat[args.time_var] = pd.to_numeric(dat[args.time_var], errors="coerce")
        dat["subject_block"] = dat["subject_block"].map(clean_missing)
        dat[args.bioproject_col] = dat[args.bioproject_col].map(clean_missing)

        # Required complete cases: y, time, subject_block, BioProject.
        before_n = dat.shape[0]
        dat = dat.loc[
            dat["__y__"].notna()
            & dat[args.time_var].notna()
            & dat["subject_block"].notna()
            & dat[args.bioproject_col].notna(),
            :,
        ].copy()
        after_required_n = dat.shape[0]

        if after_required_n < args.min_n:
            diag_rows.append({
                "outcome": outcome,
                "fit_status": "skipped",
                "skip_reason": f"too_few_complete_required_rows:{after_required_n}<min_n:{args.min_n}",
                "n_rows_input": before_n,
                "n_rows_model": after_required_n,
            })
            continue

        n_groups = dat["subject_block"].nunique()
        if n_groups < args.min_groups:
            diag_rows.append({
                "outcome": outcome,
                "fit_status": "skipped",
                "skip_reason": f"too_few_subject_groups:{n_groups}<min_groups:{args.min_groups}",
                "n_rows_input": before_n,
                "n_rows_model": after_required_n,
                "n_subject_groups": n_groups,
            })
            continue

        X, term_meta, dropped, rank_before, ncol_before, rank_note = make_fixed_effect_matrix(
            dat,
            outcome=outcome,
            time_var=args.time_var,
            bioproject_col=args.bioproject_col,
            categorical_covariates=categorical_covariates,
            numeric_covariates=numeric_covariates,
            include_bioproject=args.include_bioproject,
            missing_level=args.missing_level,
            drop_no_within_bp_variation=args.drop_no_within_bp_variation,
        )

        if not dropped.empty:
            dropped = dropped.copy()
            dropped.insert(0, "outcome", outcome)
            dropped_rows.append(dropped)

        # Align y/groups to X.
        dat2 = dat.loc[X.index].copy()
        y = dat2["__y__"].astype(float)
        groups = dat2["subject_block"].astype(str)

        if X.shape[1] < 2:
            diag_rows.append({
                "outcome": outcome,
                "fit_status": "skipped",
                "skip_reason": "too_few_fixed_effect_columns_after_filtering",
                "n_rows_input": before_n,
                "n_rows_model": dat2.shape[0],
                "n_subject_groups": groups.nunique(),
                "rank_note": rank_note,
            })
            continue

        try:
            fit, optimizer, warn_msg = fit_one_mixedlm(y, X, groups)
            params = fit.params
            bse = fit.bse
            pvalues = fit.pvalues
            conf = fit.conf_int()

            # Random intercept variance/residual variance.
            try:
                random_intercept_var = float(fit.cov_re.iloc[0, 0])
            except Exception:
                random_intercept_var = np.nan
            residual_var = float(fit.scale) if hasattr(fit, "scale") else np.nan

            meta_map = term_meta.set_index("column").to_dict(orient="index") if not term_meta.empty else {}

            for term in X.columns:
                fixed_rows.append({
                    "outcome": outcome,
                    "term": term,
                    "source_term": meta_map.get(term, {}).get("source_term", term),
                    "term_type": meta_map.get(term, {}).get("term_type", "unknown"),
                    "estimate": float(params.get(term, np.nan)),
                    "std_error": float(bse.get(term, np.nan)),
                    "z_value": float(params.get(term, np.nan) / bse.get(term, np.nan)) if np.isfinite(bse.get(term, np.nan)) and bse.get(term, np.nan) != 0 else np.nan,
                    "p_value": float(pvalues.get(term, np.nan)),
                    "p_significance_label": p_to_stars(pvalues.get(term, np.nan)),
                    "ci_lower": float(conf.loc[term, 0]) if term in conf.index else np.nan,
                    "ci_upper": float(conf.loc[term, 1]) if term in conf.index else np.nan,
                    "n_rows_model": int(dat2.shape[0]),
                    "n_subject_groups": int(groups.nunique()),
                    "n_bioprojects": int(dat2[args.bioproject_col].nunique()),
                    "optimizer": optimizer,
                    "model": f"{outcome} ~ BioProject + {args.time_var} + covariates + (1|subject_block)",
                })

            diag_rows.append({
                "outcome": outcome,
                "fit_status": "ok",
                "skip_reason": "",
                "n_rows_input": int(before_n),
                "n_rows_after_required_complete_case": int(after_required_n),
                "n_rows_model": int(dat2.shape[0]),
                "n_rows_dropped_required_complete_case": int(before_n - after_required_n),
                "n_subject_groups": int(groups.nunique()),
                "n_bioprojects": int(dat2[args.bioproject_col].nunique()),
                "n_fixed_effect_columns_before_rank_filter": int(ncol_before),
                "rank_before_filter": int(rank_before),
                "n_fixed_effect_columns_model": int(X.shape[1]),
                "rank_model": int(np.linalg.matrix_rank(X.values)),
                "rank_note": rank_note,
                "converged": bool(getattr(fit, "converged", False)),
                "optimizer": optimizer,
                "log_likelihood": float(fit.llf),
                "aic": float(fit.aic) if hasattr(fit, "aic") else np.nan,
                "bic": float(fit.bic) if hasattr(fit, "bic") else np.nan,
                "random_intercept_variance": random_intercept_var,
                "residual_variance": residual_var,
                "warnings": warn_msg,
                "time_var": args.time_var,
                "include_bioproject": bool(args.include_bioproject),
            })

        except Exception as e:
            diag_rows.append({
                "outcome": outcome,
                "fit_status": "failed",
                "skip_reason": f"{type(e).__name__}:{str(e)}",
                "n_rows_input": int(before_n),
                "n_rows_after_required_complete_case": int(after_required_n),
                "n_rows_model": int(dat2.shape[0]) if 'dat2' in locals() else np.nan,
                "n_subject_groups": int(groups.nunique()) if 'groups' in locals() else np.nan,
                "n_bioprojects": int(dat[args.bioproject_col].nunique()) if args.bioproject_col in dat.columns else np.nan,
                "time_var": args.time_var,
            })
            continue

    fixed_df = pd.DataFrame(fixed_rows)
    diag_df = pd.DataFrame(diag_rows)
    dropped_df = pd.concat(dropped_rows, ignore_index=True) if dropped_rows else pd.DataFrame(columns=["outcome", "term", "term_type", "reason"])

    if not fixed_df.empty:
        fixed_df["q_value_all_fixed_effects"] = bh_fdr(fixed_df["p_value"])
        fixed_df["q_significance_label_all_fixed_effects"] = fixed_df["q_value_all_fixed_effects"].map(p_to_stars)
        time_mask = fixed_df["term"].eq(args.time_var)
        fixed_df["q_value_sampling_time_effects"] = np.nan
        if time_mask.any():
            fixed_df.loc[time_mask, "q_value_sampling_time_effects"] = bh_fdr(fixed_df.loc[time_mask, "p_value"])
        fixed_df["q_significance_label_sampling_time_effects"] = fixed_df["q_value_sampling_time_effects"].map(p_to_stars)

    time_df = fixed_df.loc[fixed_df["term"].eq(args.time_var), :].copy() if not fixed_df.empty else pd.DataFrame()

    fixed_df.to_csv(outputs["fixed_effects"], sep="\t", index=False)
    time_df.to_csv(outputs["time_effects"], sep="\t", index=False)
    diag_df.to_csv(outputs["diagnostics"], sep="\t", index=False)
    dropped_df.to_csv(outputs["dropped_terms"], sep="\t", index=False)

    summary = {
        "input": str(input_fp),
        "sampling_time": str(time_fp),
        "out_dir": str(out_dir),
        "outcomes_requested": parse_list(args.outcomes),
        "outcomes_fitted_or_attempted": outcomes,
        "time_var": args.time_var,
        "categorical_covariates_requested": parse_list(args.categorical_covariates),
        "categorical_covariates_available": categorical_covariates,
        "numeric_covariates_requested": parse_list(args.numeric_covariates),
        "numeric_covariates_available": numeric_covariates,
        "include_bioproject": bool(args.include_bioproject),
        "random_effect": "random intercept for subject_block = BioProject__subject_id",
        "n_fixed_effect_rows": int(fixed_df.shape[0]),
        "n_time_effect_rows": int(time_df.shape[0]),
        "n_diagnostics_rows": int(diag_df.shape[0]),
        "outputs": {k: str(v) for k, v in outputs.items()},
    }
    with open(outputs["summary_json"], "w") as f:
        json.dump(summary, f, indent=2)

    print("\nDone.", flush=True)
    print("Outputs:", flush=True)
    for k, v in outputs.items():
        print(f"  {k}: {v}", flush=True)

    if not time_df.empty:
        show_cols = ["outcome", "term", "estimate", "std_error", "z_value", "p_value", "q_value_sampling_time_effects", "n_rows_model", "n_subject_groups", "converged"]
        show_cols = [c for c in show_cols if c in time_df.columns]
        print("\nSampling time effects:", flush=True)
        print(time_df[show_cols].to_string(index=False), flush=True)
    else:
        print("\nNo sampling time effect rows were produced.", flush=True)


if __name__ == "__main__":
    main()

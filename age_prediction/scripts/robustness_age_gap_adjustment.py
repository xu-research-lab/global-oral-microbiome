#!/usr/bin/env python
# =====================================================================
# Robustness checks for the disease excess-age-gap figure
# (plot_age_acceleration_1row.py). The main figure reports the unadjusted
# median age-gap difference vs HC, which mixes (i) disease effects,
# (ii) age-distribution differences between groups and (iii) study-
# specific (batch/sequencing) effects.  This script recomputes the same
# site x status contrasts with three adjustments and writes a
# side-by-side comparison table:
#
#   1. agestd   - direct standardization of the HC age composition to the
#                 disease group's (10-year strata, median gap per stratum,
#                 bootstrap CI) -> removes age-composition confounding.
#   2. ols_age  - OLS gap ~ status + age + age^2 with HC1 robust SE
#                 (mean-based, continuous age control).
#   3. instudy  - within-study HC reference: for disease studies that
#                 also contain >=5 same-site HC samples, the gap
#                 difference is computed inside the same study, so
#                 sequencing/processing batch cancels out.  Statuses
#                 living in disease-only studies are not identifiable
#                 here and get NaN (study effect and status effect are
#                 perfectly collinear for them; said explicitly in the
#                 "single_study" / notes columns).
#
# Output (next to the existing figure stats):
#   result/figures/age_gap_robustness_stats.csv           (one row per site x status)
#   result/figures/age_gap_robustness_instudy_detail.csv  (per-study within-study diffs)
# =====================================================================
from pathlib import Path
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

BASE = Path(os.environ.get("GLOBAL_ORAL_ROOT") or Path(__file__).resolve().parents[2])
RESULT = Path(__file__).resolve().parents[1] / "result"
DATA_DIR = RESULT / "leave_studies_out_rf"
METADATA_PATH = (BASE / "data" / "metadata-and-biom-v2-qc"
                 / "metadata_oral_qc_dedup.tsv")
OUT_DIR = RESULT / "figures"

ORAL_SITES = ["saliva", "tooth", "tongue", "oral mucosa"]
MIN_N = 10          # same filter as the original figure
MIN_N_STRAT = 5     # min HC samples per age stratum for standardization
MIN_N_STRAT_DIS = 3 # min disease samples per age stratum
N_BOOT = 1000
SEED = 42


def decade_bin(age):
    d = (np.floor(np.asarray(age, float) / 10) * 10)
    return np.clip(d, 0, 70)  # 70+ as one open-top bin


def ols_hc1(y, X):
    """OLS with HC1 heteroskedastic-robust SE. Returns coefs, ses, p (z-based)."""
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    resid = y - X @ beta
    n, p = X.shape
    meat = X.T @ (X * (resid ** 2)[:, None])
    V = XtX_inv @ meat @ XtX_inv * (n / max(n - p, 1))
    se = np.sqrt(np.clip(np.diag(V), 0, None))
    from scipy.stats import norm
    z = np.divide(beta, se, out=np.zeros_like(beta), where=se > 0)
    pv = 2 * norm.sf(np.abs(z))
    return beta, se, pv


def main():
    rng = np.random.default_rng(SEED)
    metadata = pd.read_csv(METADATA_PATH, sep="\t", dtype={"sample_id": str})

    rows, detail_rows = [], []
    for j in ORAL_SITES:
        files = sorted(glob.glob(str(DATA_DIR /
                f"age_predict_test_*_{j}_genus_range_1_80.csv")))
        if not files:
            print(f"[{j}] no prediction CSVs found, skipping")
            continue
        preds = pd.concat([pd.read_csv(f, dtype={"sample_id": str}) for f in files],
                          ignore_index=True)
        meta_pick = metadata.loc[metadata.site_cleaned == j].set_index("sample_id")
        sid = np.intersect1d(meta_pick.index.values, preds.sample_id.values)
        if len(sid) == 0:
            print(f"[{j}] no sample overlap with metadata, skipping")
            continue
        df = preds.set_index("sample_id").loc[sid].copy()
        df["status"] = meta_pick.loc[sid, "status_cleaned"].values
        df["bio"] = meta_pick.loc[sid, "bioproject_id"].astype(str).values
        df["gap"] = df.pred_age.values - df.true_age.values

        counts = df["status"].value_counts()
        keep = [s for s in counts.index if s == "HC" or counts[s] >= MIN_N]
        df = df[df["status"].isin(keep)]
        hc_all = df.loc[df.status == "HC"]
        hc_med = np.median(hc_all.gap.values)

        for s in [x for x in keep if x != "HC"]:
            dis = df.loc[df.status == s]
            n_d, n_hc = len(dis), len(hc_all)
            unadj = np.median(dis.gap.values) - hc_med

            # ---- 1. direct age standardization (HC composition -> disease's)
            d_bin, h_bin = decade_bin(dis.true_age.values), decade_bin(hc_all.true_age.values)
            strata = np.intersect1d(np.unique(d_bin), np.unique(h_bin))
            w, md, mh = [], [], []
            cov_n = 0
            for t in strata:
                d_t, h_t = dis.gap.values[d_bin == t], hc_all.gap.values[h_bin == t]
                if len(d_t) >= MIN_N_STRAT_DIS and len(h_t) >= MIN_N_STRAT:
                    w.append(len(d_t)); md.append(np.median(d_t)); mh.append(np.median(h_t))
                    cov_n += len(d_t)
            w, md, mh = np.array(w, float), np.array(md), np.array(mh)
            if len(w) == 0:
                astd = ast_lo = ast_hi = np.nan
                cover = 0.0
            else:
                cover = cov_n / n_d
                astd = float(np.sum(w / w.sum() * (md - mh)))
                # per-stratum sample arrays for bootstrap resampling
                strata_dis = {t: dis.gap.values[d_bin == t] for t in strata}
                strata_hc = {t: hc_all.gap.values[h_bin == t] for t in strata}
                used = [t for t in strata
                        if len(strata_dis[t]) >= MIN_N_STRAT_DIS
                        and len(strata_hc[t]) >= MIN_N_STRAT]
                boots = []
                for _ in range(N_BOOT):
                    num = np.empty(len(used))
                    for k, t in enumerate(used):
                        dd = rng.choice(strata_dis[t], len(strata_dis[t]), True)
                        hh = rng.choice(strata_hc[t], len(strata_hc[t]), True)
                        num[k] = np.median(dd) - np.median(hh)
                    ww = np.array([np.sum(d_bin == t) for t in used], float)
                    boots.append(np.sum(ww / ww.sum() * num))
                ast_lo, ast_hi = np.percentile(boots, [2.5, 97.5])

            # ---- 2. OLS gap ~ status + age + age^2
            sub = df  # all kept statuses, HC as reference
            y = sub.gap.values.astype(float)
            ages = sub.true_age.values.astype(float)
            cols = [np.ones(len(y)), ages, ages ** 2]
            names = ["HC(ref)", "age", "age2"]
            for st in [x for x in keep if x != "HC"]:
                cols.append((sub.status.values == st).astype(float))
                names.append(st)
            X = np.column_stack(cols)
            beta, se, pv = ols_hc1(y, X)
            idx = names.index(s)
            ols_c, ols_se, ols_p = beta[idx], se[idx], pv[idx]
            from scipy.stats import norm
            z = norm.ppf(.975)
            ols_lo, ols_hi = ols_c - z * ols_se, ols_c + z * ols_se

            # ---- 3. within-study HC reference
            hc_counts = hc_all.groupby("bio").size()
            qual = [b for b in dis.bio.unique() if hc_counts.get(b, 0) >= MIN_N_STRAT]
            if qual:
                d_in = dis.loc[dis.bio.isin(qual)]
                h_in = hc_all.loc[hc_all.bio.isin(qual)]
                instudy = np.median(d_in.gap.values) - np.median(h_in.gap.values)
                in_n = len(d_in)
                for b in qual:
                    db = dis.loc[dis.bio == b]
                    detail_rows.append(dict(site=j, status=s, bioproject=b,
                        n_dis=len(db), n_hc=int(hc_counts[b]),
                        med_gap_dis=float(np.median(db.gap.values)),
                        med_gap_hc=float(np.median(hc_all.loc[hc_all.bio == b, "gap"].values)),
                        diff=float(np.median(db.gap.values) -
                                   np.median(hc_all.loc[hc_all.bio == b, "gap"].values))))
                in_note = ""
            else:
                instudy, in_n, in_note = np.nan, 0, \
                    "no disease study contains >=5 same-site HC"

            rows.append(dict(
                site=j, status=s, n=n_d, n_studies=dis.bio.nunique(),
                single_study=dis.bio.nunique() == 1,
                med_age_dis=float(np.median(dis.true_age)),
                med_age_hc=float(np.median(hc_all.true_age)),
                hc_med_gap=float(hc_med),
                unadj_delta=float(unadj),
                astd_delta=astd, astd_lo=ast_lo, astd_hi=ast_hi,
                astd_coverage=cover,
                ols_age_adj=ols_c, ols_lo=ols_lo, ols_hi=ols_hi, ols_p=ols_p,
                instudy_delta=instudy, instudy_n=in_n, instudy_note=in_note,
            ))

    stats = pd.DataFrame(rows)
    out = OUT_DIR / "age_gap_robustness_stats.csv"
    stats.to_csv(out, index=False)
    pd.DataFrame(detail_rows).to_csv(
        OUT_DIR / "age_gap_robustness_instudy_detail.csv", index=False)

    pd.set_option("display.width", 200)
    for j in ORAL_SITES:
        sub = stats[stats.site == j]
        if sub.empty:
            continue
        print(f"\n===== {j} =====")
        cols = ["status", "n", "n_studies", "med_age_dis",
                "unadj_delta", "astd_delta", "astd_lo", "astd_hi",
                "astd_coverage", "ols_age_adj", "instudy_delta", "instudy_n"]
        print(sub[cols].round(2).to_string(index=False))
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
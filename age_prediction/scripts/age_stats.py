"""Disease age-gap comparisons with healthy controls."""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from age_common import ORAL_SITES, RESULT, site_predictions, status_groups

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

def boot_ci_median_diff(dis, hc, rng, n_boot=2000):
    d = np.array([np.median(rng.choice(dis, len(dis), True)) -
                  np.median(rng.choice(hc, len(hc), True)) for _ in range(n_boot)])
    return np.percentile(d, [2.5, 97.5])

def stars(q):
    return "***" if q < .001 else "**" if q < .01 else "*" if q < .05 else "ns"


def disease_age_gap_stats(metadata, data_dir=None, min_n=10, n_boot=2000):
    data_dir = RESULT / "leave_studies_out_rf" if data_dir is None else data_dir
    stat_rows = []
    rng = np.random.default_rng(0)
    for site in ORAL_SITES:
        frame = site_predictions(site, metadata, data_dir, sort_ids=True).rename(
            columns={"status_cleaned": "status"})
        frame, keep = status_groups(frame, min_n)
        hc = frame.loc[frame.status == "HC", "delta_age"].to_numpy()
        if not len(hc):
            continue
        rows = []
        for status in [s for s in keep if s != "HC"]:
            values = frame.loc[frame.status == status, "delta_age"].to_numpy()
            p = mannwhitneyu(values, hc, alternative="two-sided")[1]
            lo, hi = boot_ci_median_diff(values, hc, rng, n_boot)
            rows.append(dict(site=site, status=status, n=len(values), med=np.median(values),
                             delta=np.median(values) - np.median(hc), ci_lo=lo, ci_hi=hi,
                             cliff=cliffs_delta(values, hc), p=p))
        if rows:
            stats = pd.DataFrame(rows)
            stats["fdr"] = bh_fdr(stats["p"].to_numpy())
            stats["hc_med"] = np.median(hc)
            stats["hc_n"] = len(hc)
            stat_rows.append(stats)
    columns = ["site", "status", "n", "med", "delta", "ci_lo", "ci_hi", "cliff", "p", "fdr", "hc_med", "hc_n"]
    stats = pd.concat(stat_rows, ignore_index=True) if stat_rows else pd.DataFrame(columns=columns)
    stats["sig"] = stats["fdr"].map(stars)
    return stats

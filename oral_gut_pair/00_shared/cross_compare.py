#!/usr/bin/env python3
"""Compare harmonized pair metrics across cross-sectional and longitudinal designs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

from merge_pair_metrics import METRICS


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().sort_values()
    if valid.empty:
        return result
    n = valid.shape[0]
    adjusted = valid * n / np.arange(1, n + 1)
    adjusted = adjusted.iloc[::-1].cummin().iloc[::-1].clip(upper=1)
    result.loc[adjusted.index] = adjusted
    return result


def cliffs_delta(left: np.ndarray, right: np.ndarray) -> float:
    if left.size == 0 or right.size == 0:
        return np.nan
    # Rank-based identity: delta = 2U/(n1*n2)-1.
    u = mannwhitneyu(left, right, alternative="two-sided").statistic
    return float(2 * u / (left.size * right.size) - 1)


def describe(values: pd.Series) -> dict[str, float | int]:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    return {
        "n": int(numeric.shape[0]),
        "mean": float(numeric.mean()) if not numeric.empty else np.nan,
        "sd": float(numeric.std()) if numeric.shape[0] > 1 else np.nan,
        "median": float(numeric.median()) if not numeric.empty else np.nan,
        "q25": float(numeric.quantile(0.25)) if not numeric.empty else np.nan,
        "q75": float(numeric.quantile(0.75)) if not numeric.empty else np.nan,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(input_path, sep="\t")
    required = {"analysis_set", "bioproject", "pair_id", *METRICS}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Input missing columns: {missing}")

    descriptive_rows = []
    comparison_rows = []
    for metric in METRICS:
        for analysis_set, subset in table.groupby("analysis_set"):
            row = {"metric": metric, "analysis_set": analysis_set}
            row.update(describe(subset[metric]))
            row["n_bioproject"] = int(
                subset.loc[subset[metric].notna(), "bioproject"].nunique()
            )
            descriptive_rows.append(row)

        single = pd.to_numeric(
            table.loc[table["analysis_set"].eq("single_pair"), metric],
            errors="coerce",
        ).dropna()
        time = pd.to_numeric(
            table.loc[table["analysis_set"].eq("time_series"), metric],
            errors="coerce",
        ).dropna()
        if single.empty or time.empty:
            u_stat = p_between = delta = np.nan
        else:
            test = mannwhitneyu(single, time, alternative="two-sided")
            u_stat = float(test.statistic)
            p_between = float(test.pvalue)
            delta = cliffs_delta(single.to_numpy(), time.to_numpy())

        row = {
            "metric": metric,
            "n_single_pair": int(single.shape[0]),
            "n_time_series": int(time.shape[0]),
            "median_single_pair": float(single.median()) if not single.empty else np.nan,
            "median_time_series": float(time.median()) if not time.empty else np.nan,
            "mann_whitney_u": u_stat,
            "p_between_designs": p_between,
            "cliffs_delta_single_minus_time": delta,
            "direction_agrees": bool(
                not single.empty
                and not time.empty
                and np.sign(single.median()) == np.sign(time.median())
            ),
        }
        for analysis_set, values in [("single_pair", single), ("time_series", time)]:
            nonzero = values.loc[values.ne(0)]
            if nonzero.empty:
                statistic = p_value = np.nan
            else:
                test = wilcoxon(nonzero, alternative="two-sided")
                statistic, p_value = float(test.statistic), float(test.pvalue)
            row[f"wilcoxon_stat_{analysis_set}"] = statistic
            row[f"p_vs_zero_{analysis_set}"] = p_value
        comparison_rows.append(row)

    descriptives = pd.DataFrame(descriptive_rows)
    comparisons = pd.DataFrame(comparison_rows)
    comparisons["q_between_designs"] = bh_adjust(comparisons["p_between_designs"])
    for analysis_set in ("single_pair", "time_series"):
        comparisons[f"q_vs_zero_{analysis_set}"] = bh_adjust(
            comparisons[f"p_vs_zero_{analysis_set}"]
        )
    comparisons["concordance_class"] = np.select(
        [
            comparisons["direction_agrees"]
            & comparisons["cliffs_delta_single_minus_time"].abs().lt(0.33),
            comparisons["direction_agrees"],
        ],
        ["agree_small_design_effect", "agree_with_design_shift"],
        default="direction_disagrees",
    )

    descriptives_path = out_dir / "metric_descriptives.tsv"
    comparisons_path = out_dir / "metric_concordance.tsv"
    summary_path = out_dir / "cross_comparison_summary.json"
    descriptives.to_csv(descriptives_path, sep="\t", index=False)
    comparisons.to_csv(comparisons_path, sep="\t", index=False)
    summary = {
        "input": str(input_path),
        "n_pairs": int(table.shape[0]),
        "n_bioprojects": int(table["bioproject"].nunique()),
        "metrics": METRICS,
        "n_metrics_agree_small_design_effect": int(
            comparisons["concordance_class"].eq("agree_small_design_effect").sum()
        ),
        "outputs": {
            "descriptives": str(descriptives_path),
            "concordance": str(comparisons_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

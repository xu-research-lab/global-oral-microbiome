#!/usr/bin/env python3
import datetime
import itertools
import os
import re
import sys

import pandas as pd
import skbio
from skbio import DistanceMatrix
from skbio.stats.distance import permanova


def log_msg(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)


def parse_version(version_str):
    parts = re.findall(r"\d+", version_str)
    return tuple(int(part) for part in parts[:3])


def require_minimum_versions():
    min_python = (3, 8, 0)
    min_pandas = (1, 0, 0)
    min_skbio = (0, 5, 8)

    python_version = sys.version_info[:3]
    pandas_version = parse_version(pd.__version__)
    skbio_version = parse_version(skbio.__version__)

    log_msg(
        "Runtime versions | "
        f"Python={sys.version.split()[0]} | "
        f"pandas={pd.__version__} | "
        f"scikit-bio={skbio.__version__}"
    )
    log_msg(
        "Minimum supported versions | "
        f"Python={'.'.join(map(str, min_python))} | "
        f"pandas={'.'.join(map(str, min_pandas))} | "
        f"scikit-bio={'.'.join(map(str, min_skbio))}"
    )

    errors = []
    if python_version < min_python:
        errors.append(
            f"Python {sys.version.split()[0]} is too old; need >= "
            f"{'.'.join(map(str, min_python))}"
        )
    if pandas_version < min_pandas:
        errors.append(
            f"pandas {pd.__version__} is too old; need >= "
            f"{'.'.join(map(str, min_pandas))}"
        )
    if skbio_version < min_skbio:
        errors.append(
            f"scikit-bio {skbio.__version__} is too old; need >= "
            f"{'.'.join(map(str, min_skbio))}"
        )

    if errors:
        raise RuntimeError("Environment version check failed: " + "; ".join(errors))


def clean_metadata(meta):
    meta = meta.loc[:, ~meta.columns.duplicated()].copy()

    if "bioproject_id" in meta.columns and "bioproject" not in meta.columns:
        meta = meta.rename(columns={"bioproject_id": "bioproject"})

    if "age_cleaned_group" in meta.columns:
        meta["age_cleaned_group"] = (
            meta["age_cleaned_group"]
            .str.replace("＞", ">", regex=False)
            .str.replace("＜", "<", regex=False)
        )

    if "bmi_cleaned_group" in meta.columns:
        meta["bmi_cleaned_group"] = (
            meta["bmi_cleaned_group"]
            .str.replace("＞", ">", regex=False)
            .str.replace("＜", "<", regex=False)
        )

    if "sample_id" not in meta.columns:
        raise KeyError("Metadata is missing required column: sample_id")

    return meta.set_index("sample_id", drop=True)


def load_distance_matrix(dm_path):
    try:
        return DistanceMatrix.read(dm_path, format="lsmat")
    except Exception as native_error:
        log_msg(
            f"Failed to read distance matrix with skbio lsmat parser. "
            f"Trying pandas fallback... Error: {native_error}"
        )

        df_dm = pd.read_csv(dm_path, sep="\t", index_col=0, dtype=float, engine="c")
        df_dm.index = df_dm.index.astype(str)
        df_dm.columns = df_dm.columns.astype(str)

        if len(df_dm.columns) != len(df_dm.index):
            raise ValueError(
                f"Distance matrix is not square after fallback load: {df_dm.shape}"
            )

        if list(df_dm.columns) != list(df_dm.index):
            if set(df_dm.columns) != set(df_dm.index):
                raise ValueError(
                    "Distance matrix row IDs and column IDs do not match in fallback load."
                )
            df_dm = df_dm.loc[df_dm.index, df_dm.index]

        dm = DistanceMatrix(df_dm.to_numpy(), ids=df_dm.index.tolist())
        del df_dm
        return dm


def compute_r2_from_permanova(result):
    f_value = float(result["test statistic"])
    n = int(result["sample size"])
    g = int(result["number of groups"])

    denominator = (f_value * (g - 1)) + (n - g)
    if denominator <= 0:
        return float("nan")

    return (f_value * (g - 1)) / denominator


def append_result(out_file, metric_name, var, g1, g2, sample_size, r2, p_value, f_stat):
    with open(out_file, "a", encoding="utf-8") as handle:
        handle.write(
            f"{metric_name}\t{var}\t{g1}\t{g2}\t{sample_size}\t"
            f"{r2:.4f}\t{p_value:.4f}\t{f_stat:.4f}\n"
        )


def main():
    require_minimum_versions()

    meta_path = "/home/hexincheng/projects/hxc/oral_project/data/metadata_oral_qc_dedup.tsv"
    dm_path = "/beegfs/dongbiao/oral/beta_oral_qc_dedup/unifrac/unifrac_unweight.dm"
    out_path = (
        "/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/"
        "beta_pcoa_out_20260909/unifrac_unweight/"
        "02-unifrac-unweight-distance-matrix-permanova-20260909.tsv"
    )

    variables = [
        "age_cleaned_group",
        "bmi_cleaned_group",
        "sex_cleaned",
        "platform_cleaned",
        "sequencing_region_cleaned",
        "continent_group",
        "site_cleaned",
    ]

    permutations = 999
    min_group_size = 2
    min_total_size = 5
    success_count = 0
    error_count = 0
    skipped_count = 0
    metric_name = "Unweighted_UniFrac"

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(
            "Metric\tVariable\tGroup1\tGroup2\tSample_Size\tR2\tp_value\tF_statistic\n"
        )

    log_msg("Reading metadata...")
    log_msg(
        f"Run configuration | metric={metric_name} | "
        f"variables={variables} | permutations={permutations} | output={out_path}"
    )
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    meta = clean_metadata(meta)

    log_msg(f"=== Starting processing for {metric_name} ===")
    log_msg(f"Loading {dm_path} into memory. This may take substantial RAM and time.")
    dm = load_distance_matrix(dm_path)
    log_msg(f"Matrix loaded successfully. Dimension: {dm.shape}")

    meta_ids = set(meta.index)
    common_ids = [sample_id for sample_id in dm.ids if sample_id in meta_ids]
    log_msg(
        f"Found {len(common_ids)} intersecting samples between metadata and distance matrix."
    )

    if not common_ids:
        raise RuntimeError(f"No intersecting samples found for {metric_name}.")

    sub_meta = meta.loc[common_ids].copy()

    for var in variables:
        if var not in sub_meta.columns:
            log_msg(f"Variable {var} not found in metadata. Skipping.")
            skipped_count += 1
            continue

        log_msg(f"  -> Processing combinations for variable: {var}")
        valid_meta = sub_meta.dropna(subset=[var]).copy()
        valid_meta = valid_meta[valid_meta[var] != ""]

        groups = valid_meta[var].unique().tolist()
        pairs = list(itertools.combinations(groups, 2))

        if not pairs:
            log_msg(f"    No valid pairwise comparisons for {var}.")
            skipped_count += 1
            continue

        for g1, g2 in pairs:
            pair_meta = valid_meta[valid_meta[var].isin([g1, g2])].copy()
            group_sizes = pair_meta[var].value_counts()
            pair_ids = pair_meta.index.tolist()

            if len(group_sizes) < 2:
                log_msg(f"    Skipping {g1} vs {g2} (fewer than 2 non-empty groups).")
                skipped_count += 1
                continue

            if group_sizes.min() < min_group_size:
                log_msg(
                    f"    Skipping {g1} vs {g2} "
                    f"(a group has fewer than {min_group_size} samples: {group_sizes.to_dict()})."
                )
                skipped_count += 1
                continue

            if len(pair_ids) < min_total_size:
                log_msg(
                    f"    Skipping {g1} vs {g2} "
                    f"(total sample size < {min_total_size}: {len(pair_ids)})."
                )
                skipped_count += 1
                continue

            try:
                pair_dm = dm.filter(pair_ids)
                result = permanova(
                    pair_dm,
                    pair_meta,
                    column=var,
                    permutations=permutations,
                )

                sample_size = int(result["sample size"])
                f_stat = float(result["test statistic"])
                p_value = float(result["p-value"])
                r2 = compute_r2_from_permanova(result)

                append_result(
                    out_path,
                    metric_name,
                    var,
                    g1,
                    g2,
                    sample_size,
                    r2,
                    p_value,
                    f_stat,
                )
                success_count += 1
                log_msg(
                    f"    [Success] {g1} vs {g2} | "
                    f"n={sample_size} | F={f_stat:.4f} | R2={r2:.4f} | P-val={p_value:.4f}"
                )
            except Exception as error:
                error_count += 1
                log_msg(f"    [Error] Failed on {g1} vs {g2}: {error}")

    del dm
    log_msg(f"=== Finished processing {metric_name} ===")
    log_msg(
        "Processing finished. "
        f"Successful comparisons: {success_count}; "
        f"Skipped comparisons: {skipped_count}; "
        f"Failed comparisons: {error_count}. "
        f"Output saved to: {out_path}"
    )


if __name__ == "__main__":
    main()

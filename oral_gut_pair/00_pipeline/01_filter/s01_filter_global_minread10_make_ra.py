#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
s03_filter_global_minread10_make_ra.py

Integrate previous s3-s6 filtering/RA preparation into one clean step for the
32 non-time-series oral-gut paired BioProjects.

Main design:
  1. Read one merged BIOM table and one metadata table.
  2. Match samples between metadata and BIOM.
  3. Filter ASVs globally across all matched samples:
       keep ASV if total_count_across_all_samples > min_read_count
     Default: min_read_count = 10, operator = ">".
  4. Write one filtered BIOM table, not per-BioProject BIOM files.
  5. Write one matched metadata file corresponding to filtered BIOM samples.
  6. Summarize retained reads and ASVs by BioProject after filtering.
  7. Convert filtered count BIOM to sample-wise relative abundance BIOM.

Default outputs under 01.filter_data:
  filter.biom
  filter_metadata.csv
  filter_ra.biom
  feature_filter_summary.tsv
  filter_summary_by_bioproject.tsv
  sample_depth_summary.tsv
  filter_run_summary.json
  samples_metadata_not_in_biom.txt
  samples_biom_not_in_metadata.txt
  features_keep_total_count_gt10.txt
  features_remove_total_count_not_gt10.txt
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


DEFAULT_METADATA = (
    "/beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp/"
    "00.data/metadata_non_time_series_32bp_in_biom.csv"
)

DEFAULT_BIOM = (
    "/beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp/"
    "00.data/non_time_series_32bp.biom"
)

DEFAULT_OUT_DIR = (
    "/beegfs/zhangcaiyi/oral_gut_pair_noTimeSeries_32bp/"
    "01.filter_data"
)


def read_metadata(path):
    path = str(path)
    if path.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def load_biom_table(path):
    try:
        from biom import load_table
    except ImportError as e:
        raise ImportError(
            "Cannot import biom-format. Install it first:\n"
            "  conda install -c bioconda biom-format\n"
            "or:\n"
            "  pip install biom-format"
        ) from e
    return load_table(str(path))


def write_biom_table(table, out_path, generated_by):
    try:
        from biom.util import biom_open
    except ImportError as e:
        raise ImportError("Cannot import biom.util.biom_open.") from e

    with biom_open(str(out_path), "w") as f:
        table.to_hdf5(f, generated_by=generated_by)


def write_list(items, path):
    with open(path, "w") as f:
        for x in items:
            f.write(str(x) + "\n")


def resolve_path(x, must_exist=False):
    p = Path(x).expanduser()
    if must_exist:
        return p.resolve(strict=True)
    return p.resolve(strict=False)


def check_required_columns(md, columns):
    missing = [c for c in columns if c not in md.columns]
    if missing:
        raise ValueError(
            "Missing required metadata columns:\n  "
            + "\n  ".join(missing)
            + "\n\nAvailable columns:\n  "
            + "\n  ".join(md.columns)
        )


def make_relative_abundance_table(table):
    """
    Convert count BIOM table to sample-wise relative abundance.

    BIOM matrix shape is observations/features x samples.
    For each sample j: RA_ij = count_ij / sample_total_j.
    Samples with zero total count are retained as all-zero columns.
    """
    from biom import Table

    mat = table.matrix_data.tocsc(copy=True).astype(np.float64)
    sample_sums = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)

    nonzero = sample_sums > 0
    inv_sums = np.zeros_like(sample_sums, dtype=np.float64)
    inv_sums[nonzero] = 1.0 / sample_sums[nonzero]

    mat = mat @ sparse.diags(inv_sums, offsets=0, format="csc")
    mat.eliminate_zeros()

    obs_ids = list(table.ids(axis="observation"))
    sample_ids = list(table.ids(axis="sample"))

    obs_md_iter = table.metadata(axis="observation")
    sample_md_iter = table.metadata(axis="sample")

    obs_md = None if obs_md_iter is None else list(obs_md_iter)
    sample_md = None if sample_md_iter is None else list(sample_md_iter)

    ra_table = Table(
        mat,
        observation_ids=obs_ids,
        sample_ids=sample_ids,
        observation_metadata=obs_md,
        sample_metadata=sample_md,
        table_id=(
            str(table.table_id) + "_relative_abundance"
            if table.table_id is not None
            else "relative_abundance"
        ),
    )

    return ra_table, sample_sums


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter merged 32bp BIOM globally by ASV total count across all matched samples, "
            "write filtered BIOM/metadata, summarize by BioProject, and create RA BIOM."
        )
    )

    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--biom", default=DEFAULT_BIOM)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)

    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--bioproject-col", default="bioproject")
    parser.add_argument("--site-col", default="site_group_cleaned")

    parser.add_argument(
        "--min-read-count",
        type=int,
        default=10,
        help=(
            "Global ASV total count threshold across all matched samples. "
            "Default: 10. With default operator 'gt', this keeps total_count > 10."
        ),
    )

    parser.add_argument(
        "--operator",
        choices=["gt", "ge"],
        default="gt",
        help="Use 'gt' for total_count > threshold, or 'ge' for total_count >= threshold. Default: gt.",
    )

    parser.add_argument(
        "--output-prefix",
        default="filter",
        help="Output prefix. Default creates filter.biom, filter_metadata.csv, and filter_ra.biom.",
    )

    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files.")

    args = parser.parse_args()

    metadata_fp = resolve_path(args.metadata, must_exist=True)
    biom_fp = resolve_path(args.biom, must_exist=True)
    out_dir = resolve_path(args.out_dir, must_exist=False)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_biom = out_dir / f"{args.output_prefix}.biom"
    out_metadata = out_dir / f"{args.output_prefix}_metadata.csv"
    out_ra_biom = out_dir / f"{args.output_prefix}_ra.biom"

    feature_summary_tsv = out_dir / "feature_filter_summary.tsv"
    bp_summary_tsv = out_dir / "filter_summary_by_bioproject.tsv"
    sample_depth_tsv = out_dir / "sample_depth_summary.tsv"
    run_summary_json = out_dir / "filter_run_summary.json"

    keep_features_fp = out_dir / f"features_keep_total_count_{args.operator}{args.min_read_count}.txt"
    remove_features_fp = out_dir / f"features_remove_total_count_not_{args.operator}{args.min_read_count}.txt"

    samples_metadata_not_in_biom_fp = out_dir / "samples_metadata_not_in_biom.txt"
    samples_biom_not_in_metadata_fp = out_dir / "samples_biom_not_in_metadata.txt"

    output_files = [
        out_biom,
        out_metadata,
        out_ra_biom,
        feature_summary_tsv,
        bp_summary_tsv,
        sample_depth_tsv,
        run_summary_json,
        keep_features_fp,
        remove_features_fp,
        samples_metadata_not_in_biom_fp,
        samples_biom_not_in_metadata_fp,
    ]

    if not args.overwrite:
        existing = [p for p in output_files if p.exists()]
        if existing:
            raise FileExistsError(
                "Output file exists. Use --overwrite to replace existing files. "
                f"First existing file: {existing[0]}"
            )

    print("Resolved paths:", flush=True)
    print(f"  metadata: {metadata_fp}", flush=True)
    print(f"  biom:     {biom_fp}", flush=True)
    print(f"  out_dir:  {out_dir}", flush=True)

    print("[1/7] Reading metadata ...", flush=True)
    md = read_metadata(metadata_fp)

    check_required_columns(md, [args.sample_col, args.bioproject_col, args.site_col])

    md[args.sample_col] = md[args.sample_col].astype(str).str.strip()
    md[args.bioproject_col] = md[args.bioproject_col].astype(str).str.strip()
    md[args.site_col] = md[args.site_col].astype(str).str.strip()

    duplicated = (
        md.loc[md[args.sample_col].duplicated(keep=False), args.sample_col]
        .sort_values()
        .unique()
        .tolist()
    )

    if duplicated:
        dup_path = out_dir / "duplicated_sample_ids_in_metadata.txt"
        write_list(duplicated, dup_path)
        raise RuntimeError(f"Duplicated sample IDs in metadata: {len(duplicated)}. See {dup_path}")

    print("[2/7] Loading BIOM ...", flush=True)
    table = load_biom_table(biom_fp)

    biom_samples_all = [str(x).strip() for x in table.ids(axis="sample")]
    biom_features_all = [str(x).strip() for x in table.ids(axis="observation")]

    biom_sample_set = set(biom_samples_all)
    md_sample_set = set(md[args.sample_col])

    samples_in_both = [s for s in biom_samples_all if s in md_sample_set]
    samples_metadata_not_in_biom = sorted(md_sample_set - biom_sample_set)
    samples_biom_not_in_metadata = sorted(biom_sample_set - md_sample_set)

    write_list(samples_metadata_not_in_biom, samples_metadata_not_in_biom_fp)
    write_list(samples_biom_not_in_metadata, samples_biom_not_in_metadata_fp)

    if not samples_in_both:
        raise RuntimeError("No matched samples between metadata and BIOM.")

    print(f"BIOM samples before matching: {len(biom_samples_all)}", flush=True)
    print(f"Metadata samples: {md[args.sample_col].nunique()}", flush=True)
    print(f"Matched samples used: {len(samples_in_both)}", flush=True)
    print(f"Features before filtering: {len(biom_features_all)}", flush=True)

    print("[3/7] Matching BIOM samples and metadata ...", flush=True)

    table_matched = table.filter(samples_in_both, axis="sample", inplace=False)

    sample_ids = [str(x).strip() for x in table_matched.ids(axis="sample")]
    feature_ids = [str(x).strip() for x in table_matched.ids(axis="observation")]

    md_indexed = md.set_index(args.sample_col, drop=False)
    md_matched = md_indexed.loc[sample_ids].reset_index(drop=True)

    print("[4/7] Global ASV filtering across all matched samples ...", flush=True)

    mat = table_matched.matrix_data.tocsr()
    feature_total_counts = np.asarray(mat.sum(axis=1)).ravel().astype(np.float64)
    sample_total_counts_before = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)

    if args.operator == "gt":
        keep_mask = feature_total_counts > args.min_read_count
        criterion_text = f"total_count_across_all_samples > {args.min_read_count}"
    else:
        keep_mask = feature_total_counts >= args.min_read_count
        criterion_text = f"total_count_across_all_samples >= {args.min_read_count}"

    keep_features = [fid for fid, keep in zip(feature_ids, keep_mask) if keep]
    remove_features = [fid for fid, keep in zip(feature_ids, keep_mask) if not keep]

    print(f"Filtering criterion: {criterion_text}", flush=True)
    print(f"Features retained: {len(keep_features)} / {len(feature_ids)}", flush=True)

    if len(keep_features) == 0:
        raise RuntimeError("No ASVs retained after filtering. Check threshold/operator.")

    table_filtered = table_matched.filter(keep_features, axis="observation", inplace=False)

    filtered_mat = table_filtered.matrix_data.tocsr()
    sample_total_counts_after = np.asarray(filtered_mat.sum(axis=0)).ravel().astype(np.float64)
    retained_total_reads = float(sample_total_counts_after.sum())

    write_list(keep_features, keep_features_fp)
    write_list(remove_features, remove_features_fp)

    feature_summary = pd.DataFrame({
        "feature_id": feature_ids,
        "total_count_across_all_samples": feature_total_counts.astype(np.int64),
        "keep": keep_mask.astype(bool),
        "filter_criterion": criterion_text,
    }).sort_values(
        ["keep", "total_count_across_all_samples", "feature_id"],
        ascending=[False, False, True],
    )

    feature_summary.to_csv(feature_summary_tsv, sep="\t", index=False)

    print("[5/7] Writing filtered BIOM and matched metadata ...", flush=True)

    write_biom_table(
        table_filtered,
        out_biom,
        generated_by=(
            "s03_filter_global_minread10_make_ra.py; "
            f"{criterion_text}; matched samples only"
        ),
    )

    md_matched.to_csv(out_metadata, index=False)

    print("[6/7] Summarizing retained reads/features by BioProject ...", flush=True)

    bp_arr = md_matched[args.bioproject_col].astype(str).values
    site_arr = md_matched[args.site_col].astype(str).values

    bp_rows = []
    sample_depth_rows = []

    for sample_id, bp, site, before, after in zip(
        sample_ids,
        bp_arr,
        site_arr,
        sample_total_counts_before,
        sample_total_counts_after,
    ):
        sample_depth_rows.append({
            "sample_id": sample_id,
            "bioproject": bp,
            "site": site,
            "sample_reads_before_filter": int(before),
            "sample_reads_after_filter": int(after),
            "fraction_reads_retained": (float(after) / float(before)) if before > 0 else np.nan,
        })

    sample_depth_df = pd.DataFrame(sample_depth_rows)
    sample_depth_df.to_csv(sample_depth_tsv, sep="\t", index=False)

    for bp in sorted(set(bp_arr)):
        bp_sample_mask = bp_arr == bp
        sample_indices = np.where(bp_sample_mask)[0]

        bp_mat_before = mat[:, sample_indices]
        bp_feature_totals_before = np.asarray(bp_mat_before.sum(axis=1)).ravel()

        bp_mat_after = filtered_mat[:, sample_indices]
        bp_feature_totals_after = np.asarray(bp_mat_after.sum(axis=1)).ravel()

        total_reads_before_bp = float(bp_feature_totals_before.sum())
        total_reads_after_bp = float(bp_feature_totals_after.sum())

        n_asv_before_nonzero_bp = int((bp_feature_totals_before > 0).sum())
        n_asv_after_nonzero_bp = int((bp_feature_totals_after > 0).sum())

        site_counts = pd.Series(site_arr[bp_sample_mask]).value_counts().to_dict()

        bp_rows.append({
            "bioproject": bp,
            "n_samples": int(bp_sample_mask.sum()),
            "n_oral_samples": int(site_counts.get("Oral", 0)),
            "n_gut_samples": int(site_counts.get("Gut", 0)),
            "n_asv_before_filter_nonzero_in_bp": n_asv_before_nonzero_bp,
            "n_asv_after_filter_nonzero_in_bp": n_asv_after_nonzero_bp,
            "n_global_asv_before_filter": int(len(feature_ids)),
            "n_global_asv_after_filter": int(len(keep_features)),
            "total_reads_before_filter_bp": int(total_reads_before_bp),
            "total_reads_after_filter_bp": int(total_reads_after_bp),
            "reads_removed_by_filter_bp": int(total_reads_before_bp - total_reads_after_bp),
            "fraction_reads_retained_bp": (
                total_reads_after_bp / total_reads_before_bp if total_reads_before_bp > 0 else np.nan
            ),
            "filter_criterion": criterion_text,
        })

    bp_summary = pd.DataFrame(bp_rows).sort_values("bioproject")
    bp_summary.to_csv(bp_summary_tsv, sep="\t", index=False)

    print("[7/7] Converting filtered BIOM to relative abundance ...", flush=True)

    ra_table, filtered_sample_sums = make_relative_abundance_table(table_filtered)
    write_biom_table(
        ra_table,
        out_ra_biom,
        generated_by=(
            "s03_filter_global_minread10_make_ra.py; "
            "relative_abundance = feature_count / sample_total_count within each sample"
        ),
    )

    ra_sample_sums = np.asarray(ra_table.matrix_data.sum(axis=0)).ravel()
    nonzero_filtered_samples = filtered_sample_sums > 0

    if nonzero_filtered_samples.any():
        max_abs_error_ra_sum = float(np.max(np.abs(ra_sample_sums[nonzero_filtered_samples] - 1.0)))
        min_ra_sum_nonzero = float(np.min(ra_sample_sums[nonzero_filtered_samples]))
        max_ra_sum_nonzero = float(np.max(ra_sample_sums[nonzero_filtered_samples]))
    else:
        max_abs_error_ra_sum = None
        min_ra_sum_nonzero = None
        max_ra_sum_nonzero = None

    run_summary = {
        "input_metadata": str(metadata_fp),
        "input_biom": str(biom_fp),
        "out_dir": str(out_dir),
        "outputs": {
            "filtered_biom": str(out_biom),
            "filtered_metadata": str(out_metadata),
            "filtered_ra_biom": str(out_ra_biom),
            "feature_filter_summary": str(feature_summary_tsv),
            "filter_summary_by_bioproject": str(bp_summary_tsv),
            "sample_depth_summary": str(sample_depth_tsv),
            "features_keep": str(keep_features_fp),
            "features_remove": str(remove_features_fp),
            "samples_metadata_not_in_biom": str(samples_metadata_not_in_biom_fp),
            "samples_biom_not_in_metadata": str(samples_biom_not_in_metadata_fp),
        },
        "sample_matching": {
            "n_biom_samples_before_matching": int(len(biom_samples_all)),
            "n_metadata_samples": int(md[args.sample_col].nunique()),
            "n_samples_in_both": int(len(samples_in_both)),
            "n_metadata_samples_not_in_biom": int(len(samples_metadata_not_in_biom)),
            "n_biom_samples_not_in_metadata": int(len(samples_biom_not_in_metadata)),
        },
        "filter_rule": {
            "level": "global across all matched samples in the merged BIOM",
            "criterion": criterion_text,
            "min_read_count": int(args.min_read_count),
            "operator": args.operator,
            "note": (
                "With operator='gt' and min_read_count=10, ASVs with total count "
                "across all matched samples equal to 10 are removed; only ASVs with "
                "total count > 10 are kept."
            ),
        },
        "feature_counts": {
            "n_features_before_filter": int(len(feature_ids)),
            "n_features_after_filter": int(len(keep_features)),
            "n_features_removed": int(len(remove_features)),
            "fraction_features_retained": float(len(keep_features) / len(feature_ids)),
        },
        "read_counts": {
            "total_reads_before_filter": int(sample_total_counts_before.sum()),
            "total_reads_after_filter": int(retained_total_reads),
            "reads_removed_by_filter": int(sample_total_counts_before.sum() - retained_total_reads),
            "fraction_reads_retained": (
                float(retained_total_reads / sample_total_counts_before.sum())
                if sample_total_counts_before.sum() > 0 else None
            ),
        },
        "relative_abundance": {
            "definition": "feature_count / sample_total_count within each sample after filtering",
            "n_zero_depth_samples_after_filter": int((filtered_sample_sums == 0).sum()),
            "min_ra_sum_nonzero_samples": min_ra_sum_nonzero,
            "max_ra_sum_nonzero_samples": max_ra_sum_nonzero,
            "max_abs_error_ra_sum_from_1_nonzero_samples": max_abs_error_ra_sum,
        },
        "n_bioprojects": int(bp_summary["bioproject"].nunique()),
        "bioproject_column": args.bioproject_col,
        "sample_column": args.sample_col,
        "site_column": args.site_col,
    }

    with open(run_summary_json, "w") as f:
        json.dump(run_summary, f, indent=2)

    print("\nDone.", flush=True)
    print("Main outputs:", flush=True)
    print(f"  filtered BIOM:      {out_biom}", flush=True)
    print(f"  filtered metadata:  {out_metadata}", flush=True)
    print(f"  filtered RA BIOM:   {out_ra_biom}", flush=True)
    print(f"  BP summary:         {bp_summary_tsv}", flush=True)
    print(f"  feature summary:    {feature_summary_tsv}", flush=True)
    print(f"  run summary JSON:   {run_summary_json}", flush=True)
    print("", flush=True)
    print(f"Features retained: {len(keep_features)} / {len(feature_ids)}", flush=True)
    print(f"Reads retained: {int(retained_total_reads)} / {int(sample_total_counts_before.sum())}", flush=True)


if __name__ == "__main__":
    main()

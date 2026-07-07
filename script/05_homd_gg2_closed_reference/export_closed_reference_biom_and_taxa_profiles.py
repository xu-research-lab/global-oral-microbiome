#!/usr/bin/env python3
"""Export containment-assisted GG2 closed-reference profiles to BIOM.

Outputs:
- GG2 reference-ID level BIOM from the sparse closed-reference profile.
- Phylum/class/family/genus collapsed BIOM tables, mapped-only.
- Phylum/class/family/genus collapsed BIOM tables with an additional Unmapped
  feature per sample to preserve total read counts for QC/sensitivity analyses.
- Sparse long-format TSV and feature summaries for each exported table.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
from biom import Table
from scipy import sparse


DEFAULT_PROFILE_DIR = (
    "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/"
    "containment_derep/unified_gg2_closed_reference_profile"
)
DEFAULT_TAXONOMY = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"


LEVELS = {
    "phylum": "p__",
    "class": "c__",
    "family": "f__",
    "genus": "g__",
}


def log(message: str) -> None:
    print(f"[export] {message}", flush=True)


def parse_taxon(taxon: str, prefix: str, level: str) -> str:
    if not isinstance(taxon, str) or not taxon.strip():
        return f"Unclassified_{level}"
    for part in taxon.split(";"):
        value = part.strip()
        if value.startswith(prefix):
            if value == prefix or value in {f"{prefix}s__", "s__"}:
                return f"Unclassified_{level}"
            return value
    return f"Unclassified_{level}"


def load_taxonomy(path: Path) -> pd.DataFrame:
    tax = pd.read_csv(path, sep="\t", dtype=str)
    if "Feature ID" not in tax.columns or "Taxon" not in tax.columns:
        raise ValueError(f"taxonomy file must contain 'Feature ID' and 'Taxon': {path}")
    tax = tax.rename(columns={"Feature ID": "gg2_target", "Taxon": "taxonomy"})
    tax["taxonomy"] = tax["taxonomy"].fillna("")
    for level, prefix in LEVELS.items():
        tax[level] = tax["taxonomy"].map(lambda x: parse_taxon(x, prefix, level))
    return tax[["gg2_target", "taxonomy", *LEVELS.keys()]]


def write_biom_from_long(df: pd.DataFrame, feature_col: str, out_biom: Path, generated_by: str) -> None:
    features = pd.Index(sorted(df[feature_col].astype(str).unique()))
    samples = pd.Index(sorted(df["sample_key"].astype(str).unique()))
    feature_codes = features.get_indexer(df[feature_col].astype(str))
    sample_codes = samples.get_indexer(df["sample_key"].astype(str))
    data = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(float).to_numpy()
    mat = sparse.coo_matrix((data, (feature_codes, sample_codes)), shape=(len(features), len(samples))).tocsr()
    table = Table(mat, observation_ids=features.tolist(), sample_ids=samples.tolist(), type="OTU table")
    with h5py.File(out_biom, "w") as handle:
        table.to_hdf5(handle, generated_by)


def write_outputs(df: pd.DataFrame, feature_col: str, prefix: Path, generated_by: str) -> None:
    out_long = prefix.with_suffix(".long.tsv")
    out_summary = prefix.with_suffix(".feature_summary.tsv")
    out_biom = prefix.with_suffix(".biom")

    df = df[["sample_key", feature_col, "count"]].rename(columns={feature_col: "feature_id"})
    df = df.groupby(["sample_key", "feature_id"], as_index=False, sort=False)["count"].sum()
    df.to_csv(out_long, sep="\t", index=False)

    summary = (
        df.groupby("feature_id", as_index=False)
        .agg(total_reads=("count", "sum"), sample_prevalence=("sample_key", "nunique"))
        .sort_values("total_reads", ascending=False)
    )
    summary.to_csv(out_summary, sep="\t", index=False)
    write_biom_from_long(df.rename(columns={"feature_id": feature_col}), feature_col, out_biom, generated_by)
    log(f"wrote {out_biom}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-dir", default=DEFAULT_PROFILE_DIR)
    parser.add_argument("--taxonomy", default=DEFAULT_TAXONOMY)
    args = parser.parse_args()

    profile_dir = Path(args.profile_dir)
    outdir = profile_dir / "biom_exports"
    outdir.mkdir(parents=True, exist_ok=True)
    generated_by = "containment-assisted GG2 closed-reference export"

    profile_path = profile_dir / "all_samples_gg2_closed_reference_profile.long.tsv"
    sample_fraction_path = profile_dir / "all_samples_gg2_closed_reference_mapping_fraction.tsv"

    log("loading sparse closed-reference profile")
    df = pd.read_csv(profile_path, sep="\t", usecols=["sample_key", "gg2_target", "count"], dtype={"sample_key": str, "gg2_target": str})
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0)
    log(f"profile rows={len(df):,}; samples={df.sample_key.nunique():,}; gg2 targets={df.gg2_target.nunique():,}")

    log("exporting GG2 reference-ID BIOM")
    write_outputs(df.rename(columns={"gg2_target": "feature_id"}), "feature_id", outdir / "gg2_reference_id_profile", generated_by)

    log("loading GG2 taxonomy")
    tax = load_taxonomy(Path(args.taxonomy))
    df_tax = df.merge(tax, on="gg2_target", how="left")
    for level in LEVELS:
        df_tax[level] = df_tax[level].fillna(f"Unclassified_{level}")

    log("loading per-sample unmapped reads")
    sample_fraction = pd.read_csv(sample_fraction_path, sep="\t", usecols=["sample_key", "unmapped_reads"], dtype={"sample_key": str})
    sample_fraction["unmapped_reads"] = pd.to_numeric(sample_fraction["unmapped_reads"], errors="coerce").fillna(0)
    unmapped = sample_fraction[sample_fraction["unmapped_reads"] > 0].copy()

    for level in ["phylum", "class", "family", "genus"]:
        log(f"collapsing mapped-only {level}")
        collapsed = (
            df_tax.groupby(["sample_key", level], as_index=False, sort=False)["count"]
            .sum()
        )
        write_outputs(collapsed, level, outdir / f"{level}_profile_mapped_only", generated_by)

        log(f"collapsing with-unmapped {level}")
        unmapped_level = unmapped.rename(columns={"unmapped_reads": "count"}).copy()
        unmapped_level[level] = "Unmapped"
        with_unmapped = pd.concat([collapsed[["sample_key", level, "count"]], unmapped_level[["sample_key", level, "count"]]], ignore_index=True)
        write_outputs(with_unmapped, level, outdir / f"{level}_profile_with_unmapped", generated_by)

    manifest_rows = []
    for path in sorted(outdir.glob("*")):
        if path.is_file():
            manifest_rows.append({"file": str(path), "size_bytes": path.stat().st_size})
    pd.DataFrame(manifest_rows).to_csv(outdir / "manifest.tsv", sep="\t", index=False)
    log(f"complete: {outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
da_v4 data cleaning: v2-qc BIOM -> genus counts + TSS relative abundance.

Shared steps for every depth track:
  1. drop samples with total reads < MIN_READS (1000)
  2. optional rarefaction to --depth (without replacement, seed 42);
     samples below the depth are dropped
  3. aggregate features to genus via taxonomy.tsv (g__ level)
  4. keep a genus if its prevalence is >= MIN_PREVALENCE (1%) in ANY eligible
     disease matched cohort (no site cohorts)

Tracks:
  raw     : python3 prepare_genus_tables.py
  rare1k  : python3 prepare_genus_tables.py --depth 1000 --tag _rare1k
  rare5k  : python3 prepare_genus_tables.py --depth 5000 --tag _rare5k
"""
from __future__ import print_function
import os
from pathlib import Path
import argparse
import csv
import numpy as np
import pandas as pd
import biom
from scipy import sparse

# Repository root defaults to this script's folder; raw inputs (not
# distributed, see README) default to data/ under it.
WORK_DIR = Path(os.environ.get("DA_WORK_DIR", Path(__file__).resolve().parent.parent))
DATA_DIR = Path(os.environ.get("DA_DATA_DIR", WORK_DIR / "data"))
BIOM_PATH = DATA_DIR / "biom_oral_qc_dedup.biom"
META_PATH = DATA_DIR / "metadata_oral_qc_dedup.tsv"
TAXONOMY_PATH = DATA_DIR / "taxonomy.tsv"
ELIGIBLE_PATH = WORK_DIR / "eligible_diseases.txt"

MIN_READS = 1000
RARE_SEED = 42
MIN_PREVALENCE = 0.01

SITE_MAP = {
    "saliva": "Saliva", "tongue": "Tongue", "tooth": "Tooth",
    "oral swab": "Oral_swab", "oral mucosa": "Oral_Mucosa",
    "gingiva": "Gingiva", "other": "Others",
}
STATUS_MAP = {"HC": "Healthy Control", "Healthy": "Healthy Control"}
HC = "Healthy Control"


def parse_args():
    ap = argparse.ArgumentParser(description="da_v4 genus table preparation")
    ap.add_argument(
        "--depth", type=int, default=None,
        help="rarefaction depth; omit for the non-rarefied (raw) track"
    )
    ap.add_argument(
        "--tag", default="",
        help="output suffix, e.g. _rare1k / _rare5k (default: empty = raw)"
    )
    return ap.parse_args()


def extract_genus(taxon, fallback):
    for part in str(taxon).split(";"):
        part = part.strip()
        if part.startswith("g__"):
            genus = part[3:].strip()
            return genus if genus else fallback
    return fallback


def prevalence(binary_mat, sample_idx):
    if len(sample_idx) == 0:
        return None
    return np.asarray(binary_mat[:, sample_idx].getnnz(axis=1)).ravel() / float(
        len(sample_idx)
    )


def load_metadata():
    meta = pd.read_csv(META_PATH, sep="\t", low_memory=False)
    str_cols = [
        "sample_id", "bioproject_id", "site_cleaned", "continent_group",
        "status_cleaned", "sex_cleaned", "age_cleaned_group",
        "sequencing_region_cleaned",
    ]
    for col in str_cols:
        if col in meta.columns:
            meta[col] = meta[col].astype(str).str.strip()
            meta.loc[meta[col].isin(["", "nan", "None", "NA"]), col] = pd.NA
    if "site_cleaned" in meta.columns:
        meta["site_cleaned"] = meta["site_cleaned"].replace(SITE_MAP)
    if "status_cleaned" in meta.columns:
        meta["status_cleaned"] = meta["status_cleaned"].replace(STATUS_MAP)
    if "age_cleaned_group" in meta.columns:
        meta["age_cleaned_group"] = meta["age_cleaned_group"].str.replace(
            "＞", ">", regex=False
        )
    return meta.dropna(subset=["sample_id"])


def load_eligible_diseases():
    if not ELIGIBLE_PATH.exists():
        raise SystemExit(
            "Missing %s — run compute_eligible_diseases.py first" % ELIGIBLE_PATH
        )
    diseases = [
        line.strip()
        for line in ELIGIBLE_PATH.read_text().splitlines()
        if line.strip()
    ]
    print("Eligible diseases (%d): %s" % (len(diseases), ", ".join(diseases)))
    return diseases


def rarefy_table(table, depth, seed):
    """Rarefy BIOM table; samples with library size < depth are dropped."""
    try:
        return table.subsample(
            depth, axis="sample", by_id=False,
            with_replacement=False, seed=seed
        )
    except TypeError:
        print("biom.Table.subsample has no seed=; using NumPy rarefaction")
        return rarefy_numpy(table, depth, seed)


def rarefy_numpy(table, depth, seed):
    """Per-sample rarefaction without replacement (NumPy fallback)."""
    rng = np.random.RandomState(seed)
    obs_ids = np.array(table.ids(axis="observation"), dtype=str)
    sample_ids = np.array(table.ids(axis="sample"), dtype=str)
    mat = table.matrix_data.tocsc()
    totals = np.asarray(mat.sum(axis=0)).ravel()
    keep_s = np.where(totals >= depth)[0]
    print("  NumPy rarefy: keeping %d / %d samples" % (len(keep_s), len(sample_ids)))

    data, rows, cols = [], [], []
    for new_j, j in enumerate(keep_s):
        col = mat.getcol(j)
        expanded = np.repeat(col.indices, col.data.astype(np.int64))
        chosen = rng.choice(expanded, size=depth, replace=False)
        uniq, cnt = np.unique(chosen, return_counts=True)
        rows.extend(uniq.tolist())
        cols.extend([new_j] * len(uniq))
        data.extend(cnt.tolist())

    out = sparse.coo_matrix(
        (data, (rows, cols)), shape=(len(obs_ids), len(keep_s)), dtype=np.float64
    ).tocsr()
    return biom.Table(
        out, observation_ids=obs_ids.tolist(),
        sample_ids=sample_ids[keep_s].tolist()
    )


def disease_cohort(meta, disease):
    """Matched-BP cohort for one disease; returns None if the filter fails."""
    disease_meta = meta[
        meta["status_cleaned"].notna()
        & ~meta["status_cleaned"].isin(["Other", HC])
        & meta["bioproject_id"].notna()
    ]
    if disease not in set(disease_meta["status_cleaned"]):
        return None
    disease_bp = disease_meta.loc[
        disease_meta["status_cleaned"] == disease, "bioproject_id"
    ].unique()
    sub = meta[
        meta["status_cleaned"].isin([HC, disease])
        & meta["bioproject_id"].isin(disease_bp)
    ]
    if sub.empty:
        return None
    bp_tbl = pd.crosstab(sub["bioproject_id"], sub["status_cleaned"])
    if disease not in bp_tbl.columns or HC not in bp_tbl.columns:
        return None
    valid_bp = bp_tbl.index[(bp_tbl[disease] > 0) & (bp_tbl[HC] > 0)]
    sub = sub[sub["bioproject_id"].isin(valid_bp)]
    n_case = int(sub["status_cleaned"].eq(disease).sum())
    n_hc = int(sub["status_cleaned"].eq(HC).sum())
    n_bp = int(sub["bioproject_id"].nunique())
    passes = n_case >= 30 and n_hc >= 30 and n_bp >= 2
    return sub, n_case, n_hc, n_bp, passes


def main():
    args = parse_args()
    tag = args.tag
    depth = args.depth
    track = "raw" if depth is None else "rare%d" % depth
    WORK_DIR.mkdir(parents=True, exist_ok=True)

    counts_path = WORK_DIR / ("genus_table_counts%s.tsv" % tag)
    relabund_path = WORK_DIR / ("genus_table_relabund%s.tsv" % tag)
    report_path = WORK_DIR / ("data_cleaning_report%s.csv" % tag)
    print("=== da_v4 prepare genus tables: track=%s tag='%s' ===" % (track, tag))

    eligible_diseases = load_eligible_diseases()

    print("Loading taxonomy:", TAXONOMY_PATH)
    tax = pd.read_csv(TAXONOMY_PATH, sep="\t")
    tax = tax.rename(columns={"Feature ID": "FeatureID", "Taxon": "Taxon"})
    tax_map = dict(zip(tax["FeatureID"].astype(str), tax["Taxon"].astype(str)))

    print("Loading BIOM:", BIOM_PATH)
    table = biom.load_table(str(BIOM_PATH))
    n_feat_in, n_samp_in = table.shape
    print("BIOM as loaded: %d features x %d samples" % (n_feat_in, n_samp_in))

    # --- Step 1: shared sample filter, total reads >= MIN_READS ---
    totals = np.asarray(table.sum(axis="sample")).ravel()
    print(
        "Library size: median=%.0f p10=%.0f min=%.0f max=%.0f"
        % (
            float(np.median(totals)),
            float(np.percentile(totals, 10)),
            float(totals.min()),
            float(totals.max()),
        )
    )
    ids_all = np.array(table.ids(axis="sample"), dtype=str)
    keep_ids = ids_all[totals >= MIN_READS]
    n_after_minreads = len(keep_ids)
    print(
        "Step1 reads >= %d: keep %d / %d samples (drop %d)"
        % (MIN_READS, n_after_minreads, n_samp_in, n_samp_in - n_after_minreads)
    )
    table = table.filter(set(keep_ids), axis="sample", inplace=False)

    # --- Step 2: optional rarefaction ---
    if depth is None:
        print("Step2 rarefaction: skipped (raw track)")
    else:
        print("Step2 rarefying to depth=%d seed=%d ..." % (depth, RARE_SEED))
        table = rarefy_table(table, depth, RARE_SEED)
        chk = np.asarray(table.sum(axis="sample")).ravel()
        print(
            "  after rarefy: %d samples, depth min=%.0f max=%.0f"
            % (table.shape[1], chk.min(), chk.max())
        )
    n_samp_out = table.shape[1]
    print("Samples entering genus aggregation: %d" % n_samp_out)

    # --- Step 3: genus aggregation ---
    obs_ids = np.array(table.ids(axis="observation"), dtype=str)
    sample_ids = np.array(table.ids(axis="sample"), dtype=str)
    mat = table.matrix_data.tocsr()
    genera = np.array(
        [extract_genus(tax_map.get(obs, ""), obs) for obs in obs_ids], dtype=object
    )
    unique_genera, inverse = np.unique(genera, return_inverse=True)
    group = sparse.csr_matrix(
        (np.ones(len(obs_ids), dtype=np.float64), (inverse, np.arange(len(obs_ids)))),
        shape=(len(unique_genera), len(obs_ids)),
    )
    genus_counts = (group @ mat).tocsr()
    binary = genus_counts.copy()
    binary.data = np.ones_like(binary.data)
    print("Step3 unique genera before prevalence filter: %d" % len(unique_genera))

    # --- Step 4: prevalence >= MIN_PREVALENCE in any eligible disease cohort ---
    meta = load_metadata()
    meta["sample_id"] = meta["sample_id"].astype(str)
    sample_to_idx = {sid: i for i, sid in enumerate(sample_ids)}
    meta = meta[meta["sample_id"].isin(sample_to_idx)].copy()
    meta["_idx"] = meta["sample_id"].map(sample_to_idx).astype(int)
    print("Samples with metadata: %d" % len(meta))

    keep = np.zeros(len(unique_genera), dtype=bool)
    report_rows = [{
        "scope": "overall",
        "track": track,
        "tag": tag,
        "Disease": "",
        "min_reads": MIN_READS,
        "rare_depth": "" if depth is None else depth,
        "rare_seed": "" if depth is None else RARE_SEED,
        "min_prevalence": MIN_PREVALENCE,
        "n_samples_biom": n_samp_in,
        "n_samples_after_min_reads": n_after_minreads,
        "n_samples_final": n_samp_out,
        "n_case": "",
        "n_hc": "",
        "n_bp": "",
        "passes_filter": "",
        "n_genera_before_filter": len(unique_genera),
        "n_genera_kept": "",
    }]

    for disease in eligible_diseases:
        res = disease_cohort(meta, disease)
        if res is None:
            print("  Disease '%s': not present in this track" % disease)
            continue
        sub, n_case, n_hc, n_bp, passes = res
        if not passes:
            print(
                "  Disease '%s': SKIP for prevalence (case=%d HC=%d BP=%d)"
                % (disease, n_case, n_hc, n_bp)
            )
        else:
            prev = prevalence(binary, sub["_idx"].to_numpy())
            keep |= prev >= MIN_PREVALENCE
            print(
                "  Disease '%s': n=%d (case=%d HC=%d BP=%d), genera >=%.0f%%: %d"
                % (
                    disease, len(sub), n_case, n_hc, n_bp,
                    MIN_PREVALENCE * 100, int((prev >= MIN_PREVALENCE).sum()),
                )
            )
        report_rows.append({
            "scope": "disease",
            "track": track,
            "tag": tag,
            "Disease": disease,
            "min_reads": MIN_READS,
            "rare_depth": "" if depth is None else depth,
            "rare_seed": "" if depth is None else RARE_SEED,
            "min_prevalence": MIN_PREVALENCE,
            "n_samples_biom": "",
            "n_samples_after_min_reads": "",
            "n_samples_final": len(sub),
            "n_case": n_case,
            "n_hc": n_hc,
            "n_bp": n_bp,
            "passes_filter": passes,
            "n_genera_before_filter": "",
            "n_genera_kept": "",
        })

    kept = np.where(keep)[0]
    print("Step4 keeping %d / %d genera" % (len(kept), len(unique_genera)))
    genus_counts = genus_counts[kept, :]
    kept_genera = unique_genera[kept]
    report_rows[0]["n_genera_kept"] = len(kept_genera)

    nnz = genus_counts.getnnz()
    sparsity = 1.0 - nnz / float(max(1, genus_counts.shape[0] * genus_counts.shape[1]))
    report_rows[0]["sparsity"] = "%.4f" % sparsity
    print("Sparsity (fraction of zeros): %.4f" % sparsity)

    print("Writing cleaning report:", report_path)
    fieldnames = list(report_rows[0].keys())
    with report_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report_rows)

    genus_dense = np.asarray(genus_counts.toarray(), dtype=np.int64)
    print("Writing counts:", counts_path)
    with counts_path.open("w") as fh:
        fh.write("Genus\t" + "\t".join(sample_ids) + "\n")
        for genus, row in zip(kept_genera, genus_dense):
            fh.write(genus + "\t" + "\t".join(str(int(x)) for x in row) + "\n")

    sample_totals = np.asarray(genus_counts.sum(axis=0)).ravel().astype(np.float64)
    sample_totals[sample_totals == 0] = 1.0
    genus_rel = (genus_counts @ sparse.diags(1.0 / sample_totals)).toarray()
    print("Writing relabund:", relabund_path)
    with relabund_path.open("w") as fh:
        fh.write("Genus\t" + "\t".join(sample_ids) + "\n")
        for genus, row in zip(kept_genera, genus_rel):
            fh.write(genus + "\t" + "\t".join("%.8g" % x for x in row) + "\n")

    print(
        "Done track=%s: %d genera x %d samples"
        % (track, len(kept_genera), len(sample_ids))
    )


if __name__ == "__main__":
    main()

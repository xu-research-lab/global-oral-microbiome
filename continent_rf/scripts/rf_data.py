"""Metadata alignment and sparse feature preparation for continent RF."""

import numpy as np
import pandas as pd
from scipy import sparse

REGIONS = ["North America", "Asia", "Europe", "Oceania", "South America", "Africa"]
ORAL_SITES = ["Saliva", "Tongue", "Tooth"]


def parse_genus(taxon):
    if isinstance(taxon, str):
        for part in taxon.split(";"):
            if part.strip().startswith("g__"):
                return part.strip()[3:].strip() or None
    return None


def load_genus_map(path):
    tax = pd.read_csv(path, sep="\t", usecols=["Feature ID", "Taxon"])
    tax["genus"] = tax["Taxon"].map(parse_genus)
    tax = tax.dropna(subset=["genus"])
    return dict(zip(tax["Feature ID"].astype(str), tax["genus"]))


def aggregate_genus(table, genus_map):
    """Sum named genera, discard unassigned features, then normalize samples."""
    labels = [genus_map.get(str(fid)) for fid in table.ids(axis="observation")]
    known = np.array([label is not None for label in labels])
    genera = sorted({label for label in labels if label is not None})
    if not genera:
        raise ValueError("No BIOM features have a named genus in the taxonomy")
    lookup = {genus: i for i, genus in enumerate(genera)}
    columns = [lookup[label] for label in labels if label is not None]
    mapping = sparse.csr_matrix(
        (np.ones(len(columns)), (np.arange(len(columns)), columns)),
        shape=(len(columns), len(genera)),
    )
    matrix = table.matrix_data.T.tocsr()[:, known] @ mapping
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    # Preserve zero rows (samples whose features are all unassigned), as before.
    matrix = sparse.diags(1.0 / np.where(totals > 0, totals, 1)) @ matrix
    print(f"Genus aggregation: {known.sum():,}/{len(labels):,} features -> {len(genera):,} genera")
    return matrix.tocsr()


def load_dataset(dataset, metadata_path, biom_path, *, split_mode="project",
                 min_study_samples=20, genus_map=None):
    """Read BIOM once and retain its sample order for metadata and matrix."""
    import biom

    if dataset not in {"oral", "gut"} or split_mode not in {"project", "unit"}:
        raise ValueError("Invalid dataset or split mode")
    meta = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    rename = ({"continent_group": "region", "site_cleaned": "site",
               "bioproject_id": "project"} if dataset == "oral" else {})
    required = {"sample_id", *(rename if rename else ["region", "project"])}
    missing = sorted(required.difference(meta.columns))
    if missing:
        raise ValueError(f"{dataset} metadata missing columns: {', '.join(missing)}")
    # Select the source columns before renaming to avoid duplicate canonical names.
    meta = meta[list(sorted(required))].rename(columns=rename)
    if dataset == "gut":
        meta["site"] = "Gut"
    text_columns = ["sample_id", "region", "project"] + (["site"] if dataset == "oral" else [])
    for column in text_columns:
        meta[column] = meta[column].astype("string").str.strip()
    meta["site"] = meta["site"].replace({s.lower(): s for s in ORAL_SITES})
    meta = meta[
        meta["region"].isin(REGIONS)
        & meta["site"].isin(ORAL_SITES if dataset == "oral" else ["Gut"])
        & meta["project"].notna() & meta["sample_id"].notna()
    ].copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    if meta["sample_id"].duplicated().any():
        raise ValueError(f"{dataset} metadata contains duplicate sample_id values")
    meta["unit"] = (meta["region"] + "||" + meta["project"]
                    if split_mode == "unit" else meta["project"])
    keys = ["site", "region", "project"]
    sizes = meta.groupby(keys).size().reset_index(name="n")
    valid = sizes.loc[sizes["n"] >= min_study_samples, keys]
    meta = meta.merge(valid, on=keys, how="inner")
    print(f"{dataset}: keep {len(meta):,} samples from {len(valid):,} site-region-project groups")

    table = biom.load_table(str(biom_path))
    wanted = set(meta["sample_id"])
    shared = [sid for sid in table.ids(axis="sample") if str(sid) in wanted]
    if not shared:
        raise ValueError(f"No {dataset} sample IDs overlap between metadata and BIOM")
    table = table.filter(shared, axis="sample", inplace=False)
    table.remove_empty(inplace=True)
    if not table.shape[1]:
        raise ValueError(f"No nonempty {dataset} samples remain")
    sample_ids = np.asarray(table.ids(axis="sample"), dtype=str)
    meta = meta.set_index("sample_id", drop=False).loc[sample_ids].reset_index(drop=True)
    meta["_matrix_index"] = np.arange(len(meta), dtype=int)
    matrix = (aggregate_genus(table, genus_map) if genus_map is not None else
              table.norm(axis="sample", inplace=False).matrix_data.T.tocsr())
    return meta, matrix

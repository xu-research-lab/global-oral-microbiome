"""Genus aggregation, cohort eligibility and leave-one-project-out splits."""
import numpy as np
import pandas as pd

from disease_common import DATA_DIR, INPUT_DIR, SITES, TAXONOMY_PATH, fold_name, load_metadata, manifest_path, split_path

MIN_GROUP_TOTAL = 30
MIN_MATCHED_PROJECTS = 2
FEATURE_THRESHOLD = 0.01


def parse_genus(taxon, include_unassigned=False):
    if pd.isna(taxon):
        return None
    return next((p.strip() for p in str(taxon).split(";")
                 if p.strip().startswith("g__") and (include_unassigned or len(p.strip()) > 3)), None)


def aggregate_to_genus(table, taxonomy, scheme):
    genus = taxonomy.Taxon.map(lambda x: parse_genus(x, include_unassigned=(scheme == "mean")))
    mapping = dict(zip(taxonomy.loc[genus.notna(), "Feature ID"].astype(str), genus.dropna()))
    keep = [f for f in table.ids(axis="observation") if f in mapping]
    if not keep:
        raise ValueError("No features map to a genus")
    table = table.filter(keep, axis="observation", inplace=False)
    return table.collapse(lambda feature, _: mapping[feature], axis="observation")


def matched_projects(meta, disease):
    counts = meta.groupby(["bioproject_id", "status_cleaned"]).sample_id.count().reset_index(name="n")
    sub = counts[counts.status_cleaned.isin(["HC", disease])]
    matched = sub.groupby("bioproject_id").status_cleaned.nunique()
    matched = matched[matched == 2].index.tolist()
    sub = sub[sub.bioproject_id.isin(matched)]
    totals = sub.groupby("status_cleaned").n.sum()
    if len(matched) < MIN_MATCHED_PROJECTS or any(totals.get(g, 0) < MIN_GROUP_TOTAL for g in ("HC", disease)):
        return []
    return matched


def select_genera(table, meta, scheme, threshold=FEATURE_THRESHOLD):
    features = np.asarray(table.ids(axis="observation"))
    samples = np.asarray(table.ids(axis="sample"))
    values = table.matrix_data.toarray()
    cohort = meta.set_index("sample_id").loc[samples, "bioproject_id"].astype(str).values
    if scheme == "prevalence":
        values = (values > 0).astype(np.float32)
    elif scheme != "mean":
        raise ValueError(f"Unknown filter: {scheme}")
    means = pd.DataFrame(values.T, index=cohort, columns=features).groupby(level=0).mean()
    return features[(means.max(axis=0) >= threshold).values]


def prepare_splits(scheme):
    import biom

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    table = aggregate_to_genus(biom.load_table(INPUT_DIR / "biom_oral_qc_dedup.biom"),
                               pd.read_csv(TAXONOMY_PATH, sep="\t"), scheme)
    meta = load_metadata()
    table.filter(meta.sample_id.values, axis="sample")
    table.remove_empty()
    samples = table.ids(axis="sample")
    table = table.filter(samples[table.sum(axis="sample") > 1000], axis="sample", inplace=False)
    table.remove_empty()
    meta = meta[meta.sample_id.isin(table.ids(axis="sample")) &
                meta.status_cleaned.notna() & (meta.status_cleaned != "Other") &
                meta.site_cleaned.isin(SITES)]
    folds = []
    for disease in sorted(meta.loc[meta.status_cleaned != "HC", "status_cleaned"].unique()):
        for site in SITES:
            local = meta[meta.site_cleaned == site]
            projects = matched_projects(local, disease)
            if not projects:
                continue
            picked = local[local.bioproject_id.isin(projects) & local.status_cleaned.isin(["HC", disease])]
            subset = table.filter(picked.sample_id.values, axis="sample", inplace=False)
            subset.norm(axis="sample", inplace=True)
            features = select_genera(subset, picked, scheme)
            if not len(features):
                continue
            subset = subset.filter(features, axis="observation", inplace=False)
            subset.norm(axis="sample", inplace=True)
            for project in projects:
                row = pd.Series(dict(site=site.replace(" ", "_"), disease=disease.replace(" ", "_"), study=project))
                name = fold_name(row)
                for partition in ("train", "test"):
                    keep = picked.bioproject_id == project
                    if partition == "train":
                        keep = ~keep
                    split = subset.filter(picked.loc[keep, "sample_id"].values, axis="sample", inplace=False)
                    with biom.util.biom_open(str(split_path(scheme, partition, name)), "w") as handle:
                        split.to_hdf5(handle, "disease classification")
                folds.append(row.to_dict())
                print(f"{scheme}: {name}; {len(features)} genera", flush=True)
    if not folds:
        raise ValueError("No eligible disease/site/project comparisons")
    pd.DataFrame(folds).to_csv(manifest_path(scheme), index=False)

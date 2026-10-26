"""Shared paths, metadata and fold names for disease classification."""
from pathlib import Path
import os

import pandas as pd

MODULE_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(os.environ.get("GLOBAL_ORAL_ROOT", MODULE_DIR.parent))
INPUT_DIR = PROJECT_ROOT / "data/metadata-and-biom-v2-qc"
RESULT_DIR = Path(os.environ.get("DISEASE_CLASSIFICATION_RESULT_DIR", MODULE_DIR / "result"))
DATA_DIR = Path(os.environ.get("DISEASE_CLASSIFICATION_DATA_DIR", MODULE_DIR / "data"))
TAXONOMY_PATH = Path(os.environ.get("GLOBAL_ORAL_TAXONOMY", "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"))
SITES = ("saliva", "tooth")
PLOT_SITES = SITES
RF_OPTIONS = dict(n_estimators=500, max_depth=5, max_features="log2", random_state=42)


def load_metadata(path=None):
    meta = pd.read_csv(path or INPUT_DIR / "metadata_oral_qc_dedup.tsv", sep="\t", low_memory=False)
    if meta.sample_id.duplicated().any():
        raise ValueError("Duplicate sample IDs in metadata")
    return meta


def load_labels():
    meta = load_metadata().set_index("sample_id")
    status = meta.status_cleaned
    valid = status.notna() & (status != "Other")
    return (status[valid] != "HC").astype(int).rename("group")


def fold_name(row):
    return f"{row.site}_{row.disease}_{row.study}"


def split_path(scheme, partition, name):
    prefix = "lso_genus" if scheme == "mean" else "lsoPgenus"
    return DATA_DIR / f"{prefix}_{partition}_{name}.biom"


def manifest_path(scheme):
    return DATA_DIR / f"run_info_{scheme}.csv"


def read_manifest(scheme):
    path = manifest_path(scheme)
    manifest = pd.read_csv(path, dtype=str)
    if manifest[["site", "disease", "study"]].isna().any().any() or manifest.duplicated().any():
        raise ValueError(f"Invalid or duplicate folds in {path}")
    return manifest


def predictions_path(scheme):
    suffix = "lso" if scheme == "mean" else "lsoP"
    return RESULT_DIR / f"RF_results_{suffix}.csv"

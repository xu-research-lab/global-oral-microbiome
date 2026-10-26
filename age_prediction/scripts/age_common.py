"""Shared paths, metadata, and prediction-table readers for age analyses."""

import os
from pathlib import Path
import re

import pandas as pd

MODULE_DIR = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("GLOBAL_ORAL_ROOT") or MODULE_DIR.parent)
RESULT = Path(os.environ.get("AGE_PREDICTION_RESULT_DIR") or MODULE_DIR / "result")
METADATA_PATH = ROOT / "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"
GENUS_TABLE_PATH = ROOT / "data/metadata-and-biom-v2-qc/genus_table_counts.tsv"
FIGURE_DIR = RESULT / "figures"
ORAL_SITES = ["saliva", "tooth", "tongue", "oral mucosa"]
SCHEMES = {
    "project": ("leave_studies_out_rf", "genus"),
    "random": ("random_split_rf", "random"),
}
PREDICTION_RE = re.compile(r"age_predict_test_(\d+)_(.+?)_(genus|random)_range_1_80\.csv")

DISEASE_ABBR = {'Aphthous': 'AS',
 'Atopic Dermatitis': 'AtD',
 'Axial Spondyloarthritis': 'axSpA',
 'Chagas Disease': 'CD',
 'Celiac Disease': 'Celiac Disease',
 'Colorectal Cancer': 'CRC',
 'COVID-19': 'COVID-19',
 'Dental Caries': 'DC',
 'Depression': 'Depression',
 'Gastroesophageal Reflux Disease': 'GERD',
 'Head and Neck Cancer': 'HNC',
 'HIV': 'HIV',
 'IBD': 'IBD',
 'ICU': 'ICU',
 'Lung Cancer': 'LC',
 'Obesity': 'Obesity',
 'Oral Lichen Planus': 'OLP',
 'Oral Mucositis': 'OM',
 'Oral Squamous Cell Carcinoma': 'OSCC',
 'Oral Submucous Fibrosis': 'OSMF',
 'Oral Verrucous Hyperplasia': 'OVH',
 'OSCC': 'OSCC',
 'Periodontal Disease': 'PD'}


def load_metadata(path=METADATA_PATH):
    meta = pd.read_csv(path, sep="\t", dtype={"sample_id": str}, low_memory=False)
    required = {"sample_id", "site_cleaned", "age_cleaned", "status_cleaned", "bioproject_id"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata missing columns: {sorted(missing)}")
    if meta["sample_id"].isna().any() or meta["sample_id"].duplicated().any():
        raise ValueError("Metadata sample_id values must be present and unique")
    return meta


def prediction_files(directory, token, site=None):
    files = []
    for path in sorted(Path(directory).glob(f"age_predict_test_*_{token}_range_1_80.csv")):
        match = PREDICTION_RE.fullmatch(path.name)
        if match and (site is None or match[2].replace("_", " ") == site):
            files.append(path)
    return files


def read_predictions(directory, token, site=None):
    frames = []
    for path in prediction_files(directory, token, site):
        match = PREDICTION_RE.fullmatch(path.name)
        frame = pd.read_csv(path, dtype={"sample_id": str})
        required = {"sample_id", "true_age", "pred_age"}
        if not required.issubset(frame.columns):
            raise ValueError(f"{path} missing columns: {sorted(required - set(frame.columns))}")
        frame["site"] = match[2].replace("_", " ")
        frame["fold"] = int(match[1])
        frames.append(frame)
    result = (pd.concat(frames, ignore_index=True) if frames else
              pd.DataFrame(columns=["sample_id", "true_age", "pred_age", "site", "fold"]))
    if result.duplicated(["site", "sample_id"]).any():
        raise ValueError(f"Repeated out-of-fold sample IDs in {directory}")
    return result


def site_predictions(site, metadata, directory, token="genus", *, sort_ids=False):
    predictions = read_predictions(directory, token, site)
    meta = metadata.loc[metadata["site_cleaned"] == site,
                        ["sample_id", "status_cleaned", "bioproject_id"]]
    frame = predictions.merge(meta, on="sample_id", how="inner", validate="many_to_one")
    if sort_ids:
        frame = frame.sort_values("sample_id").reset_index(drop=True)
    frame["delta_age"] = frame["pred_age"] - frame["true_age"]
    return frame


def status_groups(frame, min_n=10):
    counts = frame["status"].value_counts()
    keep = [status for status in counts.index if status == "HC" or counts[status] >= min_n]
    return frame[frame["status"].isin(keep)].copy(), keep

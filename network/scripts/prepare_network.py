"""Prepare paired project/status count tables and run cohort SparCC."""
from pathlib import Path
import os
import shutil
import subprocess

import numpy as np
import pandas as pd

from network_common import MODULE_DIR, RESULT_DIR, GROUPS

STATUS = {"HC": "Healthy Control", "PD": "Periodontal Disease", "DC": "Dental Caries", "OSCC": "OSCC"}


def select_projects(metadata):
    counts = (metadata.loc[metadata.site_cleaned.notna()]
              .groupby(["bioproject_id", "status_cleaned"]).size().unstack(fill_value=0)
              .reindex(columns=STATUS.values(), fill_value=0))
    return counts.loc[(counts[STATUS["HC"]] > 0) & (counts[list(STATUS.values())[1:]].sum(axis=1) > 0)]


def filter_unit(counts):
    if counts.shape[1] < 3:
        return None
    selected = counts.loc[(counts > 0).sum(axis=1) >= 5]
    selected = selected.loc[:, selected.sum(axis=0) > 0]
    return selected if selected.shape[0] >= 2 and selected.shape[1] >= 3 else None


def prepare_units(metadata_path=None, counts_path=None):
    if counts_path is None:
        raise ValueError("prepare requires --counts: the source TSV with bare genus IDs; the source is not bundled")
    metadata = pd.read_csv(metadata_path or MODULE_DIR / "data/metadata_after_cleaning.tsv", sep="\t", dtype=str)
    if metadata.sample_id.duplicated().any():
        raise ValueError("Duplicate metadata sample IDs")
    projects = select_projects(metadata)
    selected = metadata.loc[metadata.bioproject_id.isin(projects.index) & metadata.site_cleaned.notna()
                            & metadata.status_cleaned.isin(STATUS.values())]
    header = pd.read_csv(counts_path, sep="\t", nrows=0).columns
    samples = set(selected.sample_id)
    counts = pd.read_csv(counts_path, sep="\t", index_col=0,
                         usecols=[header[0], *[column for column in header[1:] if column in samples]])
    if counts.index.duplicated().any() or counts.index.astype(str).str.contains(";").any():
        raise ValueError("Counts require unique bare genus IDs, not taxonomy lineages")
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError("Counts must be finite nonnegative integers")
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    projects.to_csv(RESULT_DIR / "paired_projects.csv", sep="\t")
    records, study_lines = [], []
    for study in sorted(projects.index):
        for group in GROUPS:
            ids = selected.loc[(selected.bioproject_id == study) & (selected.status_cleaned == STATUS[group]), "sample_id"]
            if ids.empty:
                continue
            unit = filter_unit(counts.loc[:, [sample for sample in ids if sample in counts.columns]])
            path = RESULT_DIR / "biom_paired" / f"{group}_genus" / f"{study}.tsv"
            correlation = RESULT_DIR / "per_site_sparcc" / f"{group}_genus" / f"{study}.csv"
            if unit is not None:
                path.parent.mkdir(parents=True, exist_ok=True)
                text = unit.astype(np.int64).to_csv(sep="\t")
                if correlation.exists() and (not path.exists() or path.read_text() != text):
                    correlation.unlink()
                path.write_text(text)
                study_lines.append(f"{group}_genus {study}")
            elif path.exists():
                path.unlink()
                correlation.unlink(missing_ok=True)
            records.append(dict(unit_id=f"{group}_genus", category=group, study=study,
                                n_samples=len(ids), n_genera=0 if unit is None else len(unit),
                                kept=unit is not None,
                                path=str(path.relative_to(RESULT_DIR)) if unit is not None else ""))
    pd.DataFrame(records).to_csv(RESULT_DIR / "units.csv", index=False)
    (RESULT_DIR / "study_list_genus.txt").write_text("\n".join(study_lines) + "\n")
    print(f"[prepare] {len(projects)} paired projects; {len(study_lines)} units")


def run_sparcc(unit_index=None):
    lines = (RESULT_DIR / "study_list_genus.txt").read_text().splitlines()
    if unit_index is not None:
        if not 1 <= unit_index <= len(lines):
            raise ValueError(f"Unit index must be between 1 and {len(lines)}")
        lines = [lines[unit_index - 1]]
    rscript = os.environ.get("RSCRIPT", shutil.which("Rscript") or "Rscript")
    for line in lines:
        group, study = line.split()
        counts = RESULT_DIR / "biom_paired" / group / f"{study}.tsv"
        output = RESULT_DIR / "per_site_sparcc" / group / f"{study}.csv"
        if output.exists():
            print(f"[sparcc] cached: {group}/{study}")
            continue
        subprocess.run([rscript, str(Path(__file__).with_name("sparcc_cohort.R")),
                        str(counts), str(output)], check=True)

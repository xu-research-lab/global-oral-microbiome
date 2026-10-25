#!/usr/bin/env python3
"""Merge oral-only and gut-only project-level RF CV outputs, then refresh plots/report.

Expected layout:
  results/rf_oral_only_backup_20260728/project/
  continent_rf/result/oral/  # gut-only after --datasets gut
"""

from __future__ import annotations

import base64
import os
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("GLOBAL_ORAL_ROOT") or Path(__file__).resolve().parents[2])
MLRES = Path(__file__).resolve().parents[1] / "result"
SCRIPT = Path(__file__).resolve().parent
ORAL_BK = ROOT / "results/rf_oral_only_backup_20260728"
PROJ = MLRES / "oral"
GUT_DIR = MLRES / "gut"
PY = Path(sys.executable)

SITE_ORDER = ["Saliva", "Tongue", "Tooth", "Gut"]
REGION_ORDER = [
    "North America",
    "Asia",
    "Europe",
    "Oceania",
    "South America",
    "Africa",
]
SUMMARY_FILES = [
    "region_ova_summary.tsv",
    "region_ova_fold_metrics.tsv",
    "region_ova_predictions.tsv",
    "region_ovo_summary.tsv",
    "region_ovo_fold_metrics.tsv",
    "region_ovo_predictions.tsv",
]


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _write(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)


def merge_pair(oral_dir: Path, gut_dir: Path, out_dir: Path, label: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in SUMMARY_FILES:
        oral = _read(oral_dir / name)
        gut = _read(gut_dir / name)
        if "dataset" not in oral.columns or "dataset" not in gut.columns:
            raise ValueError(f"{name} missing dataset column")
        oral = oral[oral["dataset"] == "Oral"].copy()
        gut = gut[gut["dataset"] == "Gut"].copy()
        if oral.empty:
            raise ValueError(f"{oral_dir / name} has no Oral rows")
        if gut.empty:
            raise ValueError(f"{gut_dir / name} has no Gut rows")
        merged = pd.concat([oral, gut], ignore_index=True)
        _write(merged, out_dir / name)
        print(f"[{label}] {name}: oral={len(oral):,} gut={len(gut):,} -> {len(merged):,}")

    oral_cfg = json.loads((oral_dir / "run_config.json").read_text())
    gut_cfg = json.loads((gut_dir / "run_config.json").read_text())
    cfg = dict(oral_cfg)
    cfg["datasets"] = ["oral", "gut"]
    cfg["merged_from"] = {
        "oral_backup": str(oral_dir),
        "gut_run": str(gut_dir),
        "oral_config_datasets": oral_cfg.get("datasets"),
        "gut_config_datasets": gut_cfg.get("datasets"),
    }
    (out_dir / "run_config.json").write_text(json.dumps(cfg, indent=2, sort_keys=True) + "\n")


def run_plots() -> None:
    env = {
        **dict(**{k: v for k, v in __import__("os").environ.items()}),
        "MPLCONFIGDIR": "/tmp/mplconfig_xuhuan",
        "MPLBACKEND": "Agg",
    }
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    for script in ["plot_oral_gut_project_cv.py"]:
        print(f"Running {script}")
        subprocess.check_call([str(PY), str(SCRIPT / script)], cwd=str(ROOT), env=env)


def main() -> None:
    oral_proj = ORAL_BK / "project"
    if not (oral_proj / "region_ova_summary.tsv").exists():
        raise FileNotFoundError(f"Missing oral backup under {oral_proj}")
    if not (PROJ / "region_ova_summary.tsv").exists():
        raise FileNotFoundError(f"Missing gut/project summary under {PROJ}")

    # Stage gut outputs aside, then write merged tables into official dirs.
    gut_stage = GUT_DIR
    gut_proj = gut_stage
    gut_proj.mkdir(parents=True, exist_ok=True)

    for name in SUMMARY_FILES + ["run_config.json"]:
        src_p = PROJ / name
        if src_p.exists():
            (gut_proj / name).write_bytes(src_p.read_bytes())

    # Validate gut stage actually contains Gut.
    gcheck = _read(gut_proj / "region_ova_summary.tsv")
    if set(gcheck["dataset"].unique()) != {"Gut"}:
        # Maybe already merged; still ok if Gut present.
        if "Gut" not in set(gcheck["dataset"].unique()):
            raise ValueError(f"Expected Gut rows in {gut_proj}/region_ova_summary.tsv, got {set(gcheck['dataset'])}")

    merge_pair(oral_proj, gut_proj, PROJ, "project")
    run_plots()
    print("Merge + plots complete.")


if __name__ == "__main__":
    main()

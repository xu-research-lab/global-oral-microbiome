#!/bin/bash
#SBATCH -J ts_st2_pair_oral
#SBATCH -p cu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --exclude=cu01,cu02
#SBATCH -o logs/sourcetracker/st2.%A_%a.out
#SBATCH -e logs/sourcetracker/st2.%A_%a.err

set -euo pipefail

PROJECT_ROOT="/home/zhangcaiyi/project/oral_gut_pair_TimeSeries_11bp"
cd "$PROJECT_ROOT"

echo "Job started at: $(date)"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Node: $(hostname)"
echo "Workdir: $(pwd)"

MANIFEST="${MANIFEST:-$PROJECT_ROOT/04.sourcetracking/sourcetracker/st2_tasks_manifest.tsv}"

SOURCE_RAREFACTION_DEPTH="${SOURCE_RAREFACTION_DEPTH:-1000}"
SINK_RAREFACTION_DEPTH="${SINK_RAREFACTION_DEPTH:-1000}"
OVERWRITE="${OVERWRITE:-0}"
TASK_ID_OFFSET="${TASK_ID_OFFSET:-0}"

export MANIFEST
export TASK_ID_OFFSET
export SOURCE_RAREFACTION_DEPTH
export SINK_RAREFACTION_DEPTH
export OVERWRITE

echo "MANIFEST: $MANIFEST"
echo "SOURCE_RAREFACTION_DEPTH: $SOURCE_RAREFACTION_DEPTH"
echo "SINK_RAREFACTION_DEPTH: $SINK_RAREFACTION_DEPTH"
echo "OVERWRITE: $OVERWRITE"
echo "TASK_ID_OFFSET: $TASK_ID_OFFSET"

if [ ! -s "$MANIFEST" ]; then
  echo "ERROR: manifest not found or empty: $MANIFEST"
  exit 1
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID is empty. Submit this script with --array."
  exit 1
fi

module load miniconda/4.9.2 || true
source activate sourcetracker2 || true

echo "Python: $(which python)"
python --version
echo "SourceTracker2:"
which sourcetracker2 || true
sourcetracker2 --help | head -n 5 || true

eval "$(
python - <<'PY'
import os
import pandas as pd
from pathlib import Path
import shlex

manifest = os.environ["MANIFEST"]
task_id = int(os.environ["SLURM_ARRAY_TASK_ID"]) + int(os.environ.get("TASK_ID_OFFSET", "0"))

df = pd.read_csv(manifest, sep="\t")

if "task_id" not in df.columns:
    raise RuntimeError("Manifest missing column: task_id")

row = df.loc[pd.to_numeric(df["task_id"], errors="coerce") == task_id, :]

if row.empty:
    raise RuntimeError(f"No manifest row found for task_id={task_id}")

r = row.iloc[0]

required = ["input_biom", "map_tsv", "st2_out_dir"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise RuntimeError("Manifest missing columns: " + ",".join(missing))

def emit(name, value):
    print(f'{name}={shlex.quote(str(value))}')

emit("INPUT_BIOM", r["input_biom"])
emit("MAP_TSV", r["map_tsv"])
emit("ST2_OUT_DIR", r["st2_out_dir"])
emit("BIOPROJECT", r.get("bioproject", "NA"))
emit("PAIR_ID", r.get("pair_id", "NA"))
emit("GUT_SAMPLE_ID", r.get("gut_sample_id", "NA"))
PY
)"

echo "Task info:"
echo "  BIOPROJECT: $BIOPROJECT"
echo "  PAIR_ID: $PAIR_ID"
echo "  GUT_SAMPLE_ID: $GUT_SAMPLE_ID"
echo "  INPUT_BIOM: $INPUT_BIOM"
echo "  MAP_TSV: $MAP_TSV"
echo "  ST2_OUT_DIR: $ST2_OUT_DIR"

if [ ! -s "$INPUT_BIOM" ]; then
  echo "ERROR: input BIOM not found: $INPUT_BIOM"
  exit 1
fi

if [ ! -s "$MAP_TSV" ]; then
  echo "ERROR: map TSV not found: $MAP_TSV"
  exit 1
fi

# SourceTracker2 requires output_dir to not already exist.
if [ -s "$ST2_OUT_DIR/mixing_proportions.txt" ] && [ "$OVERWRITE" != "1" ]; then
  echo "mixing_proportions.txt already exists; skip task."
  echo "Existing output: $ST2_OUT_DIR/mixing_proportions.txt"
  exit 0
fi

if [ -d "$ST2_OUT_DIR" ]; then
  if [ "$OVERWRITE" = "1" ]; then
    echo "OVERWRITE=1, removing old ST2 output directory..."
    rm -rf "$ST2_OUT_DIR"
  else
    echo "Existing ST2 output directory found but no mixing_proportions.txt."
    echo "Removing stale/failed directory so SourceTracker2 can recreate it:"
    echo "  $ST2_OUT_DIR"
    rm -rf "$ST2_OUT_DIR"
  fi
fi

echo "Running SourceTracker2..."

sourcetracker2 gibbs \
  -i "$INPUT_BIOM" \
  -m "$MAP_TSV" \
  -o "$ST2_OUT_DIR" \
  --source_rarefaction_depth "$SOURCE_RAREFACTION_DEPTH" \
  --sink_rarefaction_depth "$SINK_RAREFACTION_DEPTH"

echo "Checking output..."
if [ -s "$ST2_OUT_DIR/mixing_proportions.txt" ]; then
  echo "SUCCESS: $ST2_OUT_DIR/mixing_proportions.txt"
else
  echo "ERROR: SourceTracker2 finished but mixing_proportions.txt was no
  echo "ST2
  ls -lh "$ST2_OUT_DIR" || 
  exit 1
fi



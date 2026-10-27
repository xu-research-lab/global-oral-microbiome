#!/usr/bin/env bash
#SBATCH --job-name=oral_gut_st2
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --exclude=cu01,cu02
#SBATCH --output=/home/zhangcaiyi/project/oral_gut_pair/logs/slurm/st2.%A_%a.out
#SBATCH --error=/home/zhangcaiyi/project/oral_gut_pair/logs/slurm/st2.%A_%a.err

set -euo pipefail

ROOT=/home/zhangcaiyi/project/oral_gut_pair
: "${MANIFEST:?MANIFEST must point to st2_tasks_manifest.tsv}"
: "${SLURM_ARRAY_TASK_ID:?Run through sbatch --array}"

if [[ ! "$(hostname)" =~ ^cu(0[3-9]|10)$ ]]; then
  echo "Refusing to run ST2 outside cu03-cu10: $(hostname)" >&2
  exit 2
fi
if [[ ! -s "$MANIFEST" ]]; then
  echo "Manifest missing or empty: $MANIFEST" >&2
  exit 2
fi

module load miniconda/4.9.2 >/dev/null 2>&1 || true

TASK_ID_OFFSET="${TASK_ID_OFFSET:-0}"
eval "$(
  python - "$MANIFEST" "$SLURM_ARRAY_TASK_ID" "$TASK_ID_OFFSET" <<'PY'
import shlex
import sys
import pandas as pd

manifest = sys.argv[1]
task_id = int(sys.argv[2]) + int(sys.argv[3])
table = pd.read_csv(manifest, sep="\t")
required = {"task_id", "input_biom", "map_tsv", "st2_out_dir"}
missing = required - set(table.columns)
if missing:
    raise RuntimeError(f"Manifest missing columns: {sorted(missing)}")
row = table.loc[pd.to_numeric(table["task_id"], errors="coerce").eq(task_id)]
if row.shape[0] != 1:
    raise RuntimeError(f"Expected one manifest row for task_id={task_id}; got {row.shape[0]}")
row = row.iloc[0]
for name, column in [
    ("INPUT_BIOM", "input_biom"),
    ("MAP_TSV", "map_tsv"),
    ("ST2_OUT_DIR", "st2_out_dir"),
]:
    print(f"{name}={shlex.quote(str(row[column]))}")
PY
)"

if [[ ! -s "$INPUT_BIOM" || ! -s "$MAP_TSV" ]]; then
  echo "Task input missing for array index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

if [[ -s "$ST2_OUT_DIR/mixing_proportions.txt" ]]; then
  echo "Existing successful output; skipping: $ST2_OUT_DIR"
  exit 0
fi
if [[ -d "$ST2_OUT_DIR" ]]; then
  rm -rf "$ST2_OUT_DIR"
fi

conda run -n sourcetracker2 sourcetracker2 gibbs \
  -i "$INPUT_BIOM" \
  -m "$MAP_TSV" \
  -o "$ST2_OUT_DIR" \
  --source_rarefaction_depth "${SOURCE_RAREFACTION_DEPTH:-1000}" \
  --sink_rarefaction_depth "${SINK_RAREFACTION_DEPTH:-1000}"

test -s "$ST2_OUT_DIR/mixing_proportions.txt"
echo "ST2_SUCCESS task=${SLURM_ARRAY_TASK_ID} node=$(hostname)"

#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/zhangcaiyi/project/oral_gut_pair
RUNNER="$ROOT/scripts/slurm/run_st2_array_task.sh"
MAX_CONCURRENT="${MAX_CONCURRENT:-100}"
MAX_ARRAY_TASKS="${MAX_ARRAY_TASKS:-1000}"

if [[ "$(hostname)" != "cu01" ]]; then
  echo "ST2 arrays must be submitted from cu01; current node: $(hostname)" >&2
  exit 2
fi

mkdir -p "$ROOT/logs/slurm"

entries=(
  "single_pair|$ROOT/results/single_pair/04_source_tracking/sourcetracker/st2_tasks_manifest.tsv"
  "time_series|$ROOT/results/time_series/04_source_tracking/sourcetracker/st2_tasks_manifest.tsv"
  "tooth_vs_saliva|$ROOT/results/single_pair/tooth_vs_saliva/04_source_tracking/st2_tasks_manifest.tsv"
)
for entry in "${entries[@]}"; do
  cohort="${entry%%|*}"
  manifest="${entry#*|}"
  if [[ ! -s "$manifest" ]]; then
    echo "Missing manifest for $cohort: $manifest" >&2
    exit 2
  fi
  n_tasks="$(
    module load miniconda/4.9.2 >/dev/null 2>&1 || true
    python - "$manifest" <<'PY'
import sys
import pandas as pd
table = pd.read_csv(sys.argv[1], sep="\t")
print(table["task_id"].nunique())
PY
  )"
  if [[ "$n_tasks" -lt 1 ]]; then
    echo "No ST2 tasks for $cohort" >&2
    exit 2
  fi
  offset=0
  while [[ "$offset" -lt "$n_tasks" ]]; do
    remaining=$((n_tasks - offset))
    chunk_size="$MAX_ARRAY_TASKS"
    if [[ "$remaining" -lt "$chunk_size" ]]; then
      chunk_size="$remaining"
    fi
    job_id="$(
      sbatch --parsable \
        --job-name="st2_${cohort}" \
        --partition=cu \
        --nodes=1 \
        --exclude=cu01,cu02 \
        --array="1-${chunk_size}%${MAX_CONCURRENT}" \
        --export="ALL,MANIFEST=${manifest},TASK_ID_OFFSET=${offset}" \
        "$RUNNER"
    )"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$(date --iso-8601=seconds)" "$cohort" "$n_tasks" "$offset" "$chunk_size" "$job_id" \
      >> "$ROOT/logs/slurm/st2_submissions.tsv"
    echo "$cohort: submitted task IDs $((offset + 1))-$((offset + chunk_size)) as job $job_id"
    offset=$((offset + chunk_size))
  done
done

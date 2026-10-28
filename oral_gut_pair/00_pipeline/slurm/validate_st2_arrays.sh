#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/zhangcaiyi/project/oral_gut_pair

if [[ "$(hostname)" != "cu01" ]]; then
  echo "ST2 validation must run on cu01; current node: $(hostname)" >&2
  exit 2
fi

module load miniconda/4.9.2 >/dev/null 2>&1 || true
outputs=(
  "$ROOT/results/single_pair/04_source_tracking/sourcetracker"
  "$ROOT/results/time_series/04_source_tracking/sourcetracker"
  "$ROOT/results/single_pair/tooth_vs_saliva/04_source_tracking"
)
for out in "${outputs[@]}"; do
  python "$ROOT/scripts/python/oral_gut_pair/validate_st2_completion.py" \
    --manifest "$out/st2_tasks_manifest.tsv" \
    --marker "$out/st2_array_complete.ok" \
    --status "$out/st2_array_completion_status.tsv"
done

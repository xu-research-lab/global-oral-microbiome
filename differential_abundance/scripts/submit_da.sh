#!/usr/bin/env bash
# Submit independent method fits, then finalize only after all three succeed.
set -euo pipefail
analysis="${1:-}"
case "$analysis" in age|continent|site|sex) ;; *)
  echo "Usage: bash differential_abundance/scripts/submit_da.sh <age|continent|site|sex>" >&2
  exit 2 ;;
esac
project_root="${GLOBAL_ORAL_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)}"
export GLOBAL_ORAL_ROOT="$project_root"
cd "$project_root"
log_dir="$project_root/differential_abundance/logs"
mkdir -p "$log_dir"
job_script="$project_root/differential_abundance/scripts/run_da.sbatch"
job_ids=()
for method in linda maaslin2 ancombc2; do
  job_id=$(sbatch --parsable --export=ALL --job-name="da_${analysis}_${method}" \
    --output="$log_dir/%x_%j.out" --error="$log_dir/%x_%j.err" \
    "$job_script" "$analysis" --stage=fit --method="$method")
  job_ids+=("${job_id%%;*}")
done
dependency=$(IFS=:; echo "${job_ids[*]}")
sbatch --parsable --export=ALL --dependency="afterok:$dependency" \
  --job-name="da_${analysis}_finalize" --output="$log_dir/%x_%j.out" \
  --error="$log_dir/%x_%j.err" "$job_script" "$analysis" --stage=finalize

#!/bin/bash
# ==============================================================================
# submit_age_da_ova.sh
# One-shot submission of the parallel SALIVA age OVA DA run (each valid age
# group vs the pooled "Other" of the other valid groups, da_v2, HC SCOPE ONLY):
#   1. ANCOMBC2 array (6 tasks, one OVA comparison per task, --only within hc)
#   2. LinDA + Maaslin2 (hc scope)
#   3. Merge + consensus (hc scope) — submitted with --dependency=afterok
#      on 1-2 so it only starts once every comparison is done.
#      (figures come later once 03-plot-saliva-ova-bubble.R is finalized)
#
# Run from project root:
#   bash differential_abundance/age_differential_ova/submit_age_da_ova.sh
# ==============================================================================
set -uo pipefail
cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."

AGE=differential_abundance/age_differential_ova

echo "=== [1/3] Submitting ANCOMBC2 array (hc scope, 1 OVA comparison per task) ==="
A=$(sbatch --parsable $AGE/run_ancombc2_array.sbatch)
echo "  ancombc2_array       : $A"

echo "=== [2/3] Submitting LinDA + Maaslin2 (hc scope) ==="
B=$(sbatch --parsable $AGE/run_linda_maaslin2.sbatch)
echo "  linda_maaslin2       : $B"

echo "=== [3/3] Submitting merge + consensus (afterok:$A:$B) ==="
M=$(sbatch --parsable --dependency=afterok:$A:$B $AGE/run_merge_consensus.sbatch)
echo "  merge_consensus      : $M"

echo ""
echo "============================================================"
echo "ALL AGE OVA DA JOBS SUBMITTED (saliva, HC only; target vs Other; robust = >=2 methods)"
echo "  ANCOMBC2 array       : $A"
echo "  linda + maaslin2     : $B"
echo "  merge + consensus    : $M  (waits for $A:$B via afterok)"
echo "============================================================"
echo "Monitor:  squeue -u \$USER"
echo "Logs:     results/logs/da_ageova_*"
echo "Results:  differential_abundance/result/age_hc_ova/"

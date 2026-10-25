#!/bin/bash
# ==============================================================================
# submit_sex_da_continent.sh
# One-shot submission of the parallel sex-within-continent DA run
# (per-continent F_vs_M, saliva, HC only, da_v2):
#   1. ANCOMBC2 array (5 tasks: 5 continents, hc)
#   2. LinDA + Maaslin2 (hc)
#   3. Merge + consensus + figures (hc) — submitted with
#      --dependency=afterok on 1-2 so it only starts once every comparison
#      is done.
#
# Run from project root:
#   bash differential_abundance/sex_differential_continent/submit_sex_da_continent.sh
# ==============================================================================
set -uo pipefail
cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."

SEX=differential_abundance/sex_differential_continent

echo "=== [1/3] Submitting ANCOMBC2 array (5 tasks: 5 continents, hc) ==="
A=$(sbatch --parsable $SEX/run_ancombc2_array.sbatch)
echo "  ancombc2_array       : $A"

echo "=== [2/3] Submitting LinDA + Maaslin2 (hc) ==="
B=$(sbatch --parsable $SEX/run_linda_maaslin2.sbatch)
echo "  linda_maaslin2       : $B"

echo "=== [3/3] Submitting merge + consensus (afterok:$A:$B) ==="
M=$(sbatch --parsable --dependency=afterok:$A:$B $SEX/run_merge_consensus.sbatch)
echo "  merge_consensus      : $M"

echo ""
echo "============================================================"
echo "ALL SEX-WITHIN-CONTINENT DA JOBS SUBMITTED (saliva, hc; robust = >=2 methods)"
echo "  ANCOMBC2 array       : $A"
echo "  linda + maaslin2     : $B"
echo "  merge + consensus    : $M  (waits for $A:$B via afterok)"
echo "============================================================"
echo "Monitor:  squeue -u \$USER"
echo "Logs:     results/logs/da_sex_cont_*"
echo "Results:  results/da_v2/sex_continent_hc/  (legacy, not moved)"

#!/bin/bash
# ==============================================================================
# submit_sex_pooled.sh
# One-shot submission of the POOLED-continent sex DA (saliva F_vs_M, ALL
# continents merged, continent_group as covariate, da_v2, HC scope ONLY):
#   1. LinDA + Maaslin2 (single job)
#   2. ANCOMBC2 (single job — one pooled comparison, no array needed)
#   3. Merge + consensus — submitted with --dependency=afterok on 1-2 so it
#      only starts once both are done.
#
# Smoke first (recommended; lands in differential_abundance/result/da_v2_smoke_pooled):
#   sbatch differential_abundance/sex_differential_continent/pooled/run_smoke.sbatch
#
# Full run (from project root):
#   bash differential_abundance/sex_differential_continent/pooled/submit_sex_pooled.sh
# ==============================================================================
set -uo pipefail
cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../../.."

POOLED=differential_abundance/sex_differential_continent/pooled

echo "=== [1/3] Submitting LinDA + Maaslin2 (pooled, hc) ==="
B=$(sbatch --parsable $POOLED/run_linda_maaslin2.sbatch)
echo "  linda_maaslin2 : $B"

echo "=== [2/3] Submitting ANCOMBC2 (pooled, hc) ==="
A=$(sbatch --parsable $POOLED/run_ancombc2.sbatch)
echo "  ancombc2       : $A"

echo "=== [3/3] Submitting merge + consensus (afterok:$A:$B) ==="
M=$(sbatch --parsable --dependency=afterok:$A:$B $POOLED/run_merge_consensus.sbatch)
echo "  merge_consensus: $M"

echo ""
echo "============================================================"
echo "ALL POOLED SEX-DA JOBS SUBMITTED (saliva F_vs_M, HC only;"
echo "robust = significant in >=2 methods)"
echo "  linda + maaslin2 : $B"
echo "  ancombc2         : $A"
echo "  merge + consensus: $M  (waits for $A:$B via afterok)"
echo "============================================================"
echo "Monitor:  squeue -u \$USER"
echo "Logs:     results/logs/da_sex_pooled_*"
echo "Results:  differential_abundance/result/sex_pooled_hc/"
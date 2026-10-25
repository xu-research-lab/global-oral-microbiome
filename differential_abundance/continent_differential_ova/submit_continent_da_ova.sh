#!/bin/bash
# ==============================================================================
# submit_continent_da_ova.sh
# One-shot submission of the parallel continent OVA DA run (each valid continent
# vs the pooled "Other" of the other valid continents, da_v2, saliva only), for
# BOTH scopes (all + hc). For each scope:
#   1. ANCOMBC2 array (1 task per OVA comparison, --only within scope)
#   2. LinDA + Maaslin2
#   3. Merge + consensus + figures — submitted with --dependency=afterok
#      on 1-2 (same scope) so it only starts once every comparison is done.
#
# Run from project root:
#   bash differential_abundance/continent_differential_ova/submit_continent_da_ova.sh
# ==============================================================================
set -uo pipefail
cd "$(cd "$(command dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../.."

S=differential_abundance/continent_differential_ova

for SCOPE in all hc; do
  echo "############### SCOPE=$SCOPE ###############"
  echo "=== [$SCOPE 1/3] Submitting ANCOMBC2 array ==="
  A=$(sbatch --parsable --export=SCOPE=$SCOPE $S/run_ancombc2_array.sbatch)
  echo "  ancombc2_array ($SCOPE)  : $A"

  echo "=== [$SCOPE 2/3] Submitting LinDA + Maaslin2 ==="
  B=$(sbatch --parsable --export=SCOPE=$SCOPE $S/run_linda_maaslin2.sbatch)
  echo "  linda_maaslin2 ($SCOPE)  : $B"

  echo "=== [$SCOPE 3/3] Submitting merge + consensus + figures (afterok:$A:$B) ==="
  M=$(sbatch --parsable --export=SCOPE=$SCOPE --dependency=afterok:$A:$B $S/run_merge_consensus.sbatch)
  echo "  merge_consensus ($SCOPE): $M  (waits for $A:$B via afterok)"
  echo ""
done

echo "============================================================"
echo "ALL CONTINENT OVA DA JOBS SUBMITTED (all + hc; saliva only; target vs Other; robust = >=2 methods)"
echo "Monitor:  squeue -u \$USER"
echo "Logs:     results/logs/da_contova_*"
echo "Results:  differential_abundance/result/continent_ova/  (all) | results/da_v2/continent_hc_ova/  (hc, legacy)"
echo "============================================================"

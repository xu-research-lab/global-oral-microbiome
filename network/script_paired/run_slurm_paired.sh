#!/bin/bash
# Paired-project HC/PD/DC/OSCC all-site genus SparCC batch job (slurm array).
#
# Each array task handles one (group_genus, study) unit from
#   network/result_paired/study_list_genus.txt
# (produced by script_hc_pd_dc_oscc/prepare_paired_units.py; projects with
# >=1 disease sample AND >=1 HC sample, split into (project,status) units)
# and writes
#   network/result_paired/per_site_sparcc/<grp>_genus/<study_id>.csv
#
# Submit from network/ (array size set to the study-list length):
#   sbatch --array=1-$(wc -l < result_paired/study_list_genus.txt) \
#       script_hc_pd_dc_oscc/run_slurm_paired.sh
# Re-submitting is safe: finished matrices are skipped (cache).

#SBATCH --job-name=sparcc_paired
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p cu
#SBATCH --mem=8G
#SBATCH -o output_paired/sparcc_paired_%a.out
#SBATCH -e output_paired/sparcc_paired_%a.err
#SBATCH --exclude=cu01

RSCRIPT=${RSCRIPT:-/home/zhangyuchen/.conda/envs/r_latest/bin/Rscript}
[ -x "$RSCRIPT" ] || { echo "ERROR: conda Rscript not found: $RSCRIPT" >&2; exit 1; }

# slurm runs a spool copy of this script; use SLURM_SUBMIT_DIR, not $0
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$SCRIPT_DIR" || { echo "ERROR: cannot cd to $SCRIPT_DIR" >&2; exit 1; }
mkdir -p output_paired

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" result_paired/study_list_genus.txt | tr -d '\r\n')
[ -z "$line" ] && exit 0

unit_id=$(echo "$line" | awk '{print $1}')
study_id=$(echo "$line" | awk '{print $2}')

biom="result_paired/biom_paired/${unit_id}/${study_id}.tsv"
out="result_paired/per_site_sparcc/${unit_id}/${study_id}.csv"
[ -f "$biom" ] || { echo "missing $biom" >&2; exit 1; }
[ -f "$out" ] && { echo "exists, skip: $out"; exit 0; }

"$RSCRIPT" script_hc_pd_dc_oscc/sparcc_cohort.R "$biom" "$out"
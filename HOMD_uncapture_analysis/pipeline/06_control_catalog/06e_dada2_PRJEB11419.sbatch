#!/bin/bash
#SBATCH --job-name=dada2_P11419
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH --nodelist=cu08
#SBATCH --mem=60G
#SBATCH -o /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/dada2_PRJEB11419.out
#SBATCH -e /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/dada2_PRJEB11419.err

# ============================================================
# DADA2 on PRJEB11419 (921 control samples, 28 threads)
# 7 metagenomic outliers excluded
# ============================================================

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

# Force R to use qiime2 libraries (not ~/R/4.5 incompatible)
export R_PROFILE_USER="${HOME}/.Rprofile.qiime2"

STUDY_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/06_control_catalog/work/control_dada2/PRJEB11419"
SCRIPT_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/01_data_assembly"

mkdir -p "${STUDY_DIR}/intermediate" "${STUDY_DIR}/results" "${STUDY_DIR}/temp"

cd "${STUDY_DIR}"

echo "=== DADA2 PRJEB11419 ==="
echo "Started: $(date)"
echo "Samples: $(wc -l < run_sample.txt)"
echo "Threads: ${SLURM_CPUS_PER_TASK:-28}"

Rscript "${SCRIPT_DIR}/01b_process_forwards.R" "${STUDY_DIR}" trim

echo "=== Finished: $(date) ==="
if [ -f "results/ASVs.fa" ]; then
    N_ASV=$(wc -l < results/ASVs.fa)
    echo "ASVs detected: $((N_ASV / 2))"
    ls -lh results/
fi

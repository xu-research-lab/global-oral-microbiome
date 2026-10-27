#!/bin/bash
#SBATCH --job-name=dada2_ctrl
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p cu
#SBATCH --nodelist=cu09
#SBATCH --mem=30G
#SBATCH -o /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/dada2_ctrl_%a.out
#SBATCH -e /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/dada2_ctrl_%a.err
#SBATCH --array=1-17

# ============================================================
# Run DADA2 on control samples (per BioProject)
# Adapted from process_data.sh for control samples
# ============================================================

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2   # has R + dada2 (R 4.1.2, ~/R/ is R 4.5 incompatible)

# Create .Rprofile to force qiime2 library path
# (R's system config hardcodes ~/R/4.5 first, .Rprofile overrides)
cat > "${HOME}/.Rprofile.qiime2" << 'REOF'
.libPaths(c(
  "/home/software/anaconda3/envs/qiime2-2022.2/lib/R/library",
  "/beegfs/software/R-packages"
))
REOF
export R_PROFILE_USER="${HOME}/.Rprofile.qiime2"

BASE_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/control_dada2"
SCRIPT_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis"

study_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BASE_DIR}/study_id.txt" | tr -d '\r\n')
file_path="${BASE_DIR}/${study_id}"

echo "=== DADA2 control sample processing ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "BioProject: ${study_id}"
echo "Working dir: ${file_path}"
echo "Started at: $(date)"

mkdir -p "${file_path}/intermediate" "${file_path}/results" "${file_path}/temp"

cd "${file_path}"

# Run DADA2 (with trimLeft=24 as per original pipeline)
Rscript "${SCRIPT_DIR}/process_forwards.R" "${file_path}" trim

echo "=== DADA2 finished for ${study_id} at $(date) ==="

# Quick summary
if [ -f "${file_path}/results/ASVs.fa" ]; then
    N_ASV=$(grep -c "^>" "${file_path}/results/ASVs.fa")
    echo "ASVs detected: ${N_ASV}"
    ls -lh "${file_path}/results/"
else
    echo "WARNING: ASVs.fa not found!"
fi

#!/bin/bash
#SBATCH --job-name=dada2_extra
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH --mem=60G
#SBATCH -o /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/dada2_extra_%a.out
#SBATCH -e /beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/logs/dada2_extra_%a.err
#SBATCH --array=1-7

# ============================================================
# DADA2 on extra control samples: negative_control + reagent_water
# ============================================================

set -e

source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2
export R_PROFILE_USER="${HOME}/.Rprofile.qiime2"

FASTQ_DUMP="/home/huangxc/software/sratoolkit.2.11.0-centos_linux64/bin/fastq-dump"
SRA_DIR="/beegfs/dongbiao/oral/control_samples/sra"
BASE_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/06_control_catalog/work/control_dada2_extra"
SCRIPT_DIR="/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/01_data_assembly"

# Build task list: bioproject<TAB>sample_id
python3 -c "
import pandas as pd, os
ctrl = pd.read_csv('/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline/06_control_catalog/results/control_sample_ids.tsv', sep='\t')
target = ctrl[ctrl['category'].isin(['negative_control', 'reagent_water'])]
sra_avail = set(os.listdir('$SRA_DIR'))
tasks = []
for _, row in target.iterrows():
    sid = row['sample_id']
    if sid in sra_avail:
        tasks.append((row['bioproject'], sid, row['category']))
with open('$BASE_DIR/task_list.txt', 'w') as f:
    for bp, sid, cat in tasks:
        f.write(f'{bp}\t{sid}\t{cat}\n')
print(f'Task list: {len(tasks)} samples')
" 2>&1

# Get this array task's BioProject
study_id=$(awk -F'\t' '{print $1}' "$BASE_DIR/task_list.txt" | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p" | tr -d '\r\n')
file_path="${BASE_DIR}/${study_id}"

echo "=== DADA2 extra controls ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "BioProject: ${study_id}"
echo "Started: $(date)"

mkdir -p "${file_path}/rawdata" "${file_path}/intermediate" "${file_path}/results" "${file_path}/temp"

# ---- Step 1: Convert SRA to FASTQ ----
echo "[convert] Converting SRA files for ${study_id}..."
awk -F'\t' -v bp="$study_id" '$1 == bp {print $2}' "$BASE_DIR/task_list.txt" | while read sid; do
    sra_path="${SRA_DIR}/${sid}/${sid}.sra"
    if [ ! -f "$sra_path" ]; then
        echo "  WARNING: SRA missing for $sid, skipping"
        continue
    fi

    # Skip if FASTQ already exists
    if [ -f "${file_path}/rawdata/${sid}.fastq.gz" ] || [ -f "${file_path}/rawdata/${sid}_1.fastq.gz" ]; then
        echo "  [skip] $sid"
        continue
    fi

    echo "  [convert] $sid"
    "$FASTQ_DUMP" --split-files --gzip --skip-technical --outdir "${file_path}/rawdata" "$sra_path" 2>&1 | tail -1

    # Rename for R script compatibility
    if [ -f "${file_path}/rawdata/${sid}_1.fastq.gz" ] && [ -f "${file_path}/rawdata/${sid}_2.fastq.gz" ]; then
        # Paired-end: keep _1/_2 naming (R script expects this)
        :
    elif [ -f "${file_path}/rawdata/${sid}_1.fastq.gz" ]; then
        # Single-end with _1 suffix -> rename to .fastq.gz
        mv "${file_path}/rawdata/${sid}_1.fastq.gz" "${file_path}/rawdata/${sid}.fastq.gz"
    fi
done

# ---- Step 2: Create run_sample.txt ----
ls "${file_path}/rawdata/"*.fastq.gz 2>/dev/null | while read f; do basename "$f"; done > "${file_path}/run_sample.txt"
echo "[convert] $(wc -l < ${file_path}/run_sample.txt) FASTQ files ready"

# ---- Step 3: Run DADA2 ----
cd "${file_path}"
echo "[dada2] Starting at $(date)"
Rscript "${SCRIPT_DIR}/01b_process_forwards.R" "${file_path}" trim

echo "=== DADA2 finished for ${study_id} at $(date) ==="
if [ -f "results/ASVs.fa" ]; then
    N_ASV=$(wc -l < results/ASVs.fa)
    echo "ASVs detected: $((N_ASV / 2))"
    ls -lh results/
fi

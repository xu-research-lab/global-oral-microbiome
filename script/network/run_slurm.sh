#!/bin/bash
#SBATCH --job-name=sparcc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p cu
#SBATCH --mem=8G
#SBATCH -o output/sparcc_%a.out
#SBATCH -e output/sparcc_%a.err
#SBATCH --array=1-41
#SBATCH --exclude=cu01

module load R/4.1.3

study_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" study_list_saliva.txt | tr -d '\r\n')

Rscript run_SpiecEasi.R ${study_id}


#!/bin/bash
#SBATCH --job-name=dada2
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p cu
#SBATCH --mem=10G
#SBATCH -o output/dada2_%a.out
#SBATCH -e output/dada2_%a.err
#SBATCH --array=345
#SBATCH --exclude=cu01

module load R/4.1.3

# read the study id and its paths
study_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" study_id.txt | tr -d '\r\n')
# study_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /home/dongbiao/word_embedding_microbiome/all_data/oral/global-oral-microbiome/mismatch_project_list.txt | tr -d '\r\n')
file_path="/beegfs/dongbiao/oral/${study_id}"
path=/home/dongbiao/software/parallel-20250222/bin

# makedirs
mkdir ${file_path}/sra ${file_path}/intermediate ${file_path}/rawdata ${file_path}/results ${file_path}/temp

# download rawdata
# prefetch --option-file ${file_path}/download.txt --output-directory ${file_path}/sra

# ls ${file_path}/sra/*/*.sra | awk -F'/' '{print $8}' | awk -F'.' '{print $1}' > ${file_path}/sample.txt

# threads=28
# run fastq-dump in parallel
# cat "${file_path}/sample.txt" | ${path}/parallel -j ${threads} --bar \
#     "fastq-dump --split-3 ${file_path}/sra/{}/{}.sra --outdir ${file_path}/rawdata"

# Run DADA2
ls ${file_path}/rawdata/ > ${file_path}/run_sample.txt

Rscript process_forwards.R ${file_path} trim
# Rscript process_forwards.R ${file_path}
# rm -rf ${file_path}/rawdata ${file_path}/intermediate ${file_path}/temp ${file_path}/run_sample.txt

#!/bin/bash

#SBATCH --job-name=cluster_features
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p cu
#SBATCH --mem=30G
#SBATCH -o output/cluster_features_%a.out
#SBATCH -e output/cluster_features_%a.err
#SBATCH --array=176
#SBATCH --exclude=cu01,cu05,cu06,cu07,cu08,cu09,cu10

# module load miniconda/4.9.2
# source activate qiime2-2022.2
source /home/hexincheng/miniconda3/bin/activate qiime2-2024.2

study_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /beegfs/dongbiao/oral/project_list.txt | tr -d '\r\n')
file_path="/beegfs/dongbiao/oral/${study_id}"
results_path="${file_path}/results"

echo "Processing study: $study_id"
echo "Working directory: $file_path"
echo "Results directory: $results_path"

echo "Step 1: generate representative sequences from DADA2.biom..."
python - "${results_path}" <<'PY'
import biom
import os
import sys

results_path = sys.argv[1]
table = biom.load_table(os.path.join(results_path, "DADA2.biom"))
seqs = table.ids(axis="observation")

with open(os.path.join(results_path, "rep-seqs.fasta"), "w") as f:
    for seq in seqs:
        f.write(f">{seq}\n{seq}\n")
PY

echo "Step 2: import representative sequences into QIIME2..."
qiime tools import \
  --input-path ${results_path}/rep-seqs.fasta \
  --output-path ${results_path}/rep-seqs.qza \
  --type 'FeatureData[Sequence]'

echo "Step 3: import DADA2 feature table into QIIME2..."
qiime tools import \
  --input-path ${results_path}/DADA2.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path ${results_path}/table.qza

echo "Step 4: run closed-reference clustering..."
qiime vsearch cluster-features-closed-reference \
    --i-sequences ${results_path}/rep-seqs.qza \
    --i-table ${results_path}/table.qza \
    --i-reference-sequences /home/hexincheng/projects/hxc/oral_project/closed_reference/2024.09.backbone.full-length.fna.qza \
    --p-perc-identity 0.99 \
    --p-threads 12 \
    --p-strand both \
    --o-clustered-table ${results_path}/clustered-table.qza \
    --o-clustered-sequences ${results_path}/clustered-sequences.qza \
    --o-unmatched-sequences ${results_path}/unmatched-sequences.qza

echo "Step 5: export clustering results..."
qiime tools export \
  --input-path ${results_path}/clustered-table.qza \
  --output-path ${results_path}/extracted-feature-table

echo "Closed-reference clustering completed."
echo "Output files:"
echo "  - ${results_path}/clustered-table.qza"
echo "  - ${results_path}/clustered-sequences.qza"
echo "  - ${results_path}/unmatched-sequences.qza"
echo "  - ${results_path}/extracted-feature-table/"

echo "Study ${study_id} closed-reference workflow finished."

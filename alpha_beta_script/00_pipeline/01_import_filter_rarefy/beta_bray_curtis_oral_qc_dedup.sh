#!/bin/bash
#SBATCH --job-name=beta_oral_qc_dedup
#SBATCH -N 1
#SBATCH --cpus-per-task=28
#SBATCH -p cu
#SBATCH --mem=250G
#SBATCH -o output/beta_oral_qc_dedup_%a.out
#SBATCH -e output/beta_oral_qc_dedup_%a.err
#SBATCH --exclude=cu01

source /home/hexincheng/miniconda3/bin/activate /home/hexincheng/miniconda3/envs/rachis-qiime2-2026.4
threads=28
export OMP_NUM_THREADS=${threads}

work_dir=/beegfs/dongbiao/oral/beta_oral_qc_dedup
metadata_path=/home/hexincheng/projects/hxc/oral_project/data/metadata_oral_qc_dedup.tsv
metadata_tsv=${work_dir}/metadata.tsv
biom_path=/home/hexincheng/projects/hxc/oral_project/data/biom_oral_qc_dedup.biom

mkdir -p ${work_dir}

python - <<PY
import pandas as pd

metadata = pd.read_csv("${metadata_path}", dtype={"sample_id": str}, sep="\t", low_memory=False)
metadata = metadata.rename(columns={"sample_id": "#SampleID"})
metadata = metadata.set_index("#SampleID")
metadata.to_csv("${metadata_tsv}", sep="\t")
PY

qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path ${biom_path} \
  --input-format BIOMV210Format \
  --output-path ${work_dir}/feature-table.qza

qiime feature-table filter-samples \
  --i-table ${work_dir}/feature-table.qza \
  --m-metadata-file ${metadata_tsv} \
  --o-filtered-table ${work_dir}/feature-table-filtered.qza

qiime feature-table filter-features \
  --i-table ${work_dir}/feature-table-filtered.qza \
  --p-min-samples 5 \
  --p-min-frequency 1 \
  --o-filtered-table ${work_dir}/feature-table-filtered-prevalence5.qza

qiime tools export \
  --input-path ${work_dir}/feature-table-filtered-prevalence5.qza \
  --output-path ${work_dir}/exported_table_prevalence5

qiime feature-table rarefy \
  --i-table ${work_dir}/feature-table-filtered-prevalence5.qza \
  --p-sampling-depth 1000 \
  --o-rarefied-table ${work_dir}/rarefied-table.qza

qiime diversity beta \
  --i-table ${work_dir}/rarefied-table.qza \
  --p-metric braycurtis \
  --p-n-jobs ${threads} \
  --o-distance-matrix ${work_dir}/bray-curtis-distance-matrix.qza

qiime diversity pcoa \
  --i-distance-matrix ${work_dir}/bray-curtis-distance-matrix.qza \
  --p-number-of-dimensions 3 \
  --o-pcoa ${work_dir}/bray-curtis-pcoa.qza

qiime emperor plot \
  --i-pcoa ${work_dir}/bray-curtis-pcoa.qza \
  --m-metadata-file ${metadata_tsv} \
  --o-visualization ${work_dir}/bray-curtis-emperor.qzv

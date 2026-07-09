#!/bin/bash
#SBATCH --job-name=alpha_div
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p cu
#SBATCH --mem=10G
#SBATCH -o output/alpha_div_%a.out
#SBATCH -e output/alpha_div_%a.err
#SBATCH --array=151,310
#SBATCH --exclude=cu01

module load miniconda/4.9.2
source activate
conda activate qiime2-2022.2
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

study_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /beegfs/dongbiao/oral/project_list.txt | tr -d '\r\n')

mkdir -p /beegfs/dongbiao/oral/${study_id}/results/alpha
path=/beegfs/dongbiao/oral/${study_id}/results/alpha
merged_metadata=/home/hexincheng/projects/hxc/oral_project/data/metadata-v2-merge-cleaned.csv

python "${script_dir}/filter_metadata_for_QIIME2_merged.py" \
    ${merged_metadata} \
    /beegfs/dongbiao/oral/${study_id}/results/DADA2.biom \
    ${path} \
    ${study_id}

qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path ${path}/DADA2.biom \
  --input-format BIOMV210Format \
  --output-path ${path}/feature-table.qza

qiime feature-table rarefy \
  --i-table ${path}/feature-table.qza \
  --p-sampling-depth 1000 \
  --o-rarefied-table ${path}/rarefied-table.qza

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ${path}/rep-seqs.fasta \
  --output-path ${path}/rep-seqs.qza

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${path}/rep-seqs.qza \
  --o-alignment ${path}/aligned-rep-seqs.qza \
  --o-masked-alignment ${path}/masked-aligned-rep-seqs.qza \
  --o-tree ${path}/unrooted-tree.qza \
  --o-rooted-tree ${path}/rooted-tree.qza

qiime diversity alpha \
  --i-table ${path}/rarefied-table.qza \
  --p-metric shannon \
  --o-alpha-diversity ${path}/shannon.qza

qiime diversity alpha \
  --i-table ${path}/rarefied-table.qza \
  --p-metric observed_features \
  --o-alpha-diversity ${path}/observed-features.qza

qiime diversity alpha \
  --i-table ${path}/rarefied-table.qza \
  --p-metric chao1 \
  --o-alpha-diversity ${path}/chao1.qza

qiime diversity alpha \
  --i-table ${path}/rarefied-table.qza \
  --p-metric simpson \
  --o-alpha-diversity ${path}/simpson.qza

qiime diversity alpha-phylogenetic \
  --i-table ${path}/rarefied-table.qza \
  --i-phylogeny ${path}/rooted-tree.qza \
  --p-metric faith_pd \
  --o-alpha-diversity ${path}/faith-pd.qza

#!/bin/bash
#SBATCH --job-name=alpha_div
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p cu
#SBATCH --mem=10G
#SBATCH -o output/alpha_div_%a.out
#SBATCH -e output/alpha_div_%a.err
#SBATCH --array=1
#SBATCH --exclude=cu01

source aactivate QIIME2_ENV

study_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../data/project_list.txt | tr -d '\r\n')

### mkdir alpha file
mkdir /beegfs/dongbiao/oral/${study_id}/results/alpha
path=/beegfs/dongbiao/oral/${study_id}/results/alpha

### getting oral metadata and biom file for qiiem2 input
filter_metadata_for_QIIME2.py \
    ../../data/projects/${study_id}/metadata_v2.csv \
    /beegfs/dongbiao/oral/${study_id}/results/DADA2.biom \
    ${path}

### import biom to QIIME2 and metadata
qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path ${path}/DADA2.biom \
  --input-format BIOMV210Format \
  --output-path ${path}/feature-table.qza
qiime feature-table rarefy \
  --i-table ${path}/feature-table.qza \
  --p-sampling-depth 1000 \
  --o-rarefied-table ${path}/rarefied-table.qza

### build tree
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ${path}/rep-seqs.fasta \
  --output-path ${path}/rep-seqs.qza
  
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${path}rep-seqs.qza \
  --o-alignment ${path}/aligned-rep-seqs.qza \
  --o-masked-alignment ${path}/masked-aligned-rep-seqs.qza \
  --o-tree ${path}/unrooted-tree.qza \
  --o-rooted-tree ${path}/rooted-tree.qza

### caculate alpha diversity
# Shannon
qiime diversity alpha \
  --i-table ${path}/rarefied-table.qza \
  --p-metric shannon \
  --o-alpha-diversity ${path}/shannon.qza

# Observed Features
qiime diversity alpha \
  --i-table ${path}/rarefied-table.qza \
  --p-metric observed_features \
  --o-alpha-diversity ${path}/observed-features.qza

# Simpson 
qiime diversity alpha \
  --i-table ${path}/rarefied-table.qza \
  --p-metric simpson \
  --o-alpha-diversity ${path}/simpson.qza

# Faith's PD
qiime diversity alpha-phylogenetic \
  --i-table ${path}/rarefied-table.qza \
  --i-phylogeny ${path}/rooted-tree.qza \
  --p-metric faith_pd \
  --o-alpha-diversity ${path}/faith-pd.qza


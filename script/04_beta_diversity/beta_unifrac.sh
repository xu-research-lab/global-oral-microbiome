#!/bin/bash
#SBATCH --job-name=beta_unifrac
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH --mem=250G
#SBATCH -o output/beta_unifrac.out
#SBATCH -e output/beta_unifrac.err
#SBATCH --exclude=cu01,cu04

source /home/hexincheng/miniconda3/bin/activate /home/hexincheng/miniconda3/envs/rachis-qiime2-2026.4
export OMP_NUM_THREADS=28

input_dir=/beegfs/dongbiao/oral/beta
work_dir=${input_dir}/unifrac
reference_tree_nwk=/beegfs/dongbiao/greengene2/2024.09.phylogeny_backbone.nwk
biom_qza=${input_dir}/rarefied-table.qza
reference_tree_qza=${work_dir}/reference-tree.qza
exported_table=${work_dir}/exported_table
export_tree=${work_dir}/export_tree

mkdir -p ${work_dir} ${exported_table} ${export_tree}

qiime tools import \
  --input-path ${reference_tree_nwk} \
  --output-path ${reference_tree_qza} \
  --type 'Phylogeny[Rooted]'

qiime phylogeny filter-tree \
  --i-tree ${reference_tree_qza} \
  --i-table ${biom_qza} \
  --o-filtered-tree ${work_dir}/rooted-tree.qza

qiime tools export \
  --input-path ${work_dir}/rooted-tree.qza \
  --output-path ${export_tree}

qiime tools export \
  --input-path ${biom_qza} \
  --output-path ${exported_table}

ssu -i ${exported_table}/feature-table.biom -o ${work_dir}/unifrac_unweight.dm -m unweighted -t ${export_tree}/tree.nwk
ssu -i ${exported_table}/feature-table.biom -o ${work_dir}/unifrac_weight.dm -m weighted_normalized -t ${export_tree}/tree.nwk

#ssu -i ${exported_table}/feature-table.biom -o ${work_dir}/unifrac_unweight_pcoa.hdf5 -m unweighted -t ${export_tree}/tree.nwk --pcoa 3 --format hdf5_nodist
#ssu -i ${exported_table}/feature-table.biom -o ${work_dir}/unifrac_weight_pcoa.hdf5 -m weighted_normalized -t ${export_tree}/tree.nwk --pcoa 3 --format hdf5_nodist

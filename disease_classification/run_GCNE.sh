#!/bin/bash
#SBATCH --job-name=GCNE
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p gpu
#SBATCH --mem=20G
#SBATCH -o GCNE_part1.out
#SBATCH -e GCNE_part1.err

module load miniconda/4.9.2 
source activate jupyter_notebook

model_file=/home/dongbiao/GCN
file_path=/home/dongbiao/GCN/data/Gastritis_16S
PYTHON=/home/dongbiao/.conda/envs/jupyter_notebook/bin/python


### part 1
while IFS=',' read -r site disease study; do
    if [[ "$site" == "site" ]]; then
            continue
    fi
      
    CUDA_VISIBLE_DEVICES=4 $PYTHON ${model_file}/GCN_model.py \
        --train_table data/lso_train_${site}_${disease}_${study}.biom \
        --test_table data/lso_test_${site}_${disease}_${study}.biom \
        --metadata_filename meta.tsv \
        --phylogeny_file_path data/phylogeny_normalize_loso_${site}_${disease}.nwk \
        --fold ${site}_${disease}_${study} --output_dir ./result --batch_size 32 \
        --learning_rate 0.001 --num_layers 4 --ffn True \
        --norm_type batch_norm --residual True \
        --hidden_channels 100 --dropout_rate 0.1 \
        --weight_decay 0.00001 
done < run_info.csv

### part 2
# while IFS=',' read -r site fold study; do
#     if [[ "$site" == "site" ]]; then
#             continue
#     fi
    
#     CUDA_VISIBLE_DEVICES=7 $PYTHON ${software_path}/deepphylo_classification.py \
#       --train_table data/five_fold_train_${site}_${study}_${fold}.biom \
#       --test_table data/five_fold_test_${site}_${study}_${fold}.biom \
#       --metadata_filename meta_oral_gut.tsv --pca_file_path data \
#       --epochs 50 \
#       -hs 32 64 80 -kec 7 9 -kep 2 4 -l 1e-4 1e-3 -bs 32 -act relu tanh \
#       --fold ${site}_${study}_${fold} --output_dir ./result
# done < run_info_oral_gut.csv

### part 3
# while IFS=',' read -r site disease continent test_study split_type; do
#     if [[ "$site" == "site" ]]; then
#             continue
#     fi
    
#     if [[ "$split_type" == "wc" ]]; then
#     train_biom=train_${split_type}_${site}_${disease}_${continent}_${test_study}
#     test_biom=train_${split_type}_${site}_${disease}_${continent}_${test_study}
#     fold=${split_type}_${site}_${disease}_${continent}_${test_study}
#     else
#     train_biom=train_${split_type}_${site}_${disease}_${continent}
#     test_biom=train_${split_type}_${site}_${disease}_${continent}
#     fold=${split_type}_${site}_${disease}_${continent}
#     fi
    
#     CUDA_VISIBLE_DEVICES=7 $PYTHON ${software_path}/deepphylo_classification.py \
#       --train_table data/${train_biom}.biom \
#       --test_table data/${test_biom}.biom \
#       --metadata_filename meta.tsv --pca_file_path ${pca_file_path} \
#       --epochs 50 \
#       -hs 32 64 80 -kec 7 9 -kep 2 4 -l 1e-4 1e-3 -bs 32 -act relu tanh \
#       --fold ${fold} --output_dir ./result
# done < run_info_continent.csv

### part 4
# while IFS=',' read -r site train_disease test_disease n_train n_test n_train_case n_test_case split_type; do
#     if [[ "$site" == "site" ]]; then
#             continue
#     fi
    
#     CUDA_VISIBLE_DEVICES=7 $PYTHON ${software_path}/deepphylo_classification.py \
#       --train_table data/train_cd_${site}_train_${train_disease}_test_${test_disease}.biom \
#       --test_table data/test_cd_${site}_train_${train_disease}_test_${test_disease}.biom \
#       --metadata_filename meta.tsv --pca_file_path ${pca_file_path} \
#       --epochs 50 \
#       -hs 32 64 80 -kec 7 9 -kep 2 4 -l 1e-4 1e-3 -bs 32 -act relu tanh \
#       --fold cd_${site}_train_${train_disease}_test_${test_disease} --output_dir ./result
# done < run_info_cross_disease.csv

### part 5
# while IFS=',' read -r site disease n_studies repeat subset test_study split_type; do
#     if [[ "$site" == "site" ]]; then
#             continue
#     fi

#     CUDA_VISIBLE_DEVICES=6 $PYTHON ${software_path}/deepphylo_classification.py \
#       --train_table data/train_lc_${site}_${disease}_k${n_studies}_r${repeat}_test_${test_study}.biom \
#       --test_table data/test_lc_${site}_${disease}_k${n_studies}_r${repeat}_test_${test_study}.biom \
#       --metadata_filename meta.tsv --pca_file_path ${pca_file_path} \
#       --epochs 50 \
#       -hs 32 64 80 -kec 7 9 -kep 2 4 -l 1e-4 1e-3 -bs 32 -act relu tanh \
#       --fold lc_${site}_${disease}_k${n_studies}_r${repeat}_test_${test_study} --output_dir ./result
# done < run_info_lc.csv

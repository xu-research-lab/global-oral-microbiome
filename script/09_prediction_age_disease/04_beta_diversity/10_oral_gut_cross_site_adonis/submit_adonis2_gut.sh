#!/bin/bash
#SBATCH --job-name=adonis2_gut
#SBATCH --output=adonis2_gut_%j.out
#SBATCH --error=adonis2_gut_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=cu09

cd /beegfs/huangxiaochang/global_oral/oral_gut_0127/gut

# 加载R环境（根据您的集群环境调整）
source /home/software/anaconda3/etc/profile.d/conda.sh
conda activate renv

# 运行分析
Rscript run_adonis2_gut.R

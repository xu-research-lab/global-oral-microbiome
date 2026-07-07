#!/bin/bash
#SBATCH --job-name=adonis2_oral
#SBATCH --output=adonis2_oral_%j.out
#SBATCH --error=adonis2_oral_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

cd /beegfs/huangxiaochang/global_oral/oral_gut_0127/oral

# 加载R环境（根据您的集群环境调整）
module load R/4.2.0 2>/dev/null || echo "Using system R"

# 运行分析
Rscript run_adonis2_oral.R

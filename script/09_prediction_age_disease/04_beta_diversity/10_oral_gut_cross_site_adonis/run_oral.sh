#!/bin/bash
#SBATCH --job-name=oral_adonis
#SBATCH --ntasks=1              # 1个任务
#SBATCH --cpus-per-task=28      # 申请28个CPU
#SBATCH -p cu
#SBATCH --mem=240G              # 内存240G
#SBATCH -o oral_adonis.out
#SBATCH -e oral_adonis.err
#SBATCH --exclude=cu01,cu09     # 排除故障节点

# --- 关键修复：在 R 启动前强制限制底层线程数 ---
# 这将确保每个 R 子进程严格只使用 1 个物理核，防止负载爆炸
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# 导入环境变量
export PATH="/home/hexincheng/miniconda3/envs/r_env/bin:$PATH"

cd /home/hexincheng/task/huang_teacher/gut_adonis_analysis_0115/oral_gut_0127/oral

# 运行 R 脚本
Rscript run_adonis_oral.R
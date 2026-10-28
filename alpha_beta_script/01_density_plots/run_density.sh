#!/bin/bash
#SBATCH --job-name=density
#SBATCH --partition=cu
#SBATCH --output=/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/output/density_%a.out
#SBATCH --error=/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/output/density_%a.err
#SBATCH --array=1-7%7
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --exclude=cu01

source /home/hexincheng/miniconda3/bin/activate /home/hexincheng/miniconda3/envs/r_env

mkdir -p /home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/beta_pcoa_out_20260909/pcoa_plots/density_logs
mkdir -p /home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/beta_pcoa_out_20260909/pcoa_plots/density

Rscript /home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/04-plot-density-contour-pcoa.R
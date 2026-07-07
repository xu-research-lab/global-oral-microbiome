#!/bin/bash
#SBATCH --job-name=genus-relative-collapse-260610
#SBATCH --partition=cu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/yanglun/YSD/global_oral/global_oral_v2/scripts/continent-analysis/scripts/logs/genus-relative-collapse-260610_%j.out
#SBATCH --error=/home/yanglun/YSD/global_oral/global_oral_v2/scripts/continent-analysis/scripts/logs/genus-relative-collapse-260610_%j.err
#SBATCH --exclude=cu01

set -euo pipefail

# ============================================================
# 路径定义
# ============================================================
PROJ_DIR="/home/yanglun/YSD/global_oral/global_oral_v2/scripts/continent-analysis"
LOG_DIR="${PROJ_DIR}/scripts/logs"
OUT_DIR="${PROJ_DIR}/input"
EXPORT_DIR="${OUT_DIR}/export-genus"

# 输入文件
BIOM="/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/biom_oral_qc.biom"
TAX="${PROJ_DIR}/input/taxonomy.tsv"

# 中间产物 & 最终产物
FEATURE_TABLE_QZA="${OUT_DIR}/feature-table.qza"
TAXONOMY_QZA="${OUT_DIR}/taxonomy.qza"
GENUS_TABLE_QZA="${OUT_DIR}/genus-table.qza"
GENUS_REL_QZA="${OUT_DIR}/genus-relative-table.qza"
GENUS_REL_TSV="${OUT_DIR}/genus_relative_abundance.tsv"

# 折叠到的分类层级：1=Domain 2=Phylum 3=Class 4=Order 5=Family 6=Genus 7=Species
LEVEL=6

# ============================================================
# 目录预创建
# ============================================================
mkdir -p "${LOG_DIR}" "${OUT_DIR}" "${EXPORT_DIR}"

# ============================================================
# 输入校验
# ============================================================
for f in "${BIOM}" "${TAX}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Input file not found: ${f}" >&2
        exit 1
    fi
done

# ============================================================
# 线程设置
# ============================================================
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# ============================================================
# 激活环境
# ============================================================
set +u
source /home/yanglun/miniconda3/bin/activate qiime2-2024.2
set -u

echo "=========================================="
echo "JobID:        ${SLURM_JOB_ID}"
echo "Node:         $(hostname)"
echo "Start:        $(date)"
echo "QIIME2:       $(qiime --version 2>&1 | head -1)"
echo "Input biom:   ${BIOM}"
echo "Input tax:    ${TAX}"
echo "Output dir:   ${OUT_DIR}"
echo "Collapse lvl: ${LEVEL} (genus)"
echo "=========================================="

# ============================================================
# Step 1: biom -> qza (FeatureTable[Frequency])
# ============================================================
echo ""
echo "[Step 1/5] Import biom -> feature-table.qza"
date
if [[ -f "${FEATURE_TABLE_QZA}" ]]; then
    echo "  Skip: ${FEATURE_TABLE_QZA} already exists."
else
    qiime tools import \
      --input-path "${BIOM}" \
      --type 'FeatureTable[Frequency]' \
      --input-format BIOMV210Format \
      --output-path "${FEATURE_TABLE_QZA}"
fi

# ============================================================
# Step 2: taxonomy.tsv -> qza (FeatureData[Taxonomy])
# ============================================================
echo ""
echo "[Step 2/5] Import taxonomy.tsv -> taxonomy.qza"
date
if [[ -f "${TAXONOMY_QZA}" ]]; then
    echo "  Skip: ${TAXONOMY_QZA} already exists."
else
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format TSVTaxonomyFormat \
      --input-path "${TAX}" \
      --output-path "${TAXONOMY_QZA}"
fi

# ============================================================
# Step 3: taxa collapse to genus level
# ============================================================
echo ""
echo "[Step 3/5] taxa collapse -> genus-table.qza  (level=${LEVEL})"
date
qiime taxa collapse \
  --i-table "${FEATURE_TABLE_QZA}" \
  --i-taxonomy "${TAXONOMY_QZA}" \
  --p-level ${LEVEL} \
  --o-collapsed-table "${GENUS_TABLE_QZA}"

# ============================================================
# Step 4: 转相对丰度
# ============================================================
echo ""
echo "[Step 4/5] relative-frequency -> genus-relative-table.qza"
date
qiime feature-table relative-frequency \
  --i-table "${GENUS_TABLE_QZA}" \
  --o-relative-frequency-table "${GENUS_REL_QZA}"

# ============================================================
# Step 5: 导出为 tsv
# ============================================================
echo ""
echo "[Step 5/5] export -> genus_relative_abundance.tsv"
date
qiime tools export \
  --input-path "${GENUS_REL_QZA}" \
  --output-path "${EXPORT_DIR}"

biom convert \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "${GENUS_REL_TSV}" \
  --to-tsv

# ============================================================
# 结果摘要
# ============================================================
echo ""
echo "=========================================="
echo "DONE"
echo "Final TSV:    ${GENUS_REL_TSV}"
echo "Size:         $(du -h "${GENUS_REL_TSV}" | cut -f1)"
echo "Genus count:  $(($(wc -l < "${GENUS_REL_TSV}") - 2))   (excluding 2 header lines)"
echo "Top 3 lines:"
head -3 "${GENUS_REL_TSV}" | cut -f1-5
echo "End:          $(date)"
echo "=========================================="

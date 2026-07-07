#!/usr/bin/env Rscript

# ============================================================
# MMUPHin 批次效应校正脚本
# 对属水平物种丰度表进行批次效应校正
# ============================================================

# 加载必要的包
suppressPackageStartupMessages({
  library(MMUPHin)
  library(dplyr)
})

cat("=== MMUPHin 批次效应校正 ===\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================
# 1. 加载数据
# ============================================================
cat("[1/5] 加载数据...\n")

# 加载丰度表
load("genus_data_relabund.RData")
cat("  - 丰度表维度:", nrow(genus_data), "行(菌属) x", ncol(genus_data), "列(样本)\n")

# 读取元数据
metadata <- read.delim("qiime2_metadata.tsv", stringsAsFactors = FALSE, check.names = FALSE)
cat("  - 元数据维度:", nrow(metadata), "行 x", ncol(metadata), "列\n")

# ============================================================
# 2. 数据预处理
# ============================================================
cat("\n[2/5] 数据预处理...\n")

# 提取样本ID（丰度表的列名）
sample_ids <- colnames(genus_data)

# 检查元数据中的样本ID列名
sample_id_col <- "#SampleID"
if (!(sample_id_col %in% colnames(metadata))) {
  if ("SampleID" %in% colnames(metadata)) {
    sample_id_col <- "SampleID"
  } else if ("sample_id" %in% colnames(metadata)) {
    sample_id_col <- "sample_id"
  }
}

# 确保元数据中的样本ID与丰度表一致
metadata_matched <- metadata[metadata[[sample_id_col]] %in% sample_ids, ]
cat("  - 匹配到元数据的样本数:", nrow(metadata_matched), "\n")

# 按丰度表列名顺序重新排列元数据
metadata_matched <- metadata_matched[match(sample_ids, metadata_matched[[sample_id_col]]), ]

# 确保批次变量是因子，并去除未使用的水平（关键！）
metadata_matched$BioProject <- droplevels(as.factor(metadata_matched$BioProject))

# ============================================================
# 3. 准备 MMUPHin 输入
# ============================================================
cat("\n[3/5] 准备 MMUPHin 输入...\n")

# MMUPHin 要求：特征-by-样本矩阵（特征在行，样本在列）
# 原始数据正是这个格式，无需转置
abundance_matrix <- as.matrix(genus_data)
# 修复：确保列名是字符类型（不是array），以便与rownames匹配
if (is.array(colnames(abundance_matrix))) {
  colnames(abundance_matrix) <- as.character(colnames(abundance_matrix))
}
cat("  - 丰度表矩阵维度:", nrow(abundance_matrix), "行(菌属) x", ncol(abundance_matrix), "列(样本)\n")

# 检查样本顺序是否一致
if (!all(colnames(abundance_matrix) == metadata_matched[[sample_id_col]])) {
  cat("  - 警告：样本顺序不一致，进行调整...\n")
  # 确保元数据顺序与丰度表一致
  metadata_matched <- metadata_matched[match(colnames(abundance_matrix), metadata_matched[[sample_id_col]]), ]
}

cat("  - 批次数量 (BioProject):", length(levels(metadata_matched$BioProject)), "\n")
cat("  - 批次分布 (前10个):\n")
print(head(table(metadata_matched$BioProject), 10))

# 准备需要保护的协变量
# 注意：所有协变量都与 Bioproject 高度混杂，MMUPHin 不允许同时使用
covariates_to_protect <- NULL  # 不使用协变量保护

# 检查协变量是否存在
available_covariates <- covariates_to_protect[covariates_to_protect %in% colnames(metadata_matched)]
missing_covariates <- covariates_to_protect[!(covariates_to_protect %in% colnames(metadata_matched))]

if (length(missing_covariates) > 0) {
  cat("  - 警告: 以下协变量不存在:", paste(missing_covariates, collapse = ", "), "\n")
}
cat("  - 将保护的协变量:", paste(available_covariates, collapse = ", "), "\n")

# 将协变量转换为因子（如果是字符型）
for (cov in available_covariates) {
  if (is.character(metadata_matched[[cov]])) {
    metadata_matched[[cov]] <- droplevels(as.factor(metadata_matched[[cov]]))
  }
}

# 关键：MMUPHin 要求 data 的行名必须与 feature_abd 的列名一致
# 将元数据的行名设为样本ID
rownames(metadata_matched) <- metadata_matched[[sample_id_col]]
cat("  - 元数据行名已设为样本ID (前5个):", paste(head(rownames(metadata_matched), 5), collapse = ", "), "\n")

# 检查缺失值
na_summary <- sapply(metadata_matched[, c("BioProject", available_covariates), drop = FALSE], function(x) sum(is.na(x)))
if (sum(na_summary) > 0) {
  cat("  - 批次和协变量缺失值统计:\n")
  print(na_summary[na_summary > 0])
  cat("  - 注意：MMUPHin会自动处理包含NA的样本\n")
}

# ============================================================
# 4. 批次效应校正 (使用 MMUPHin)
# ============================================================
cat("\n[4/5] 执行 MMUPHin 批次效应校正...\n")

# 根据 MMUPHin 文档：
# - feature_abd: feature-by-sample matrix (特征在行，样本在列)
# - batch: 批次变量名（字符串）
# - covariates: 协变量名向量（字符串）
# - data: 包含批次和协变量的数据框

cat("  - 调用 adjust_batch 函数...\n")
cat("    feature_abd 维度:", dim(abundance_matrix), "\n")
cat("    feature_abd 列名(样本)前5个:", paste(head(colnames(abundance_matrix), 5), collapse = ", "), "\n")
cat("    data 行名前5个:", paste(head(rownames(metadata_matched), 5), collapse = ", "), "\n")
cat("    样本名是否一致:", all(colnames(abundance_matrix) == rownames(metadata_matched)), "\n")
cat("    batch 变量: BioProject\n")
cat("    covariates:", paste(available_covariates, collapse = ", "), "\n")
cat("    data 维度:", dim(metadata_matched), "\n")

tryCatch({
  adjusted_result <- adjust_batch(
    feature_abd = abundance_matrix,
    batch = "BioProject",                    # 变量名字符串
    covariates = available_covariates,       # 变量名字符串向量
    data = metadata_matched,                 # 数据框
    control = list(
      zero_inflation = TRUE,                 # 使用零膨胀模型
      pseudo_count = NULL,                   # 自动设置伪计数
      diagnostic_plot = "mmuphin_diagnostic.pdf",
      conv = 1e-4,
      maxit = 1000,
      verbose = TRUE
    )
  )

  cat("  - 批次效应校正成功!\n")

  # 提取校正后的丰度表
  adjusted_abundance <- adjusted_result$feature_abd_adj

}, error = function(e) {
  cat("\n  ERROR: MMUPHin 校正失败\n")
  cat("  错误信息:", conditionMessage(e), "\n")

  # 打印更详细的错误信息
  cat("\n  --- 调试信息 ---\n")
  cat("  abundance_matrix 类型:", class(abundance_matrix), "\n")
  cat("  abundance_matrix 维度:", paste(dim(abundance_matrix), collapse = " x "), "\n")
  cat("  abundance_matrix 前5行3列:\n")
  print(abundance_matrix[1:min(5, nrow(abundance_matrix)), 1:min(3, ncol(abundance_matrix))])

  cat("\n  metadata_matched 行数:", nrow(metadata_matched), "\n")
  cat("  metadata_matched 列数:", ncol(metadata_matched), "\n")
  cat("  BioProject 类型:", class(metadata_matched$BioProject), "\n")
  cat("  BioProject 水平数:", length(levels(metadata_matched$BioProject)), "\n")

  cat("\n  各协变量类型:\n")
  for (cov in available_covariates) {
    cat("    ", cov, ":", class(metadata_matched[[cov]]),
        " NA数:", sum(is.na(metadata_matched[[cov]])), "\n")
  }

  # 重新抛出错误以停止执行
  stop("MMUPHin校正失败，请检查上述调试信息")
})

# ============================================================
# 5. 后处理和保存
# ============================================================
cat("\n[5/5] 后处理和保存结果...\n")

# 确保校正后的丰度非负
adjusted_abundance[adjusted_abundance < 0] <- 0

# 保持原始行名和列名
rownames(adjusted_abundance) <- rownames(genus_data)
colnames(adjusted_abundance) <- colnames(genus_data)

cat("  - 校正后丰度范围:", round(min(adjusted_abundance), 6), "-", round(max(adjusted_abundance), 6), "\n")
cat("  - 零值比例:", round(sum(adjusted_abundance == 0) / length(adjusted_abundance) * 100, 2), "%\n")

# 转换回数据框
genus_data_corrected <- as.data.frame(adjusted_abundance)

# 保存校正后的丰度表
save(genus_data_corrected, file = "genus_data_corrected_mmuphin.RData")
write.table(genus_data_corrected,
            file = "genus_data_corrected_mmuphin.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

cat("  - 结果已保存到:\n")
cat("    - genus_data_corrected_mmuphin.RData\n")
cat("    - genus_data_corrected_mmuphin.tsv\n")
cat("    - mmuphin_diagnostic.pdf (诊断图)\n")

# ============================================================
# 6. 批次效应评估
# ============================================================
cat("\n=== 批次效应评估 ===\n")

if (requireNamespace("vegan", quietly = TRUE)) {
  library(vegan)

  calc_r2 <- function(abd_matrix, batch_vec) {
    abd_matrix <- abd_matrix + 1e-6
    hel <- decostand(t(abd_matrix), method = "hellinger")
    dist_mat <- vegdist(hel, method = "bray")
    adonis_result <- adonis2(dist_mat ~ batch_vec, permutations = 99)
    return(adonis_result$R2[1])
  }

  cat("计算 PERMANOVA R² (基于Bray-Curtis距离)...\n")

  set.seed(42)
  sample_idx <- sample(1:ncol(adjusted_abundance), min(1000, ncol(adjusted_abundance)))

  r2_before <- calc_r2(as.matrix(genus_data)[, sample_idx], metadata_matched$BioProject[sample_idx])
  r2_after <- calc_r2(adjusted_abundance[, sample_idx], metadata_matched$BioProject[sample_idx])

  cat("  - 校正前批次解释方差比例 (R²):", round(r2_before, 4), "\n")
  cat("  - 校正后批次解释方差比例 (R²):", round(r2_after, 4), "\n")
  if (r2_before > 0) {
    cat("  - R²降低比例:", round((1 - r2_after/r2_before) * 100, 2), "%\n")
  }
}

cat("\n=== 完成 ===\n")
cat("结束时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#!/usr/bin/env Rscript
# Gut Adonis2 Analysis Script
# Formula: dist_matrix ~ Sex + Age_Group + Seq_Region + Geography + Healthy_Status
# Healthy_Status: YES (健康) / NO (不健康)
# by = "margin", strata = metadata$Bioproject

# --- 1. 加载库 ---
library(vegan)
library(data.table)
library(ggplot2)

# --- 2. 安全并行设置 ---
n_cores <- 8
cat("Running with Safe Cores:", n_cores, "\n")
Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

work_dir <- "/beegfs/huangxiaochang/global_oral/oral_gut_0127/gut"
setwd(work_dir)

# --- 3. 数据清洗函数 ---
clean_metadata <- function(meta) {
  for (col in names(meta)) {
    if (is.character(meta[[col]]) || is.factor(meta[[col]])) {
      meta[[col]] <- as.character(meta[[col]])
      # 空字符串或只有空格转为NA
      meta[[col]][trimws(meta[[col]]) == ""] <- NA
      # 明确的不明确值转为NA
      unclear_values <- c("Others_region", "Unknown", "unknown", "NA", "N/A", "n/a",
                          "uncategorized", "Uncategorized", "Others", "others")
      meta[[col]][meta[[col]] %in% unclear_values] <- NA
    }
  }
  return(meta)
}

# --- 4. 处理Healthy_Status函数 ---
process_healthy_status <- function(meta) {
  # 如果已经有Healthy_Status列，直接使用
  if ("Healthy_Status" %in% names(meta)) {
    cat("  Using existing Healthy_Status column\n")
    return(meta)
  }

  # 如果有Disease_study列，从中推断
  if ("Disease_study" %in% names(meta)) {
    cat("  Creating Healthy_Status from Disease_study\n")
    meta$Healthy_Status <- NA

    # 健康状态："no", "No"
    healthy_vals <- c("no", "No")
    meta$Healthy_Status[meta$Disease_study %in% healthy_vals] <- "YES"

    # 疾病状态：非"no"/"No"且非NA的其他值
    # 但排除不明确的状态
    unclear_disease <- c(NA, "", " ", "NA", "N/A", "n/a")
    is_disease <- !(meta$Disease_study %in% c(healthy_vals, unclear_disease)) &
                   !is.na(meta$Disease_study)
    meta$Healthy_Status[is_disease] <- "NO"

    cat("    YES (Healthy):", sum(meta$Healthy_Status == "YES", na.rm = TRUE), "\n")
    cat("    NO (Disease):", sum(meta$Healthy_Status == "NO", na.rm = TRUE), "\n")
    cat("    NA (Unknown):", sum(is.na(meta$Healthy_Status)), "\n")

    return(meta)
  }

  cat("  Warning: No Healthy_Status or Disease_study column found\n")
  return(meta)
}

# --- 5. 读取Metadata ---
cat("Reading Metadata...\n")
meta <- read.csv("metadata_gut_all_0127.csv", check.names = FALSE, stringsAsFactors = FALSE)

# 设置行名
if ("sample-id" %in% names(meta)) {
  rownames(meta) <- meta$"sample-id"
} else if ("SampleID" %in% names(meta)) {
  rownames(meta) <- meta$SampleID
} else {
  rownames(meta) <- meta[,1]
}

cat("Original samples:", nrow(meta), "\n")

# 清洗数据
meta <- clean_metadata(meta)

# 处理Healthy_Status
meta <- process_healthy_status(meta)

# 处理Age列（如果存在）
if ("Age" %in% names(meta)) {
  meta$Age <- suppressWarnings(as.numeric(as.character(meta$Age)))
  meta$Age[meta$Age <= 0] <- NA
}

# --- 6. 确定可用列 ---
# 目标列：BioProject(用于strata), Sex, Age_Group, Seq_Region, Geography, Healthy_Status

# Bioproject列（用于strata）
if ("Bioproject" %in% names(meta)) {
  bioproject_col <- "Bioproject"
} else if ("BioProject" %in% names(meta)) {
  bioproject_col <- "BioProject"
} else {
  stop("Error: No Bioproject/BioProject column found for strata")
}

# 构建目标列列表
target_cols <- c(bioproject_col, "Sex", "Age_Group", "Seq_Region", "Geography", "Healthy_Status")

cat("Target columns:", paste(target_cols, collapse = ", "), "\n")

# 检查哪些列实际存在
available_cols <- target_cols[target_cols %in% names(meta)]
missing_cols <- target_cols[!target_cols %in% names(meta)]

if (length(missing_cols) > 0) {
  cat("Warning: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
}

# --- 7. 选择列并清洗 ---
meta_clean <- meta[, available_cols, drop = FALSE]
meta_clean <- na.omit(meta_clean)
# 去除因子中的空水平
for (col in names(meta_clean)) {
  if (is.factor(meta_clean[[col]])) {
    meta_clean[[col]] <- droplevels(meta_clean[[col]])
  }
}

cat("After NA removal:", nrow(meta_clean), "samples (dropped:", nrow(meta) - nrow(meta_clean), ")\n")

if (nrow(meta_clean) == 0) {
  stop("Error: No samples left after cleaning!")
}

# 显示各因子的水平数
cat("\nFactor levels:\n")
for (col in available_cols) {
  if (col %in% names(meta_clean)) {
    n_levels <- length(unique(meta_clean[[col]]))
    cat("  ", col, ":", n_levels, "levels")
    if (col == "Healthy_Status") {
      cat(" (", paste(unique(meta_clean[[col]]), collapse = ", "), ")")
    }
    cat("\n")
    if (n_levels < 2) {
      cat("    Warning:", col, "has only 1 level\n")
    }
  }
}

# --- 8. 读取距离矩阵 ---
cat("\nReading Distance Matrix...\n")
dist_df <- fread("unifrac_weight.dm", header = TRUE, data.table = FALSE)
rownames(dist_df) <- dist_df[,1]
dist_df <- dist_df[,-1]
cat("Distance matrix dimensions:", dim(dist_df), "\n")

# --- 9. 样本对齐 ---
common_samples <- intersect(rownames(meta_clean), rownames(dist_df))
cat("Matched samples:", length(common_samples), "\n")

if (length(common_samples) < 10) {
  stop("Error: Too few matched samples!")
}

meta_final <- meta_clean[common_samples, , drop = FALSE]
dist_final <- as.dist(dist_df[common_samples, common_samples])

rm(dist_df, meta, meta_clean); gc()

# 确保Bioproject列存在（用于strata）
if (bioproject_col == "BioProject" && !"Bioproject" %in% names(meta_final)) {
  meta_final$Bioproject <- meta_final$BioProject
}

# --- 10. 运行 Adonis2 ---
cat("\n========================================\n")
cat("Running Adonis2 with formula:\n")
cat("  dist ~ Sex + Age_Group + Seq_Region + Geography + Healthy_Status\n")
cat("by = 'margin', strata = Bioproject\n")
cat("========================================\n\n")

tryCatch({
  ado_res <- adonis2(dist_final ~ Sex + Age_Group + Seq_Region + Geography + Healthy_Status,
                     data = meta_final,
                     permutations = 999,
                     by = "margin",
                     strata = meta_final$Bioproject,
                     parallel = n_cores)

  print(ado_res)
  write.csv(as.data.frame(ado_res), "adonis2_result_gut_0327.csv")
  cat("\nResults saved to: adonis2_result_gut_0327.csv\n")

}, error = function(e) {
  cat("Error in adonis2:", conditionMessage(e), "\n")

  # 尝试简化模型
  cat("\nTrying simplified model without Geography...\n")
  tryCatch({
    ado_res <- adonis2(dist_final ~ Sex + Age_Group + Seq_Region + Healthy_Status,
                       data = meta_final,
                       permutations = 999,
                       by = "margin",
                       strata = meta_final$Bioproject,
                       parallel = n_cores)
    print(ado_res)
    write.csv(as.data.frame(ado_res), "adonis2_result_gut_0327_no_geo.csv")
    cat("\nResults saved to: adonis2_result_gut_0327_no_geo.csv\n")
  }, error = function(e2) {
    cat("Simplified model also failed:", conditionMessage(e2), "\n")
  })
})

cat("\nDone!\n")

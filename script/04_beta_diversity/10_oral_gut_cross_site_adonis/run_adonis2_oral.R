#!/usr/bin/env Rscript
# Oral Adonis2 Analysis Script
# Formula: dist_matrix ~ Sex + Age_Group + Seq_Region + Geography + Healthy_Status
# Healthy_Status: YES (健康) / NO (不健康)
# by = "margin", strata = metadata$Bioproject

# --- 1. 加载库 ---
library(vegan)
library(ggplot2)

# --- 2. 安全并行设置 ---
n_cores <- 8
cat("Running with Safe Cores:", n_cores, "\n")
Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

work_dir <- "/beegfs/huangxiaochang/global_oral/oral_gut_0127/oral"
setwd(work_dir)

# --- 3. 加载距离矩阵 ---
cat("Loading unifrac.RData...\n")
loaded_objs <- load("unifrac.RData")
dist_big <- as.matrix(get(loaded_objs[1]))
cat("Distance matrix dimensions:", dim(dist_big), "\n")

# --- 4. 定义要分析的文件和目标列 ---
sites_files <- c("metadata_Saliva_0127.csv", "metadata_Oral_Mucosa_0127.csv",
                 "metadata_Tongue_0127.csv", "metadata_Tooth_0127.csv")

# 目标列：BioProject(用于strata), Sex, Age_Group, Seq_Region, Geography, Healthy_Status
target_cols <- c("BioProject", "Sex", "Age_Group", "Seq_Region", "Geography", "Healthy_Status")

# --- 5. 数据清洗函数 ---
clean_metadata <- function(meta) {
  for (col in names(meta)) {
    if (is.character(meta[[col]]) || is.factor(meta[[col]])) {
      meta[[col]] <- as.character(meta[[col]])
      # 空字符串或只有空格转为NA
      meta[[col]][trimws(meta[[col]]) == ""] <- NA
      # 明确的不明确值转为NA
      unclear_values <- c("Others_region", "Unknown", "unknown", "NA", "N/A", "n/a",
                          "uncategorized", "Uncategorized", "Others", "others",
                          "Uncategorized")
      meta[[col]][meta[[col]] %in% unclear_values] <- NA
    }
  }
  return(meta)
}

# --- 6. 处理Healthy_Status函数 ---
process_healthy_status <- function(meta) {
  if ("Healthy_Status" %in% names(meta)) {
    # 创建新的Healthy_Status_YN列
    meta$Healthy_Status_YN <- NA

    # 健康样本标记为 YES
    meta$Healthy_Status_YN[grep("Healthy Control|Healthy", meta$Healthy_Status, ignore.case = TRUE)] <- "YES"

    # 非健康样本（有明确疾病状态且不是NA/Uncategorized的）标记为 NO
    # 排除 NA, Uncategorized, Others 等不明确状态
    unclear_status <- c(NA, "", " ", "NA", "N/A", "n/a", "Uncategorized", "uncategorized",
                        "Others", "others")

    non_unclear <- !(meta$Healthy_Status %in% unclear_status) & !is.na(meta$Healthy_Status)
    non_healthy <- non_unclear & is.na(meta$Healthy_Status_YN)
    meta$Healthy_Status_YN[non_healthy] <- "NO"

    # 统计
    cat("  Healthy_Status processing:\n")
    cat("    YES (Healthy):", sum(meta$Healthy_Status_YN == "YES", na.rm = TRUE), "\n")
    cat("    NO (Disease):", sum(meta$Healthy_Status_YN == "NO", na.rm = TRUE), "\n")
    cat("    NA (Uncategorized/Unknown):", sum(is.na(meta$Healthy_Status_YN)), "\n")

    # 用新列替换旧列用于分析
    meta$Healthy_Status <- meta$Healthy_Status_YN
    meta$Healthy_Status_YN <- NULL
  }
  return(meta)
}

# --- 7. 循环处理每个采样部位 ---
for (file_name in sites_files) {
  site_name <- gsub("metadata_|_0127.csv", "", file_name)
  cat("\n========================================\n")
  cat("Processing:", site_name, "\n")
  cat("========================================\n")

  if (!file.exists(file_name)) {
    cat("File not found:", file_name, "\n")
    next
  }

  # 读取metadata
  meta <- read.csv(file_name, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(meta)[1] <- "SampleID"
  rownames(meta) <- meta$SampleID

  cat("Original samples:", nrow(meta), "\n")

  # 清洗数据
  meta <- clean_metadata(meta)

  # 处理Healthy_Status
  meta <- process_healthy_status(meta)

  # 检查目标列
  available_cols <- target_cols[target_cols %in% names(meta)]
  missing_cols <- target_cols[!target_cols %in% names(meta)]

  if (length(missing_cols) > 0) {
    cat("Warning: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  }

  cat("Using columns:", paste(available_cols, collapse = ", "), "\n")

  # 选择列并移除NA
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
    cat("Error: No samples left after cleaning!\n")
    next
  }

  # 显示各因子的水平数
  for (col in available_cols) {
    if (col %in% names(meta_clean)) {
      n_levels <- length(unique(meta_clean[[col]]))
      cat("  ", col, ":", n_levels, "levels")
      if (col == "Healthy_Status") {
        cat(" (", paste(unique(meta_clean[[col]]), collapse = ", "), ")")
      }
      cat("\n")
      if (n_levels < 2) {
        cat("    Warning:", col, "has only 1 level, may cause issues in adonis2\n")
      }
    }
  }

  # 样本对齐
  common_samples <- intersect(rownames(meta_clean), rownames(dist_big))
  cat("Matched with Distance Matrix:", length(common_samples), "samples\n")

  if (length(common_samples) < 10) {
    cat("Skipped: Too few samples (<10)\n")
    next
  }

  # 提取最终数据
  meta_final <- meta_clean[common_samples, , drop = FALSE]
  dist_final <- as.dist(dist_big[common_samples, common_samples])

  # 确保Bioproject作为分层变量可用
  if (!"Bioproject" %in% names(meta_final)) {
    if ("BioProject" %in% names(meta_final)) {
      meta_final$Bioproject <- meta_final$BioProject
    } else {
      cat("Error: Bioproject column not found for strata\n")
      next
    }
  }

  # --- 8. 运行 Adonis2 ---
  cat("\nRunning Adonis2 with formula:\n")
  cat("  dist ~ Sex + Age_Group + Seq_Region + Geography + Healthy_Status\n")
  cat("by = 'margin', strata = Bioproject\n\n")

  tryCatch({
    ado_res <- adonis2(dist_final ~ Sex + Age_Group + Seq_Region + Geography + Healthy_Status,
                       data = meta_final,
                       permutations = 999,
                       by = "margin",
                       strata = meta_final$Bioproject,
                       parallel = n_cores)

    print(ado_res)

    # 保存结果
    output_file <- paste0("adonis2_result_", site_name, "_0327.csv")
    write.csv(as.data.frame(ado_res), output_file)
    cat("\nResults saved to:", output_file, "\n")
    cat("Done!\n")

  }, error = function(e) {
    cat("Error in adonis2:", conditionMessage(e), "\n")
  })

  gc()
  cat("\n")
}

cat("\n========================================\n")
cat("All sites processed!\n")
cat("========================================\n")

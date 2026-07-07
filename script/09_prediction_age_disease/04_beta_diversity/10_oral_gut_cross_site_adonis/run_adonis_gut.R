library(vegan)
library(data.table)
library(tidyverse)

# --- 1. 安全并行设置 (8核) ---
# 降低核心数以保证内存安全，防止任务被 Kill
n_cores <- 8 
cat("Running with Safe Cores:", n_cores, "\n")
Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

work_dir <- "/home/hexincheng/task/huang_teacher/gut_adonis_analysis_0115/oral_gut_0127/gut"
setwd(work_dir)

# --- 2. 数据清洗 ---
cat("Reading Metadata...\n")
meta <- read.csv("metadata_gut_all_0127.csv", check.names = FALSE, row.names = 1)

# 【关键步骤】处理 Age 列
# 1. 强制转为数字 (非数字字符变 NA)
meta$Age <- suppressWarnings(as.numeric(as.character(meta$Age)))
# 2. 剔除负值或0 (逻辑上年龄应 > 0)
meta$Age[meta$Age <= 0] <- NA

target_cols <- c("BioProject", "Age", "Sex", "Age_Group", "Seq_Region", "Geography")

meta_clean <- meta %>% 
  select(all_of(target_cols)) %>% 
  na.omit() %>%       # 剔除任何含有 NA 的样本
  droplevels()        # 清除多余因子水平

cat("Raw Samples:", nrow(meta), "-> Cleaned Samples:", nrow(meta_clean), "\n")
cat("Dropped Samples:", nrow(meta) - nrow(meta_clean), "\n")

if(nrow(meta_clean) == 0) stop("Error: No samples left! Check data cleaning.")

# --- 3. 读取大矩阵 ---
cat("Reading Distance Matrix (fread)...\n")
dist_df <- fread("unifrac_weight.dm", header = TRUE, data.table = FALSE)
rownames(dist_df) <- dist_df[,1]
dist_df <- dist_df[,-1]

# --- 4. 样本对齐 ---
common_samples <- intersect(rownames(meta_clean), rownames(dist_df))
cat("Matched Samples:", length(common_samples), "\n")

if(length(common_samples) < 10) stop("Error: Too few matched samples!")

meta_final <- meta_clean[common_samples, , drop = FALSE]
dist_final <- as.dist(dist_df[common_samples, common_samples])

rm(dist_df, meta, meta_clean); gc()

# --- 5. 运行 Adonis ---
cat("Running Adonis2...\n")
tryCatch({
  ado_res <- adonis2(dist_final ~ BioProject + Age + Sex + Age_Group + Seq_Region + Geography, 
                     data = meta_final, 
                     permutations = 999, 
                     by = "margin",
                     parallel = n_cores)
  
  print(ado_res)
  write.csv(as.data.frame(ado_res), "adonis_gut_result_0127.csv")
  cat("Success!\n")
}, error = function(e) { print(e) })
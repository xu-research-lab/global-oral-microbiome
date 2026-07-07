library(vegan)
library(ggplot2)

# --- 1. 安全并行设置 (8核) ---
n_cores <- 8
cat("Running with Safe Cores:", n_cores, "\n")
Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

work_dir <- "/home/hexincheng/task/huang_teacher/gut_adonis_analysis_0115/oral_gut_0127/oral"
setwd(work_dir)

# --- 2. 加载 RData ---
cat("Loading unifrac.RData...\n")
loaded_objs <- load("unifrac.RData")
dist_big <- as.matrix(get(loaded_objs[1]))

sites_files <- c("metadata_Saliva_0127.csv", "metadata_Oral_Mucosa_0127.csv", 
                 "metadata_Tongue_0127.csv", "metadata_Tooth_0127.csv")
target_cols <- c("BioProject", "Age", "Sex", "Age_Group", "Seq_Region", "Geography")

for (file_name in sites_files) {
  site_name <- gsub("metadata_|_0127.csv", "", file_name)
  cat("\n=== Processing:", site_name, "===\n")
  
  if(!file.exists(file_name)) next
  
  meta <- read.csv(file_name, check.names = FALSE)
  colnames(meta)[1] <- "SampleID"
  rownames(meta) <- meta$SampleID
  
  # 【关键步骤】清洗 Age
  # 1. 强制转数字 (解决 "60+", "18-60" 等问题)
  meta$Age <- suppressWarnings(as.numeric(as.character(meta$Age)))
  # 2. 标记负值 (解决 -0.13 等问题)
  if(any(meta$Age < 0, na.rm = TRUE)) {
    cat("  Warning: Found negative Age values! Setting them to NA.\n")
    meta$Age[meta$Age < 0] <- NA
  }
  
  # 清洗 NA
  meta_clean <- meta[, target_cols, drop = FALSE]
  meta_clean <- na.omit(meta_clean)
  # 去除因子中的空水平
  for (col in names(meta_clean)) {
    if (is.factor(meta_clean[[col]])) {
      meta_clean[[col]] <- droplevels(meta_clean[[col]])
    }
  }
  
  cat("  Raw:", nrow(meta), "-> Cleaned:", nrow(meta_clean), "(Dropped:", nrow(meta)-nrow(meta_clean), ")\n")
  
  if(nrow(meta_clean) == 0) { cat("  Error: No samples left after cleaning!\n"); next }
  
  common_samples <- intersect(rownames(meta_clean), rownames(dist_big))
  cat("  Matched with Distance Matrix:", length(common_samples), "\n")
  
  if(length(common_samples) < 10) { cat("  Skipped: <10 samples\n"); next }
  
  meta_final <- meta_clean[common_samples, , drop = FALSE]
  dist_final <- as.dist(dist_big[common_samples, common_samples])
  
  cat("  Running Adonis2...\n")
  tryCatch({
    ado_res <- adonis2(dist_final ~ BioProject + Age + Sex + Age_Group + Seq_Region + Geography, 
                       data = meta_final, 
                       permutations = 999, 
                       by = "margin",
                       parallel = n_cores)
    
    write.csv(as.data.frame(ado_res), paste0("adonis_result_", site_name, "_0127.csv"))
    cat("  Done.\n")
  }, error = function(e) { cat("  Failed:", conditionMessage(e), "\n") })
  
  gc()
}
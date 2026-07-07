#!/usr/bin/env Rscript
# Oral Adonis2 Analysis Script (Bray-Curtis Distance)
# Formula: dist_matrix ~ Sex + Age_Group + Geography + Region_custom + Healthy_Status
# Healthy_Status: YES (健康) / NO (不健康)
# by = "margin", strata = metadata$BioProject

# --- 1. 加载库 ---
library(vegan)

# --- 2. 安全并行设置 ---
n_cores <- 8
cat("Running with Safe Cores:", n_cores, "\n")
Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

work_dir <- "/beegfs/huangxiaochang/global_oral/genus_correction"
setwd(work_dir)

# --- 3. 加载丰度数据 ---
cat("Loading genus_data_corrected_mmuphin.RData...\n")
load("genus_data_corrected_mmuphin.RData")
cat("Abundance matrix dimensions:", dim(genus_data_corrected), "\n")

# --- 4. 加载 metadata ---
cat("Loading qiime2_metadata.tsv...\n")
meta_all <- read.table("qiime2_metadata.tsv", header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "", fill = TRUE)
cat("Metadata dimensions:", dim(meta_all), "\n")

# --- 5. 定义要分析的采样部位 ---
sites <- c("Saliva", "Oral Mucosa", "Tongue", "Tooth")

# --- 6. 数据清洗函数 ---
clean_metadata <- function(meta) {
  for (col in names(meta)) {
    if (is.character(meta[[col]]) || is.factor(meta[[col]])) {
      meta[[col]] <- as.character(meta[[col]])
      meta[[col]][trimws(meta[[col]]) == ""] <- NA
      unclear_values <- c("Others_region", "Unknown", "unknown", "NA", "N/A", "n/a",
                          "uncategorized", "Uncategorized", "Others", "others",
                          "Others_body", "Other_body_site")
      meta[[col]][meta[[col]] %in% unclear_values] <- NA
    }
  }
  return(meta)
}

# --- 7. 处理 Healthy_Status 函数 ---
process_healthy_status <- function(meta) {
  if ("Healthy_Status" %in% names(meta)) {
    meta$Healthy_Status_YN <- NA
    meta$Healthy_Status_YN[grep("Healthy Control|Healthy", meta$Healthy_Status, ignore.case = TRUE)] <- "YES"
    unclear_status <- c(NA, "", " ", "NA", "N/A", "n/a", "Uncategorized", "uncategorized", "Others", "others")
    non_unclear <- !(meta$Healthy_Status %in% unclear_status) & !is.na(meta$Healthy_Status)
    non_healthy <- non_unclear & is.na(meta$Healthy_Status_YN)
    meta$Healthy_Status_YN[non_healthy] <- "NO"
    cat("  Healthy_Status processing:\n")
    cat("    YES (Healthy):", sum(meta$Healthy_Status_YN == "YES", na.rm = TRUE), "\n")
    cat("    NO (Disease):", sum(meta$Healthy_Status_YN == "NO", na.rm = TRUE), "\n")
    cat("    NA (Uncategorized/Unknown):", sum(is.na(meta$Healthy_Status_YN)), "\n")
    meta$Healthy_Status <- meta$Healthy_Status_YN
    meta$Healthy_Status_YN <- NULL
  }
  return(meta)
}

# --- 8. 创建 Geography 列（基于 Country_Cleaned 映射到洲） ---
create_geography <- function(country_vec) {
  geography_map <- list(
    Africa = c("Algeria", "Egypt", "Ethiopia", "Ghana", "Kenya", "Libya", "Morocco",
               "Nigeria", "South Africa", "Tanzania", "Tunisia", "Uganda", "Zimbabwe",
               "Burkina Faso", "Cameroon", "Cote d'Ivoire", "Democratic Republic of the Congo",
               "Mali", "Niger", "Senegal", "Somalia", "Sudan", "Rwanda", "Mauritius"),
    Asia = c("China", "Japan", "India", "Indonesia", "Pakistan", "Bangladesh", "Philippines",
             "Vietnam", "Turkey", "Iran", "Thailand", "Myanmar", "South Korea", "Iraq",
             "Afghanistan", "Saudi Arabia", "Uzbekistan", "Malaysia", "Nepal", "Yemen",
             "North Korea", "Sri Lanka", "Kazakhstan", "Syria", "Cambodia", "Israel",
             "Jordan", "Laos", "Lebanon", "Oman", "Kuwait", "Georgia", "Mongolia",
             "Armenia", "Qatar", "Bahrain", "Timor-Leste", "United Arab Emirates",
             "Kyrgyzstan", "Turkmenistan", "Singapore", "Tajikistan"),
    Europe = c("Germany", "United Kingdom", "France", "Italy", "Spain", "Poland", "Ukraine",
               "Netherlands", "Belgium", "Greece", "Portugal", "Czech Republic", "Romania",
               "Hungary", "Sweden", "Austria", "Serbia", "Switzerland", "Bulgaria", "Denmark",
               "Finland", "Slovakia", "Norway", "Ireland", "Croatia", "Lithuania", "Slovenia",
               "Latvia", "Estonia", "Moldova", "Albania", "North Macedonia", "Bosnia and Herzegovina",
               "Luxembourg", "Malta", "Iceland", "Montenegro", "Cyprus", "Belarus", "Russia"),
    North_America = c("United States", "Mexico", "Canada", "Guatemala", "Cuba", "Haiti",
                      "Dominican Republic", "Honduras", "Nicaragua", "El Salvador", "Costa Rica",
                      "Panama", "Jamaica", "Trinidad and Tobago", "Belize", "Barbados",
                      "Bahamas", "Saint Lucia", "Grenada", "Saint Vincent and the Grenadines"),
    South_America = c("Brazil", "Colombia", "Argentina", "Peru", "Venezuela", "Chile", "Ecuador",
                      "Bolivia", "Paraguay", "Uruguay", "Guyana", "Suriname", "French Guiana"),
    Oceania = c("Australia", "Papua New Guinea", "New Zealand", "Fiji", "Solomon Islands",
                "Vanuatu", "Samoa", "Kiribati", "Tonga", "Micronesia", "Palau", "Tuvalu",
                "Nauru", "Marshall Islands")
  )
  country_to_geo <- unlist(lapply(names(geography_map), function(geo) {
    setNames(rep(geo, length(geography_map[[geo]])), geography_map[[geo]])
  }))
  geography <- sapply(country_vec, function(country) {
    if (is.na(country) || country == "") return(NA)
    if (country %in% names(country_to_geo)) {
      return(country_to_geo[[country]])
    }
    return(NA)
  })
  return(geography)
}

# --- 9. 目标列 ---
target_cols <- c("BioProject", "Sex", "Age_Group", "Geography", "Region_custom", "Healthy_Status")

# --- 10. 循环处理每个采样部位 ---
results_list <- list()

for (site_name in sites) {
  cat("\n========================================\n")
  cat("Processing:", site_name, "\n")
  cat("========================================\n")

  # 过滤当前部位的样本
  meta <- meta_all[meta_all$Site_Cleaned == site_name, ]

  if (nrow(meta) == 0) {
    cat("No samples found for", site_name, "\n")
    next
  }

  cat("Original samples:", nrow(meta), "\n")

  # 创建 Geography 列
  cat("Creating Geography column from Country_Cleaned...\n")
  meta$Geography <- create_geography(meta$Country_Cleaned)

  # 清洗数据
  meta <- clean_metadata(meta)

  # 处理 Healthy_Status
  meta <- process_healthy_status(meta)

  # 添加 SampleName 用于匹配（在 clean 之后，因为 clean 可能修改列）
  meta$SampleName <- meta$X.SampleID

  # 检查目标列
  available_cols <- target_cols[target_cols %in% names(meta)]
  missing_cols <- target_cols[!target_cols %in% names(meta)]

  if (length(missing_cols) > 0) {
    cat("Warning: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  }

  cat("Using columns:", paste(available_cols, collapse = ", "), "\n")

  # 选择列并移除 NA（保留 SampleName 用于匹配）
  cols_to_keep <- c("SampleName", available_cols)
  meta_clean <- meta[, cols_to_keep, drop = FALSE]
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

  # 样本对齐 - 与丰度数据中的样本
  abundance_data <- genus_data_corrected
  common_samples <- intersect(meta_clean$SampleName, colnames(abundance_data))
  cat("Matched with Abundance Matrix:", length(common_samples), "samples\n")

  if (length(common_samples) < 10) {
    cat("Skipped: Too few samples (<10)\n")
    next
  }

  # 提取最终数据
  meta_matched <- meta_clean[meta_clean$SampleName %in% common_samples, , drop = FALSE]
  meta_matched <- meta_matched[match(common_samples, meta_matched$SampleName), , drop = FALSE]
  rownames(meta_matched) <- meta_matched$SampleName
  meta_final <- meta_matched[, available_cols, drop = FALSE]

  abundance_final <- abundance_data[, common_samples, drop = FALSE]

  # 转置使行为样本，列为物种
  abundance_final <- t(abundance_final)

  # --- 11. 计算 Bray-Curtis 距离 ---
  cat("Calculating Bray-Curtis distance...\n")
  dist_final <- vegdist(abundance_final, method = "bray")

  # 确保 BioProject 作为分层变量可用
  if (!"BioProject" %in% names(meta_final)) {
    cat("Error: BioProject column not found for strata\n")
    next
  }

  # --- 12. 运行 Adonis2 ---
  cat("\nRunning Adonis2 with formula:\n")
  cat("  dist ~ Sex + Age_Group + Geography + Region_custom + Healthy_Status\n")
  cat("by = 'margin', strata = BioProject\n\n")

  tryCatch({
    ado_res <- adonis2(dist_final ~ Sex + Age_Group + Geography + Region_custom + Healthy_Status,
                       data = meta_final,
                       permutations = 999,
                       by = "margin",
                       strata = meta_final$BioProject,
                       parallel = n_cores)

    print(ado_res)

    # 保存结果到列表
    results_list[[site_name]] <- ado_res

    # 保存结果（处理空格）
    output_file <- paste0("adonis2_result_", gsub(" ", "_", site_name), "_brays.csv")
    write.csv(as.data.frame(ado_res), output_file)
    cat("\nResults saved to:", output_file, "\n")
    cat("Done!\n")

  }, error = function(e) {
    cat("Error in adonis2:", conditionMessage(e), "\n")
  })

  gc()
  cat("\n")
}

# 保存所有结果
saveRDS(results_list, "adonis2_results_all_sites.rds")
cat("\nAll site results saved to: adonis2_results_all_sites.rds\n")

cat("\n========================================\n")
cat("All sites processed!\n")
cat("========================================\n")

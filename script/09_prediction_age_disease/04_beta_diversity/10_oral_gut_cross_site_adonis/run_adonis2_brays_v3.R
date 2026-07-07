#!/usr/bin/env Rscript
library(vegan)

n_cores <- 8
cat("Running with Safe Cores:", n_cores, "\n")
Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

setwd("/beegfs/huangxiaochang/global_oral/genus_correction")

cat("Loading genus_data_corrected_mmuphin.RData...\n")
load("genus_data_corrected_mmuphin.RData")
cat("Abundance matrix dimensions:", dim(genus_data_corrected), "\n")

cat("Loading qiime2_metadata.tsv...\n")
meta_all <- read.table("qiime2_metadata.tsv", header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "", fill = TRUE)
cat("Metadata dimensions:", dim(meta_all), "\n")

# 使用第一列（#SampleID）作为样本 ID（与丰度数据列名匹配）
colnames(meta_all)[1] <- "SampleID"
meta_all$SampleName <- meta_all$SampleID

sites <- c("Saliva", "Oral Mucosa", "Tongue", "Tooth")

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

target_cols <- c("BioProject", "Sex", "Age_Group", "Geography", "Sequencing_Region", "Healthy_Status")
results_list <- list()

for (site_name in sites) {
  cat("\n========================================\n")
  cat("Processing:", site_name, "\n")
  cat("========================================\n")

  # 过滤样本
  meta <- meta_all[meta_all$Site_Cleaned == site_name, ]

  if (nrow(meta) == 0) {
    cat("No samples found for", site_name, "\n")
    next
  }

  cat("Original samples:", nrow(meta), "\n")

  # 创建 Geography
  cat("Creating Geography column...\n")
  meta$Geography <- create_geography(meta$Country_Cleaned)

  # 按顺序处理：先 clean，然后 process_healthy_status
  # SampleName 已经在 meta_all 中添加了，过滤后会保留
  cat("SampleName in meta before clean:", "SampleName" %in% names(meta), "\n")

  cat("Running clean_metadata...\n")
  meta <- clean_metadata(meta)

  cat("Running process_healthy_status...\n")
  meta <- process_healthy_status(meta)

  # 检查 SampleName 是否存在
  cat("SampleName after processing:", "SampleName" %in% names(meta), "\n")

  # 如果 SampleName 不在，重新添加
  if (!"SampleName" %in% names(meta)) {
    cat("Re-adding SampleName...\n")
    meta$SampleName <- meta_all$SampleID[meta_all$Site_Cleaned == site_name]
    # 重新对齐行
    meta <- meta[order(rownames(meta)), ]
  }

  available_cols <- target_cols[target_cols %in% names(meta)]
  cat("Available cols:", paste(available_cols, collapse = ", "), "\n")

  cols_to_keep <- c("SampleName", available_cols)
  cat("Cols to keep:", paste(cols_to_keep, collapse = ", "), "\n")

  cat("Selecting columns...\n")
  meta_clean <- meta[, cols_to_keep, drop = FALSE]
  cat("Success! meta_clean dim:", dim(meta_clean), "\n")

  meta_clean <- na.omit(meta_clean)
  cat("After NA removal:", nrow(meta_clean), "samples\n")

  for (col in names(meta_clean)) {
    if (is.factor(meta_clean[[col]])) {
      meta_clean[[col]] <- droplevels(meta_clean[[col]])
    }
  }

  abundance_data <- genus_data_corrected
  common_samples <- intersect(meta_clean$SampleName, colnames(abundance_data))
  cat("Matched with Abundance Matrix:", length(common_samples), "samples\n")

  if (length(common_samples) < 10) {
    cat("Skipped: Too few samples (<10)\n")
    next
  }

  meta_matched <- meta_clean[meta_clean$SampleName %in% common_samples, , drop = FALSE]
  meta_matched <- meta_matched[match(common_samples, meta_matched$SampleName), , drop = FALSE]
  rownames(meta_matched) <- meta_matched$SampleName
  meta_final <- meta_matched[, available_cols, drop = FALSE]

  abundance_final <- abundance_data[, common_samples, drop = FALSE]
  abundance_final <- t(abundance_final)

  cat("Calculating Bray-Curtis distance...\n")
  dist_final <- vegdist(abundance_final, method = "bray")

  if (!"BioProject" %in% names(meta_final)) {
    cat("Error: BioProject column not found for strata\n")
    next
  }

  cat("Running Adonis2...\n")
  cat("  Formula: dist ~ Sex + Age_Group + Geography + Sequencing_Region + Healthy_Status\n")
  cat("  strata = BioProject\n")

  tryCatch({
    ado_res <- adonis2(dist_final ~ Sex + Age_Group + Geography + Sequencing_Region + Healthy_Status,
                       data = meta_final,
                       permutations = 999,
                       by = "margin",
                       strata = meta_final$BioProject,
                       parallel = n_cores)

    print(ado_res)
    results_list[[site_name]] <- ado_res

    output_file <- paste0("adonis2_result_", gsub(" ", "_", site_name), "_brays.csv")
    write.csv(as.data.frame(ado_res), output_file)
    cat("Results saved to:", output_file, "\n")
    cat("Done!\n")

  }, error = function(e) {
    cat("Error in adonis2:", conditionMessage(e), "\n")
  })

  gc()
}

saveRDS(results_list, "adonis2_results_all_sites.rds")
cat("\nAll site results saved!\n")

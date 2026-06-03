# ============================================================
# Tidyverse + Functionalized: Filter oral-gut pair metadata
# 说明：
# - metadata_oral_gut_pair_v2_merge.csv 和
#   metadata_pair_v2_status_unknown_merge.csv 中，
#   除了保留 sample_id 命中的 Oral 样本外，
#   还要把这些 Oral 样本对应的 pair_id 下的 Gut 样本一起保留。
# - metadata_pair_v2_additional_oral_site_merge.csv 保持原逻辑不变。
# - 新增：三个 _qc.csv 文件中的 Gut 样本还需通过 metadata_all_qc.tsv 的二次校验，
#   未通过的样本及其同 pair_id 的行将被删除。
# ============================================================

library(dplyr)
library(readr)

# ----------------------- 配置路径 -----------------------
out_dir <- "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/filter/filter-pair-0602/qc-pair-0603/results/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ref_file     <- "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/metadata_oral_qc.tsv"
ref_all_file <- "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/metadata_all_qc.tsv"

files <- list(
  pair = "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-gut-pair/260602-update/metadata_oral_gut_pair_v2_merge.csv",
  add  = "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-gut-pair/260602-update/metadata_pair_v2_additional_oral_site_merge.csv",
  unk  = "/home/yanglun/YSD/global_oral/metadata_v2_merge/oral-gut-pair/260602-update/metadata_pair_v2_status_unknown_merge.csv"
)

# ------------------ 读取参照文件 ------------------
valid_ids <- read_tsv(ref_file, show_col_types = FALSE) %>%
  pull(sample_id) %>%
  unique()

# 新增：metadata_all_qc.tsv 中 site_group_cleaned == "Gut" 的 sample_id
valid_gut_ids <- read_tsv(ref_all_file, show_col_types = FALSE) %>%
  filter(site_group_cleaned == "Gut") %>%
  pull(sample_id) %>%
  unique()

# ------------------ 定义处理函数 ------------------
# 函数1：打第一层标签
annotate_df <- function(filepath, valid_ids) {
  read_csv(filepath, show_col_types = FALSE) %>%
    mutate(filter_reason = ifelse(!(sample_id %in% valid_ids), "QC-filtered", NA_character_))
}

# 函数2：输出 annotation 文件
write_annotation <- function(df, base_name) {
  write_csv(df, file.path(out_dir, paste0(base_name, "_qc_annotation.csv")))
}

# 函数3：根据 valid_gut_ids 二次过滤
#   - 取出 site_group_cleaned == "Gut" 但 sample_id 不在 valid_gut_ids 中的行
#   - 将这些行对应的 pair_id 全部从该文件中移除
filter_by_gut <- function(df, valid_gut_ids) {
  invalid_pair_ids <- df %>%
    filter(site_group_cleaned == "Gut", !(sample_id %in% valid_gut_ids)) %>%
    pull(pair_id) %>%
    unique()

  df %>% filter(!(pair_id %in% invalid_pair_ids))
}

# 函数4：输出最终的清洁 _qc.csv 文件
write_qc <- function(df, base_name) {
  write_csv(df, file.path(out_dir, paste0(base_name, "_qc.csv")))
}

# ------------------ 执行主逻辑 ------------------
# 1. 标注所有文件
df_pair <- annotate_df(files$pair, valid_ids)
df_add  <- annotate_df(files$add, valid_ids)
df_unk  <- annotate_df(files$unk, valid_ids)

# 2. 识别"被保留的 Oral 样本"对应的 pair_id
rescued_pair_ids <- df_pair %>%
  filter(is.na(filter_reason), site_group_cleaned == "Oral") %>%
  pull(pair_id) %>%
  unique()

# 3. 对主文件和 unknown 文件执行 pair_id 级别救回
df_pair <- df_pair %>%
  mutate(filter_reason = ifelse(
    pair_id %in% rescued_pair_ids,
    NA_character_,
    filter_reason
  ))

df_unk <- df_unk %>%
  mutate(filter_reason = ifelse(
    is.na(status_cleaned),
    NA_character_,
    filter_reason
  ))

# 4. additional oral site 文件保持原逻辑

# 5. 输出 annotation 文件
write_annotation(df_pair, "metadata_oral_gut_pair_v2_merge")
write_annotation(df_add,  "metadata_pair_v2_additional_oral_site_merge")
write_annotation(df_unk,  "metadata_pair_v2_status_unknown_merge")

# 6. 生成清洁版（去掉 filter_reason 列）
df_pair_clean <- df_pair %>% filter(is.na(filter_reason)) %>% select(-filter_reason)
df_add_clean  <- df_add  %>% filter(is.na(filter_reason)) %>% select(-filter_reason)
df_unk_clean  <- df_unk  %>% filter(is.na(filter_reason)) %>% select(-filter_reason)

# 7. 新增：对清洁版进行 Gut 样本二次校验
df_pair_clean <- filter_by_gut(df_pair_clean, valid_gut_ids)
df_add_clean  <- filter_by_gut(df_add_clean,  valid_gut_ids)
df_unk_clean  <- filter_by_gut(df_unk_clean,  valid_gut_ids)

# 8. 输出最终的 _qc.csv 文件
write_qc(df_pair_clean, "metadata_oral_gut_pair_v2_merge")
write_qc(df_add_clean,  "metadata_pair_v2_additional_oral_site_merge")
write_qc(df_unk_clean,  "metadata_pair_v2_status_unknown_merge")

# ==========================================================
# 9. 生成 project_oral_gut_pair_list_qc.txt（基于二次过滤后的结果）
# ==========================================================
unique_bioprojects <- df_pair_clean %>%
  pull(pair_bioproject) %>%
  unique() %>%
  na.omit()

bioproject_list_path <- file.path(out_dir, "project_oral_gut_pair_list_qc.txt")
writeLines(as.character(unique_bioprojects), bioproject_list_path)
# ==========================================================

# ------------------ 生成结构化日志 ------------------
# 日志记录最终三个 _qc.csv 文件的数据统计情况
log_path <- file.path(out_dir, "qc_summary_log.txt")

summarize_df <- function(df, name) {
  n_total   <- nrow(df)
  n_oral    <- sum(df$site_group_cleaned == "Oral", na.rm = TRUE)
  n_gut     <- sum(df$site_group_cleaned == "Gut",  na.rm = TRUE)
  n_pair    <- length(unique(df$pair_id))
  n_proj    <- length(unique(na.omit(df$pair_bioproject)))
  paste0(
    "----- ", name, " -----\n",
    "  Total rows        : ", n_total, "\n",
    "  Oral samples      : ", n_oral,  "\n",
    "  Gut  samples      : ", n_gut,   "\n",
    "  Unique pair_id    : ", n_pair,  "\n",
    "  Unique bioproject : ", n_proj,  "\n"
  )
}

cat("==================================================\n", file = log_path)
cat(" Filter Log Generated on:", as.character(Sys.time()), "\n",
    file = log_path, append = TRUE)
cat(" Final QC files statistics (after Gut secondary validation)\n",
    file = log_path, append = TRUE)
cat("==================================================\n\n",
    file = log_path, append = TRUE)

cat(summarize_df(df_pair_clean, "metadata_oral_gut_pair_v2_merge_qc.csv"),
    "\n", file = log_path, append = TRUE)
cat(summarize_df(df_add_clean,  "metadata_pair_v2_additional_oral_site_merge_qc.csv"),
    "\n", file = log_path, append = TRUE)
cat(summarize_df(df_unk_clean,  "metadata_pair_v2_status_unknown_merge_qc.csv"),
    "\n", file = log_path, append = TRUE)

cat("==================================================\n",
    file = log_path, append = TRUE)
cat(" Bioproject list file: ", bioproject_list_path, "\n",
    file = log_path, append = TRUE)
cat(" Unique bioprojects   : ", length(unique_bioprojects), "\n",
    file = log_path, append = TRUE)
cat("==================================================\n",
    file = log_path, append = TRUE)

cat("处理完成：三个 _qc.csv 文件已通过 metadata_all_qc.tsv 的 Gut 样本二次校验。\n")

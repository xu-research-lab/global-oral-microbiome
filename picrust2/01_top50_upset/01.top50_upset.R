###############################################################################
# 01.top50_upset.R — 文档 §5（Top50 + UpSet），按指南模板用 ComplexUpset
# 用法: Rscript 01.top50_upset.R <pathway|EC|KO>
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(ComplexUpset); library(ggplot2)
})

type <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(type)) stop("Usage: Rscript 01.top50_upset.R <pathway|EC|KO>")

proj <- "/home/hexincheng/projects/hxc/oral_project/picrust2"
TOP_N <- 50

type_cfg <- list(
  pathway = list(abun = "pathways_out/path_abun_unstrat.tsv.gz",  short = "pathway"),
  EC      = list(abun = "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", short = "EC"),
  KO      = list(abun = "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", short = "KO")
)
cfg <- type_cfg[[type]]
if (is.null(cfg)) stop("type must be pathway, EC, or KO")

meta_file <- file.path(proj, "00.data/metadata_oral_qc_dedup.tsv")
abun_file <- file.path(proj, cfg$abun)
outdir    <- file.path(proj, paste0("downstream_", type, "_nsti2/01.top50_upset"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 辅助函数（与指南模板一致）----
read_function_table <- function(fp) {
  x <- fread(fp, data.table = FALSE)
  rownames(x) <- x[[1]]; x[[1]] <- NULL
  as.matrix(x)
}
tss <- function(mat) { cs <- colSums(mat); cs[cs == 0] <- 1; sweep(mat, 2, cs, "/") }

# ---- 1. 读取 ----
mat  <- read_function_table(abun_file)
meta <- fread(meta_file, data.table = FALSE, select = c("sample_id", "continent_group"))
rownames(meta) <- meta$sample_id

# ---- 2. 样本对齐 ----
common <- intersect(colnames(mat), rownames(meta))
mat  <- mat[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]

# ---- 3. 过滤 continent NA ----
keep <- !is.na(meta$continent_group) & meta$continent_group != "NA"
mat  <- mat[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
cat("Samples:", ncol(mat), " Features:", nrow(mat), "\n")

# ---- 4. TSS ----
mat_ra <- tss(mat)

# ---- 5. 每 continent mean RA ----
continents <- sort(unique(as.character(meta$continent_group)))
mean_list <- lapply(continents, function(cont) {
  s <- rownames(meta)[meta$continent_group == cont]
  data.frame(feature_id = rownames(mat_ra), continent_group = cont,
             mean_RA = rowMeans(mat_ra[, s, drop = FALSE]), stringsAsFactors = FALSE)
})
mean_df <- bind_rows(mean_list)

# ---- 6. 每 continent TopN ----
top_df <- mean_df %>%
  group_by(continent_group) %>%
  arrange(desc(mean_RA), .by_group = TRUE) %>%
  slice_head(n = TOP_N) %>% mutate(rank = row_number()) %>% ungroup()
fwrite(top_df, file.path(outdir, paste0("top", TOP_N, "_", cfg$short, "_by_continent.tsv")), sep = "\t")

# ---- 7. membership matrix ----
membership <- top_df %>%
  mutate(in_top = TRUE) %>%
  select(feature_id, continent_group, in_top) %>%
  pivot_wider(names_from = continent_group, values_from = in_top, values_fill = FALSE)
fwrite(membership, file.path(outdir, paste0("top", TOP_N, "_", cfg$short, "_membership_matrix.tsv")), sep = "\t")

# ---- 8. shared all / unique one ----
cont_cols <- setdiff(colnames(membership), "feature_id")
membership$n_continents <- rowSums(membership[, cont_cols])
shared_all <- membership %>% filter(n_continents == length(cont_cols)) %>% select(-n_continents)
unique_one <- membership %>% filter(n_continents == 1) %>% select(-n_continents)
fwrite(shared_all, file.path(outdir, paste0("top", TOP_N, "_", cfg$short, "_shared_all_continents.tsv")), sep = "\t")
fwrite(unique_one, file.path(outdir, paste0("top", TOP_N, "_", cfg$short, "_continent_unique.tsv")), sep = "\t")

# ---- 8b. shared 2-5 continents (intermediate overlap layers) ----
for (n_cont in 2:(length(cont_cols) - 1)) {
  shared_n <- membership %>% filter(n_continents == n_cont) %>% select(-n_continents)
  fwrite(shared_n, file.path(outdir, paste0("top", TOP_N, "_", cfg$short, "_shared_", n_cont, "_continents.tsv")), sep = "\t")
  cat("Shared", n_cont, "continents:", nrow(shared_n), "\n")
}

# ---- 9. UpSet plot (ComplexUpset, 与指南模板一致) ----
plot_df <- membership %>% select(feature_id, all_of(cont_cols))

p <- upset(
  plot_df,
  intersect = cont_cols,
  name = paste0("Top", TOP_N, " ", type),
  width_ratio = 0.2
) +
  ggtitle(paste0("Top", TOP_N, " ", type, " overlap by continent"))

ggsave(file.path(outdir, paste0("top", TOP_N, "_", cfg$short, "_UpSet.pdf")), p, width = 10, height = 6)
ggsave(file.path(outdir, paste0("top", TOP_N, "_", cfg$short, "_UpSet.png")), p, width = 10, height = 6, dpi = 300, type = "cairo", bg = "white")

cat("Shared all:", nrow(shared_all), " | Unique one:", nrow(unique_one), "\n")
cat("=== 01.top50_upset (", type, ") done ===\n", sep = "")

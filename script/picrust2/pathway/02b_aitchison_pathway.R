###############################################################################
# 02b_aitchison_pathway.R  —  §6.3 (Aitchison)
# CLR → Euclidean(Aitchison) → 距离矩阵(RDS) + PCoA(3轴)  [PERMANOVA已拆分]
# 内存优化同 02a
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(vegan); library(ggplot2); library(ggsci)
})

mem <- function(label) {
  g <- gc(reset = FALSE)
  cat(sprintf("[MEM] %-25s Vcells=%.0fMB Ncells=%.0fMB\n", label, g["Vcells","used"], g["Ncells","used"]))
}

pathway_file <- "/home/hexincheng/projects/hxc/oral_project/picrust2/pathways_out/path_abun_unstrat.tsv.gz"
meta_file    <- "/home/hexincheng/projects/hxc/oral_project/picrust2/00.data/metadata_oral_qc_dedup.tsv"
outdir       <- "/home/hexincheng/projects/hxc/oral_project/picrust2/downstream_pathway_nsti2/02.global_structure"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 读数据 ----
mat <- fread(pathway_file, data.table = FALSE)
rownames(mat) <- mat[[1]]; mat[[1]] <- NULL
mat <- as.matrix(mat)

meta <- fread(meta_file, data.table = FALSE,
              select = c("sample_id","continent_group","site_cleaned",
                         "status_cleaned","sequencing_region_cleaned","platform_cleaned"))
rownames(meta) <- meta$sample_id

# ---- status_broad ----
meta$status_broad <- ifelse(meta$status_cleaned == "HC", "HC",
                     ifelse(is.na(meta$status_cleaned) | meta$status_cleaned %in% c("NA",""),
                            NA, "Disease"))
meta$site_cleaned[meta$site_cleaned == "NA" | is.na(meta$site_cleaned)] <- "Other"
meta$sequencing_region_cleaned[meta$sequencing_region_cleaned == "NA" | is.na(meta$sequencing_region_cleaned)] <- "Unknown"
meta$platform_cleaned[meta$platform_cleaned == "NA" | is.na(meta$platform_cleaned)] <- "Unknown"

# ---- 样本对齐 + 过滤 ----
common <- intersect(colnames(mat), rownames(meta))
mat  <- mat[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
keep <- !is.na(meta$continent_group) & meta$continent_group != "NA"
mat  <- mat[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]

cov_cols <- c("continent_group","site_cleaned","status_broad",
               "sequencing_region_cleaned","platform_cleaned")
cc <- complete.cases(meta[, cov_cols])
mat  <- mat[, cc, drop = FALSE]
meta <- meta[cc, , drop = FALSE]
cat("Samples after complete-case:", ncol(mat), "| Pathways:", nrow(mat), "\n")
mem("After read & filter")

# 去零和样本
keep_zero <- colSums(mat) > 0
mat  <- mat[, keep_zero, drop = FALSE]
meta <- meta[keep_zero, , drop = FALSE]
cat("After removing zero-sum:", ncol(mat), "samples\n")
mem("After zero-sum filter")

# ---- CLR ----
clr <- function(m) {
  pc <- min(m[m > 0], na.rm = TRUE) / 2
  m[m == 0] <- pc
  lm <- log(m); sweep(lm, 2, colMeans(lm), "-")
}
mat_clr <- clr(mat)
rm(mat); gc()
mem("After CLR")

# ---- Aitchison distance (Euclidean on CLR) ----
dm_file <- file.path(outdir, "dist_aitchison.rds")
if (file.exists(dm_file)) {
  cat("--- Loading cached Aitchison distance ---\n")
  dm <- readRDS(dm_file); rm(mat_clr); gc()
} else {
  cat("--- Computing Aitchison distance ---\n")
  dm <- dist(t(mat_clr), method = "euclidean")
  rm(mat_clr); gc(); saveRDS(dm, dm_file, compress = FALSE)
  cat("--- Distance saved to RDS ---\n")
}
mem("After dist")

# ---- PCoA (skip if coords exist) ----
coords_file <- file.path(outdir, "pcoa_aitchison_coords.tsv")
if (file.exists(coords_file)) {
  cat("--- PCoA coords exist, skipping ---\n")
  pco_df <- fread(coords_file, data.table = FALSE)
  ve1 <- ve2 <- ve3 <- NA
} else {

cat("--- PCoA ---\n")
n <- attr(dm, "Size")
pco <- cmdscale(dm, k = 3, eig = FALSE)
eig1 <- sum(pco[,1,drop=FALSE]^2)
eig2 <- sum(pco[,2,drop=FALSE]^2)
eig3 <- sum(pco[,3,drop=FALSE]^2)
total_inertia <- sum(dm^2) / n
ve1 <- eig1 / total_inertia * 100
ve2 <- eig2 / total_inertia * 100
ve3 <- eig3 / total_inertia * 100
cat(sprintf("PCoA1: %.2f%%  PCoA2: %.2f%%  PCoA3: %.2f%%\n", ve1, ve2, ve3))

# 保存坐标
pco_df <- data.frame(sample_id = rownames(meta),
                     PC1 = pco[,1], PC2 = pco[,2], PC3 = pco[,3],
                     continent_group = meta$continent_group)
fwrite(pco_df, file.path(outdir, "pcoa_aitchison_coords.tsv"), sep = "\t")




}
rm(pco_df); gc()
mem("After PCoA plot")

# 距离矩阵已保存: dist_aitchison.rds
# PERMANOVA 请单独运行: 02_permanova_aitchison_pathway.R
cat("=== Distance + PCoA done (aitchison, pathway) ===\n")

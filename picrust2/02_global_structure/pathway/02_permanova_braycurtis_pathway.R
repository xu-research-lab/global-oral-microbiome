###############################################################################
# 02_permanova_braycurtis_pathway.R  —  PERMANOVA + betadisper only
# 加载已保存的距离矩阵 RDS，单独运行 PERMANOVA + betadisper
# 前提: dist_braycurtis.rds 已由距离计算脚本生成
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(vegan)
})

mem <- function(label) {
  g <- gc(reset = FALSE)
  cat(sprintf("[MEM] %-25s Vcells=%.0fMB Ncells=%.0fMB\n", label, g["Vcells","used"], g["Ncells","used"]))
}

meta_file <- "/home/hexincheng/projects/hxc/oral_project/picrust2/00.data/metadata_oral_qc_dedup.tsv"
outdir    <- "/home/hexincheng/projects/hxc/oral_project/picrust2/downstream_pathway_nsti2/02.global_structure"
dm_file   <- file.path(outdir, "dist_braycurtis.rds")

if (!file.exists(dm_file)) stop("Distance matrix not found: ", dm_file)

# ---- 读 metadata + 过滤（与距离脚本完全一致）----
meta <- fread(meta_file, data.table = FALSE,
              select = c("sample_id","continent_group","site_cleaned",
                         "status_cleaned","sequencing_region_cleaned","platform_cleaned"))
rownames(meta) <- meta$sample_id

meta$status_broad <- ifelse(meta$status_cleaned == "HC", "HC",
                     ifelse(is.na(meta$status_cleaned) | meta$status_cleaned %in% c("NA",""),
                            NA, "Disease"))
meta$site_cleaned[meta$site_cleaned == "NA" | is.na(meta$site_cleaned)] <- "Other"
meta$sequencing_region_cleaned[meta$sequencing_region_cleaned == "NA" | is.na(meta$sequencing_region_cleaned)] <- "Unknown"
meta$platform_cleaned[meta$platform_cleaned == "NA" | is.na(meta$platform_cleaned)] <- "Unknown"

cov_cols <- c("continent_group","site_cleaned","status_broad",
               "sequencing_region_cleaned","platform_cleaned")
cc <- complete.cases(meta[, cov_cols])
meta <- meta[cc, , drop = FALSE]
keep <- !is.na(meta$continent_group) & meta$continent_group != "NA"
meta <- meta[keep, , drop = FALSE]
cat("Metadata samples after filter:", nrow(meta), "\n")

# ---- 加载距离矩阵 ----
cat("Loading distance matrix:", dm_file, "\n")
dm <- readRDS(dm_file)
mem("After loading dist")

# ---- 对齐样本 ----
dm_labels <- labels(dm)
meta <- meta[dm_labels, , drop = FALSE]
cat("Samples for PERMANOVA:", nrow(meta), "\n")
mem("After align")

# ---- PERMANOVA ----
cat("--- PERMANOVA (Bray-Curtis) ---\n")
set.seed(42)
res <- adonis2(dm ~ continent_group + site_cleaned + status_broad +
                sequencing_region_cleaned + platform_cleaned,
               data = meta, permutations = 999, by = "margin",
               parallel = 14)
out <- data.frame(Term = rownames(res), Df = res$Df, SumOfSqs = res$SumOfSqs,
                  R2 = res$R2, F = res$F, `Pr(>F)` = res$`Pr(>F)`, check.names = FALSE)
fwrite(out, file.path(outdir, "permanova_braycurtis.tsv"), sep = "\t")
cat("=== PERMANOVA (Bray-Curtis) ===\n"); print(res); cat("\n")
rm(res, out); gc()
mem("After adonis2")

# ---- betadisper (custom: partial PCoA) ----
cat("--- betadisper (custom) ---\n")
betadisper_custom <- function(dm, group, k = 100, nperm = 999) {
  group <- as.factor(group)
  n <- attr(dm, "Size")
  groups <- levels(group)
  pco <- cmdscale(dm, k = min(k, n - 1), eig = FALSE)
  centroids <- t(sapply(groups, function(g) colMeans(pco[group == g, , drop = FALSE])))
  gidx <- as.integer(group)
  cent_exp <- centroids[gidx, , drop = FALSE]
  d2c <- sqrt(rowSums((pco - cent_exp)^2))
  obs_f <- anova(lm(d2c ~ group))$`F value`[1]
  perm_f <- numeric(nperm)
  for (i in seq_len(nperm)) {
    pg <- sample(group)
    pc <- t(sapply(groups, function(g) colMeans(pco[pg == g, , drop = FALSE])))
    pe <- pc[as.integer(pg), , drop = FALSE]
    pd <- sqrt(rowSums((pco - pe)^2))
    perm_f[i] <- anova(lm(pd ~ pg))$`F value`[1]
    rm(pe, pd, pc)
    if (i %% 200 == 0) gc()
  }
  pval <- (sum(perm_f >= obs_f, na.rm = TRUE) + 1) / (nperm + 1)
  grp_dist <- tapply(d2c, group, mean)
  list(F = obs_f, p = pval, group_dist = grp_dist,
       k = min(k, n - 1), nperm = nperm)
}

bd <- betadisper_custom(dm, meta$continent_group, k = 100, nperm = 999)
fwrite(data.frame(Groups = names(bd$group_dist), MeanDist = bd$group_dist),
       file.path(outdir, "betadisper_dist_braycurtis.tsv"), sep = "\t")
fwrite(data.frame(F = bd$F, p = bd$p, k = bd$k, nperm = bd$nperm),
       file.path(outdir, "betadisper_braycurtis.tsv"), sep = "\t")
cat(sprintf("betadisper F=%.4f, p=%.4f (k=%d axes)\n", bd$F, bd$p, bd$k))
rm(bd, dm); gc()
mem("After betadisper")

cat("=== PERMANOVA Bray-Curtis pathway done ===\n")

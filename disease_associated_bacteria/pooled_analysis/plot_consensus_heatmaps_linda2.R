#!/usr/bin/env Rscript
# Per-disease heatmaps for genera significant in all 3 methods with same Direction.
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(data.table)
})

WORK_DIR <- "/home/zhangcaiyi/project/global_oral_diseae_da/da_v2"
OUT_DIR <- file.path(WORK_DIR, "method_comparison_linda2")
HM_DIR <- file.path(OUT_DIR, "consensus_heatmaps")
dir.create(HM_DIR, showWarnings = FALSE, recursive = TRUE)

META_PATH <- paste0(
  "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/",
  "data/metadata-and-biom-v2-qc/metadata_oral_qc_dedup.tsv"
)
ABUND_PATH <- file.path(WORK_DIR, "genus_table_relabund.tsv")
OVERLAP_PATH <- file.path(OUT_DIR, "overlap_summary_3way.csv")

ov <- read.csv(OVERLAP_PATH, stringsAsFactors = FALSE)
parse_all3 <- function(s) {
  if (is.na(s) || !nzchar(s)) return(data.frame(Genus = character(), Direction = character()))
  parts <- strsplit(s, ";", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  spl <- strsplit(parts, "|", fixed = TRUE)
  data.frame(
    Genus = vapply(spl, `[`, character(1), 1),
    Direction = vapply(spl, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )
}

cons_rows <- list()
for (i in seq_len(nrow(ov))) {
  gdf <- parse_all3(ov$all3[[i]])
  if (nrow(gdf) == 0) next
  gdf$Disease <- ov$Disease[[i]]
  cons_rows[[length(cons_rows) + 1]] <- gdf[, c("Disease", "Genus", "Direction")]
}
cons_df <- if (length(cons_rows)) do.call(rbind, cons_rows) else {
  data.frame(Disease = character(), Genus = character(), Direction = character())
}
write.csv(cons_df, file.path(OUT_DIR, "consensus_all3_genera.csv"), row.names = FALSE)
cat(sprintf("Consensus genera (all3): %d rows across %d diseases\n",
            nrow(cons_df), length(unique(cons_df$Disease))))

meta <- read.delim(META_PATH, stringsAsFactors = FALSE, check.names = FALSE)
meta$sample_id <- trimws(as.character(meta$sample_id))
meta$bioproject_id <- trimws(as.character(meta$bioproject_id))
meta$status_cleaned <- trimws(as.character(meta$status_cleaned))
meta$status_cleaned[meta$status_cleaned %in% c("HC", "Healthy")] <- "Healthy Control"
keep_ids <- !is.na(meta$sample_id) & meta$sample_id != "" & !duplicated(meta$sample_id)
meta <- meta[keep_ids, , drop = FALSE]
valid <- !is.na(meta$bioproject_id) & meta$bioproject_id != "" &
  !is.na(meta$status_cleaned) & meta$status_cleaned != "" &
  meta$status_cleaned != "Other"
meta <- meta[valid, , drop = FALSE]
rownames(meta) <- meta$sample_id

clean_genus <- function(x) {
  x <- trimws(as.character(x))
  has_genus <- grepl("g__", x)
  extracted <- regmatches(x, regexpr("g__[^;]*", x))
  x[has_genus] <- sub("^g__", "", extracted)
  x <- sub("_.+$", "", x)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}

abund_raw <- as.data.frame(
  data.table::fread(ABUND_PATH, sep = "\t", header = TRUE, check.names = FALSE),
  stringsAsFactors = FALSE
)
raw_names <- abund_raw$Genus
abund_raw$Genus <- NULL
abund_mat <- as.matrix(abund_raw)
storage.mode(abund_mat) <- "numeric"
cleaned <- clean_genus(raw_names)
ok <- !is.na(cleaned) & cleaned != ""
abund_mat <- abund_mat[ok, , drop = FALSE]
cleaned <- cleaned[ok]
abund <- rowsum(abund_mat, group = cleaned, reorder = FALSE)

eligible_samples <- function(disease) {
  disease_bp <- unique(meta$bioproject_id[meta$status_cleaned == disease])
  keep <- meta$status_cleaned %in% c("Healthy Control", disease) &
    meta$bioproject_id %in% disease_bp
  sub <- meta[keep, , drop = FALSE]
  bp_tbl <- table(sub$bioproject_id, sub$status_cleaned)
  if (!("Healthy Control" %in% colnames(bp_tbl)) || !(disease %in% colnames(bp_tbl))) {
    return(NULL)
  }
  valid_bp <- rownames(bp_tbl)[bp_tbl[, "Healthy Control"] > 0 & bp_tbl[, disease] > 0]
  if (length(valid_bp) == 0) return(NULL)
  sub <- sub[sub$bioproject_id %in% valid_bp, , drop = FALSE]
  n_case <- sum(sub$status_cleaned == disease)
  n_ctrl <- sum(sub$status_cleaned == "Healthy Control")
  n_bp <- length(unique(sub$bioproject_id))
  if (n_case < 30 || n_ctrl < 30 || n_bp < 2) return(NULL)
  sub <- sub[sub$sample_id %in% colnames(abund), , drop = FALSE]
  if (nrow(sub) == 0) return(NULL)
  sub
}

diseases <- unique(cons_df$Disease)
hm_index <- list()

for (disease in diseases) {
  gdf <- cons_df[cons_df$Disease == disease, , drop = FALSE]
  genera <- unique(gdf$Genus)
  genera <- genera[genera %in% rownames(abund)]
  if (length(genera) == 0) {
    cat(sprintf("SKIP %s: no consensus genera in abundance table\n", disease))
    next
  }
  sub <- eligible_samples(disease)
  if (is.null(sub)) {
    cat(sprintf("SKIP %s: no eligible samples\n", disease))
    next
  }

  status_ord <- ifelse(sub$status_cleaned == "Healthy Control", 0L, 1L)
  ord <- order(sub$bioproject_id, status_ord, sub$sample_id)
  sub <- sub[ord, , drop = FALSE]
  sid <- sub$sample_id

  mat <- abund[genera, sid, drop = FALSE]
  mat_log <- log10(mat + 1e-6)
  is_healthy <- sub$status_cleaned == "Healthy Control"
  if (any(is_healthy) && any(!is_healthy)) {
    healthy_mean <- rowMeans(mat_log[, is_healthy, drop = FALSE], na.rm = TRUE)
  } else {
    healthy_mean <- rowMeans(mat_log, na.rm = TRUE)
  }
  mat_plot <- mat_log - healthy_mean

  dir_map <- setNames(gdf$Direction, gdf$Genus)
  row_dir <- dir_map[genera]
  row_ord <- order(row_dir, genera)
  mat_plot <- mat_plot[row_ord, , drop = FALSE]
  row_dir <- row_dir[row_ord]

  status_lab <- ifelse(sub$status_cleaned == "Healthy Control", "Healthy", "Disease")
  bp_levels <- sort(unique(sub$bioproject_id))
  set.seed(42)
  bp_cols <- setNames(
    circlize::rand_color(length(bp_levels), luminosity = "bright"),
    bp_levels
  )
  ha <- HeatmapAnnotation(
    Status = status_lab,
    BioProject = sub$bioproject_id,
    col = list(
      Status = c(Healthy = "#4DAF4A", Disease = "#E41A1C"),
      BioProject = bp_cols
    ),
    annotation_name_gp = gpar(fontsize = 9),
    simple_anno_size = unit(3, "mm"),
    show_legend = TRUE
  )
  ra <- rowAnnotation(
    Direction = row_dir,
    col = list(Direction = c("+" = "#377EB8", "-" = "#FF7F00")),
    annotation_name_gp = gpar(fontsize = 9),
    simple_anno_size = unit(3, "mm")
  )

  zmax <- max(abs(mat_plot), na.rm = TRUE)
  if (!is.finite(zmax) || zmax == 0) zmax <- 1
  col_fun <- colorRamp2(
    c(-zmax, 0, zmax),
    c("#2166AC", "#FFFFFF", "#B2182B")
  )

  disease_safe <- gsub(" ", "_", disease)
  out_png <- file.path(HM_DIR, paste0("consensus_heatmap_", disease_safe, ".png"))
  out_pdf <- file.path(HM_DIR, paste0("consensus_heatmap_", disease_safe, ".pdf"))

  n_row <- nrow(mat_plot)
  n_col <- ncol(mat_plot)
  fig_w <- min(18, max(8, 4 + n_col / 400))
  fig_h <- min(14, max(4, 2.5 + n_row * 0.28))

  ht <- Heatmap(
    mat_plot,
    name = "log10(rel)\nvs Healthy",
    col = col_fun,
    top_annotation = ha,
    right_annotation = ra,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_title = sprintf(
      "%s: all3 consensus genera (n=%d genera, %d samples; columns ordered by BioProject)",
      disease, n_row, n_col
    ),
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    use_raster = TRUE,
    raster_quality = 3,
    border = TRUE
  )

  png(out_png, width = fig_w, height = fig_h, units = "in", res = 200)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  pdf(out_pdf, width = fig_w, height = fig_h)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  hm_index[[length(hm_index) + 1]] <- data.frame(
    Disease = disease,
    n_genera = n_row,
    n_samples = n_col,
    png = file.path("consensus_heatmaps", basename(out_png)),
    stringsAsFactors = FALSE
  )
  cat(sprintf("Saved %s (%d genera x %d samples)\n", out_png, n_row, n_col))
}

hm_df <- if (length(hm_index)) do.call(rbind, hm_index) else {
  data.frame(Disease = character(), n_genera = integer(),
             n_samples = integer(), png = character())
}
write.csv(hm_df, file.path(OUT_DIR, "consensus_heatmap_index.csv"), row.names = FALSE)
cat("Done. Index: ", file.path(OUT_DIR, "consensus_heatmap_index.csv"), "\n", sep = "")

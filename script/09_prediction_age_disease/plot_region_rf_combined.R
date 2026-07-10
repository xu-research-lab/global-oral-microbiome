#!/usr/bin/env Rscript
# ============================================================================
# RF v2 — Region prediction heatmaps (Oral + Gut combined)
# Style reference: randomforest/rf_heatmap_region.R
# Uses data.table for fast processing of large gut prediction files.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/randomforest_v2")

OUT_DIR <- "."
ORAL_OVA_PRED  <- "region_ova_predictions.csv"
ORAL_OVO_PRED  <- "region_ovo_predictions.csv"
GUT_OVA_FILE   <- "/beegfs/huangxiaochang/global_oral/randomforest/gut_analysis/gut_region_RF_one_vs_all.csv"
GUT_OVO_FILE   <- "/beegfs/huangxiaochang/global_oral/randomforest/gut_analysis/gut_region_RF_one_vs_one.csv"

# ============================================================================
# Helper: fast AUC (Wilcoxon-Mann-Whitney)
# ============================================================================
fast_auc <- function(labels, scores) {
  labels <- as.integer(labels)
  if (length(unique(labels)) < 2L) return(NA_real_)
  n_pos <- sum(labels == 1L)
  n_neg <- sum(labels == 0L)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  r <- rank(scores, ties.method = "average")
  sum(r[labels == 1L]) / (n_pos * n_neg) - (n_pos + 1) / (2 * n_neg)
}

region_order <- c("North America", "Europe", "Asia", "Oceania", "South America", "Africa")

# ============================================================================
# 1. Pre-aggregate: compute fold-level AUC (cache to avoid re-reading)
# ============================================================================
ova_cache <- "region_ova_fold_auc.rds"
ovo_cache <- "region_ovo_fold_auc.rds"

if (!file.exists(ova_cache) || !file.exists(ovo_cache)) {
  cat("=== Computing fold-level AUC (this may take a moment) ===\n")

  # --- Oral OVA ---
  cat("  Oral OVA...\n")
  oral_ova <- fread(ORAL_OVA_PRED, showProgress = FALSE)
  oral_ova_fold <- oral_ova[!is.na(true_label) & !is.na(pred_prob),
    .(fold_auc = fast_auc(true_label, pred_prob)),
    by = .(site, target_region, fold)]
  oral_ova_fold <- oral_ova_fold[!is.na(fold_auc)]

  # --- Gut OVA (data.table for speed) ---
  cat("  Gut OVA...\n")
  gut_ova <- fread(GUT_OVA_FILE, showProgress = FALSE,
                   select = c("fold", "true_label", "pred_prob", "region"))
  setnames(gut_ova, old = "region", new = "target_region")
  gut_ova[, site := "Gut"]
  gut_ova_fold <- gut_ova[!is.na(true_label) & !is.na(pred_prob),
    .(fold_auc = fast_auc(true_label, pred_prob)),
    by = .(site, target_region, fold)]
  gut_ova_fold <- gut_ova_fold[!is.na(fold_auc)]
  rm(gut_ova); invisible(gc())

  # --- Combine OVA ---
  ova_fold_all <- rbindlist(list(oral_ova_fold, gut_ova_fold), use.names = TRUE)
  saveRDS(ova_fold_all, ova_cache)
  cat(sprintf("  OVA: %d fold-level AUCs saved\n", nrow(ova_fold_all)))
  rm(oral_ova, oral_ova_fold, gut_ova_fold); invisible(gc())

  # --- Oral OVO ---
  cat("  Oral OVO...\n")
  oral_ovo <- fread(ORAL_OVO_PRED, showProgress = FALSE)
  oral_ovo_fold <- oral_ovo[!is.na(true_label) & !is.na(pred_prob),
    .(fold_auc = fast_auc(true_label, pred_prob)),
    by = .(site, region_1, region_2, fold)]
  oral_ovo_fold <- oral_ovo_fold[!is.na(fold_auc)]

  # --- Gut OVO ---
  cat("  Gut OVO...\n")
  gut_ovo <- fread(GUT_OVO_FILE, showProgress = FALSE,
                   select = c("fold", "true_label", "pred_prob", "region_1", "region_2"))
  gut_ovo[, site := "Gut"]
  gut_ovo_fold <- gut_ovo[!is.na(true_label) & !is.na(pred_prob),
    .(fold_auc = fast_auc(true_label, pred_prob)),
    by = .(site, region_1, region_2, fold)]
  gut_ovo_fold <- gut_ovo_fold[!is.na(fold_auc)]
  rm(gut_ovo); invisible(gc())

  # --- Combine OVO ---
  ovo_fold_all <- rbindlist(list(oral_ovo_fold, gut_ovo_fold), use.names = TRUE)
  saveRDS(ovo_fold_all, ovo_cache)
  cat(sprintf("  OVO: %d fold-level AUCs saved\n", nrow(ovo_fold_all)))
  rm(oral_ovo, oral_ovo_fold, gut_ovo_fold); invisible(gc())
} else {
  cat("=== Loading cached fold-level AUC ===\n")
  ova_fold_all <- readRDS(ova_cache)
  ovo_fold_all <- readRDS(ovo_cache)
}

# ============================================================================
# 2. Aggregate to mean AUC
# ============================================================================
cat("\n=== Aggregating to mean AUC ===\n")

ova_combined <- ova_fold_all[, .(
  mean_auc = mean(fold_auc, na.rm = TRUE),
  n_folds   = .N
), by = .(site, region = target_region)] %>%
  mutate(
    site   = factor(site, levels = c("Saliva", "Tongue", "Tooth", "Gut")),
    region = factor(region, levels = region_order),
    label_text = sprintf("%.2f", mean_auc)
  ) %>%
  filter(!is.na(region))
cat(sprintf("  OVA: %d site×region combos\n", nrow(ova_combined)))

ovo_agg <- ovo_fold_all[, .(
  mean_auc = mean(fold_auc, na.rm = TRUE),
  n_folds   = .N
), by = .(site, region_1, region_2)]

# Make symmetric
ovo_symmetric <- rbindlist(list(
  ovo_agg,
  ovo_agg[, .(site, region_1 = region_2, region_2 = region_1, mean_auc, n_folds)]
)) %>%
  distinct() %>%
  mutate(
    site     = factor(site, levels = c("Saliva", "Tongue", "Tooth", "Gut")),
    region_1 = factor(region_1, levels = region_order),
    region_2 = factor(region_2, levels = region_order)
  ) %>%
  filter(!is.na(region_1), !is.na(region_2)) %>%
  # Remove unreliable Saliva South America–Africa pair (AUC=0.24)
  filter(!(site == "Saliva" & (
    (region_1 == "South America" & region_2 == "Africa") |
    (region_1 == "Africa" & region_2 == "South America")
  )))
cat(sprintf("  OVO: %d site×pair combos\n", nrow(ovo_agg)))

# ============================================================================
# 3. OVA heatmap (reference style: geom_tile + text)
# ============================================================================
cat("\n=== Plot OVA heatmap ===\n")

p_ova <- ggplot(ova_combined, aes(x = site, y = region, fill = mean_auc)) +
  geom_tile(color = "white", size = 0.8) +
  geom_text(aes(label = label_text), color = "black", size = 3.8, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#F7FBFF", "#D2E3F3", "#6BAED6", "#3182BD", "#08519C"),
    limits = c(0, 1), name = "Mean AUC"
  ) +
  scale_y_discrete(limits = rev(region_order)) +
  labs(title = "Geography RF — One-vs-All AUC (Oral + Gut)",
       subtitle = "5-fold cross-validation, stratified by BioProject",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"),
    axis.text.y = element_text(size = 11),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  )

ggsave(file.path(OUT_DIR, "region_ova_heatmap_combined.pdf"), p_ova,
       width = 6.5, height = 4.5)
ggsave(file.path(OUT_DIR, "region_ova_heatmap_combined.png"), p_ova,
       width = 6.5, height = 4.5, dpi = 300)
cat("  Saved region_ova_heatmap_combined.{pdf,png}\n")

# ============================================================================
# 4. OVO heatmaps: per-site + combined faceted
# ============================================================================
cat("\n=== Plot OVO heatmaps ===\n")

all_sites <- c("Saliva", "Gut")
sites_with_data <- intersect(all_sites, unique(as.character(ovo_symmetric$site)))

all_regions_ovo <- region_order[region_order %in% unique(c(
  as.character(ovo_symmetric$region_1),
  as.character(ovo_symmetric$region_2)
))]

# --- Per-site heatmaps ---
for (st in all_sites) {
  if (!st %in% sites_with_data) next

  site_df <- ovo_symmetric[ovo_symmetric$site == st]

  # Only keep regions present
  pr <- all_regions_ovo[all_regions_ovo %in% unique(c(
    as.character(site_df$region_1), as.character(site_df$region_2)))]

  grid_df <- expand.grid(
    region_1 = factor(pr, levels = pr),
    region_2 = factor(pr, levels = pr),
    stringsAsFactors = FALSE
  ) %>%
    left_join(site_df %>% as.data.frame(), by = c("region_1", "region_2")) %>%
    mutate(
      mean_auc  = ifelse(region_1 == region_2, NA_real_, mean_auc),
      label_txt = ifelse(region_1 == region_2, "",
                         ifelse(is.na(mean_auc), "",
                                sprintf("%.2f\n(%d)", mean_auc,
                                        ifelse(is.na(n_folds), 0, n_folds))))
    )

  p <- ggplot(grid_df, aes(x = region_2, y = region_1, fill = mean_auc)) +
    geom_tile(color = "white", size = 0.8) +
    geom_text(aes(label = label_txt), color = "black", size = 2.5,
              fontface = "bold", lineheight = 0.85) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0.5, limits = c(0, 1), na.value = "grey90", name = "Mean AUC"
    ) +
    labs(
      title = paste("One-vs-One Geography RF —", st),
      subtitle = "Train (row) → Predict (col); 5-fold CV stratified by BioProject",
      x = "Predicted Region", y = "Training Region",
      caption = "Mean AUC (n folds)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed()

  safe_name <- tolower(gsub(" ", "_", st))
  w <- 5 + length(pr) * 0.8
  ggsave(file.path(OUT_DIR, paste0("region_ovo_heatmap_", safe_name, ".pdf")), p,
         width = w, height = w - 0.5)
  ggsave(file.path(OUT_DIR, paste0("region_ovo_heatmap_", safe_name, ".png")), p,
         width = w, height = w - 0.5, dpi = 300)
  cat(sprintf("  Saved region_ovo_heatmap_%s.{pdf,png}\n", safe_name))
}

# --- Saliva + Gut lower-triangle, side by side ---
cat("\n=== Plot OVO: Saliva + Gut lower-triangle ===\n")

ovo_lt <- ovo_symmetric %>%
  filter(site %in% c("Saliva", "Gut")) %>%
  mutate(
    region_1 = factor(region_1, levels = all_regions_ovo),
    region_2 = factor(region_2, levels = all_regions_ovo)
  ) %>%
  filter(as.numeric(region_1) > as.numeric(region_2)) %>%
  mutate(
    label_txt = ifelse(is.na(mean_auc), "",
                       sprintf("%.2f", mean_auc))
  )

p_ovo_lt <- ggplot(ovo_lt, aes(x = region_1, y = region_2, fill = mean_auc)) +
  geom_tile(color = "white", size = 0.8) +
  geom_text(aes(label = label_txt), color = "black", size = 3.2, fontface = "bold") +
  facet_wrap(~ site, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1), na.value = "grey90", name = "Mean AUC"
  ) +
  scale_x_discrete(limits = all_regions_ovo, drop = FALSE) +
  scale_y_discrete(limits = rev(all_regions_ovo), drop = FALSE) +
  labs(
    title = "One-vs-One Geography RF — Saliva & Gut",
    subtitle = "Train (row) → Predict (col); lower triangle; 5-fold CV",
    x = "Predicted Region", y = "Training Region"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed()

n_pairs <- length(unique(ovo_lt$region_1))
ggsave(file.path(OUT_DIR, "region_ovo_heatmap_combined.pdf"), p_ovo_lt,
       width = 4.2 + n_pairs * 1.2, height = 3.0 + n_pairs * 0.55)
ggsave(file.path(OUT_DIR, "region_ovo_heatmap_combined.png"), p_ovo_lt,
       width = 4.2 + n_pairs * 1.2, height = 3.0 + n_pairs * 0.55, dpi = 300)
cat("  Saved region_ovo_heatmap_combined.{pdf,png} (Saliva + Gut, lower-triangle)\n")

# ============================================================================
# 5. Export summary tables
# ============================================================================
fwrite(ova_combined %>% as.data.table(), file.path(OUT_DIR, "region_ova_combined_summary.csv"))
fwrite(ovo_agg %>% as.data.table(), file.path(OUT_DIR, "region_ovo_combined_summary.csv"))

cat("\n=== Region RF combined plots done ===\n")

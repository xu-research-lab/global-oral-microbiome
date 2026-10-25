# ==============================================================================
# 03-plot-site-ova-bubble.R  (site OVA da_v2)
# Consensus DA figures for the site OVA pipeline. Reads the consensus tables
# produced by stage 02 (no re-run of the DA pipeline). Comparisons are OVA:
# each valid site (target) vs the pooled "Other" of the other valid sites,
# name "<site>_enrich".
#
#   Figure A: horizontal counts of >=2-method consensus genera
#   Figure B: DA bubble (rows = genera, cols = site targets;
#             size = -log10(Maaslin2 q), color = Maaslin2 coefficient),
#             at most 10 uniquely peak-abundant genera per target, selected by
#             significant positive coefficient with other-method support.
#             Group blocks have white gaps; all three panels have borders.
#             Mean relative abundance (z-score/raw) + prevalence accompany DA.
#   Figure C: top-5 consensus genera, log10(relabund) ridgelines per site.
#
# Input:  <out>/consensus_all.tsv / consensus_all_enhanced.tsv
#         data/metadata-and-biom-v2-qc/genus_table_relabund.tsv
# Default directory: <out>/figures_target_top10_globalpeak_maaslin2coef/
# Legacy --select=block55 output:
#         <out>/figures/figureA_consensus_counts_by_target.{pdf,png}
#         <out>/figures/figureB_da_bubble_{2method,3method}.{pdf,png}
#         <out>/figures/figureB_abundance_heatmap_{2method,3method}_{zscore,raw}.{pdf,png}
#         <out>/figures/figureB_prevalence_{2method,3method}.{pdf,png}
#         <out>/figures/figureB_abundance_{2method,3method}_combined{,_zscore,_raw}.{pdf,png}
#         <out>/figures/figureC_top5_abundance_{2method,3method}.{pdf,png}
#
# Usage (from project root):
#   Rscript differential_abundance/site_differential_ova/03-plot-site-ova-bubble.R
#   Rscript differential_abundance/site_differential_ova/03-plot-site-ova-bubble.R --scope=hc
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(data.table)
  library(patchwork); library(viridis); library(ggridges); library(grid)
})

source("differential_abundance/site_differential_ova/da_helpers.R")   # rd(), SITE_ORDER, parse_scope_arg, parse_ova_comparison, valid_groups
source("differential_abundance/real_lfc_helper.R")                   # real_lfc_matrix(), REAL_LFC_PS

source("differential_abundance/ova_plot_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
scope <- parse_scope_arg(args)

arg_val <- function(args, prefix, default) {
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) default else sub(prefix, "", hit[length(hit)])
}
# --select=target_topn (default): unique global mean-abundance peak,
# significant positive MaAsLin2 + other-method support; top 10 by coefficient.
# --select=block55 keeps the legacy per-target 5+5 selection.
select_mode <- arg_val(args, "--select=", "target_topn")
min_targets <- as.integer(arg_val(args, "--min-targets=", "3"))
top_n       <- as.integer(arg_val(args, "--top-n=", "-1"))
each_n      <- as.integer(arg_val(args, "--each-n=", "10"))
heatmap_mode <- arg_val(args, "--heatmap=", "zscore")

root <- Sys.getenv("DA_RESULTS_ROOT"); if (!nzchar(root)) root <- "differential_abundance/result"
out_dir <- file.path(site_root(scope), "site_da_consensus")
scope_tag <- if (scope == "hc") "site OVA (HC only)" else "site OVA (all)"

group_order <- SITE_ORDER   # canonical site order
mode_tag <- switch(select_mode,
  recurrence = paste0("_recurrence_min", min_targets,
                      if (top_n > 0) paste0("_top", top_n) else ""),
  target_topn = paste0("_target_top", each_n, "_globalpeak_maaslin2coef"),
  "")
figures_dir <- file.path(out_dir, paste0("figures", mode_tag))
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# OVA comparison "<site>_enrich" -> column label "<site> vs rest"
short_target <- function(g) paste0(g, " vs rest")

save_both <- function(path_pdf, plot, width, height, dpi = 300) {
  ggsave(path_pdf, plot, device = grDevices::cairo_pdf,
         width = width, height = height, limitsize = FALSE)
  png_path <- sub("\\.pdf$", ".png", path_pdf)
  ggsave(png_path, plot, device = ragg::agg_png,
         width = width, height = height, dpi = dpi, bg = "white",
         limitsize = FALSE)
}

cat("==========================================\n")
cat("Consensus DA plots — site OVA (each site vs Other)\n")
cat("Scope: ", scope_tag, "\n", sep="")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ==========================================================================
# 1. Load consensus tables
# ==========================================================================
cons_file <- file.path(out_dir, "consensus_all.tsv")
if (!file.exists(cons_file)) stop("Missing: ", cons_file)
cons <- fread(cons_file, data.table = FALSE)
cat("    Consensus rows:", nrow(cons), "\n")

enh_file <- file.path(out_dir, "consensus_all_enhanced.tsv")
if (!file.exists(enh_file)) stop("Missing: ", enh_file)
enh <- fread(enh_file, data.table = FALSE)
cat("    Enhanced rows:", nrow(enh), "\n")

# OVA comparison -> target (no nested stratum for site OVA; site=NA)
parsed <- lapply(enh$comparison, parse_ova_comparison)
enh$site   <- vapply(parsed, `[[`, character(1), "target")  # target = the site
enh$target <- enh$site
enh$pair_lab <- ifelse(!is.na(enh$target), short_target(enh$target), enh$comparison)

parsed_c <- lapply(cons$comparison, parse_ova_comparison)
cons$target <- vapply(parsed_c, `[[`, character(1), "target")
cons$pair_lab <- ifelse(!is.na(cons$target), short_target(cons$target), cons$comparison)

if (nrow(cons) == 0) stop("No consensus rows; check --scope.")

# ==========================================================================
# 2. Figure A: horizontal counts for >=2-method consensus
# ==========================================================================
# Derive tested targets from the full table, retaining genuine zero counts
# without treating an untested target as zero.
full_counts <- fread(file.path(out_dir, "full_table.tsv"), data.table = FALSE)
tested_targets <- unique(vapply(lapply(unique(full_counts$comparison),
                                      parse_ova_comparison), `[[`, character(1), "target"))
tested_targets <- group_order[group_order %in% tested_targets]
if (!length(tested_targets)) stop("No tested OVA targets found in full_table.tsv")
# Counts are computed below after loading the same group means as Figure B.

# ==========================================================================
# 3. Figure B: DA bubble + abundance heatmap + prevalence
# ==========================================================================
genus_label <- function(feature) {
  g <- make_genus_short(feature)
  if (nzchar(g)) return(g)
  parts <- strsplit(feature, ";", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  lab <- ""
  for (p in rev(parts)) {
    lab <- sub("^[a-z]__", "", p)
    if (nzchar(lab)) break
  }
  if (!nzchar(lab)) "unclassified" else lab
}
enh$genus_lab <- make.unique(vapply(enh$feature, genus_label, character(1)))

# per-target 5+5: target-enriched (real LFC>0) top 5 + Other-enriched (<0) top 5
# by |real LFC|; union across targets, first-come dedupe; block order = site
# order; target-enriched first within each block.
select_heatmap_features <- function(enh, min_methods, real_lfc_df,
                                    n_each = 5L) {
  keep <- enh[enh$n_methods >= min_methods & !is.na(enh$maaslin2_padj) &
              !is.na(enh$maaslin2_lfc), , drop = FALSE]
  keep$real_lfc <- real_lfc_df$real_lfc[match(
    paste(keep$feature, keep$target),
    paste(real_lfc_df$feature, real_lfc_df$target))]
  keep <- keep[!is.na(keep$real_lfc), , drop = FALSE]

  comps <- data.frame(target = unique(keep$target), stringsAsFactors = FALSE)
  comps$tgt_f <- factor(comps$target, levels = group_order)
  comps <- comps[order(comps$tgt_f), , drop = FALSE]

  ordered_rows <- character(0); block_of <- character(0)
  for (i in seq_len(nrow(comps))) {
    cc <- keep[keep$target == comps$target[i], , drop = FALSE]
    up <- cc[which(cc$real_lfc > 0), , drop = FALSE]
    dn <- cc[which(cc$real_lfc < 0), , drop = FALSE]
    up <- up[order(-up$real_lfc), , drop = FALSE]
    dn <- dn[order(dn$real_lfc), , drop = FALSE]
    take <- c(head(up$feature, n_each), head(dn$feature, n_each))
    take <- take[!take %in% ordered_rows]
    if (length(take) > 0) {
      ordered_rows <- c(ordered_rows, take)
      block_of <- c(block_of, rep(comps$target[i], length(take)))
    }
  }
  attr(ordered_rows, "block") <- block_of
  ordered_rows
}

select_recurrence_features <- function(enh, min_methods, min_targets, top_n = -1L) {
  keep <- enh[enh$n_methods >= min_methods & !is.na(enh$maaslin2_padj) &
              !is.na(enh$maaslin2_lfc), , drop = FALSE]
  if (nrow(keep) == 0) return(character(0))
  agg <- keep %>%
    group_by(feature) %>%
    summarise(n_target = n_distinct(target),
              qmin = min(maaslin2_padj, na.rm = TRUE), .groups = "drop") %>%
    filter(n_target >= min_targets) %>%
    arrange(desc(n_target), qmin)
  if (top_n > 0) agg <- head(agg, top_n)
  agg$feature
}

# target_topn selection is shared in differential_abundance/ova_plot_helpers.R.

# ---- metadata + relabund ----
meta <- load_age_metadata(hc_only = (scope == "hc"))
meta <- meta[meta[[COL_SITE]] %in% group_order, , drop = FALSE]
meta$group_lab <- meta[[COL_SITE]]
group_lab_order <- group_order[group_order %in% meta$group_lab]

read_relabund_rows <- function(feats) {
  pats <- paste0(feats, "\t")
  cmd <- paste0("head -1 ", RELABUND_TSV, "; grep -F -e '",
                paste(pats, collapse = "' -e '"), "' ", RELABUND_TSV)
  fread(cmd = cmd, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
}

site_colors <- setNames(viridis(length(group_order)), group_order)

# ---- real relabund LFC for every robust feature (top-level, shared) ----
# "Other" pool rebuilt with the SAME validity rule as the DA driver:
# valid_groups = sites with >= MIN_SAMPLES_PER_GROUP samples.
valid_sites <- valid_groups(meta, MIN_SAMPLES_PER_GROUP)
robust_feats <- unique(c(enh$feature[enh$n_methods >= 2],
                         cons$feature[cons$n_methods >= 2]))
rel_all <- read_relabund_rows(robust_feats)
samples_all <- intersect(colnames(rel_all), rownames(meta))
sdf_all <- as.data.frame(t(rel_all[, samples_all]))
colnames(sdf_all) <- rel_all[[1]]
site_vec_all <- meta[rownames(sdf_all), COL_SITE]
comps_all <- data.frame(target = unique(enh$target), stringsAsFactors = FALSE)
comps_all <- comps_all[!is.na(comps_all$target), , drop = FALSE]
comps_all$tgt_f <- factor(comps_all$target, levels = group_order)
comps_all <- comps_all[order(comps_all$tgt_f), , drop = FALSE]
# Include all plotted groups, even groups with no significant DA comparison.
mean_abundance_df <- do.call(rbind, lapply(group_lab_order, function(tgt) {
  idx <- which(meta[rownames(sdf_all), "group_lab"] == tgt)
  if (length(idx) == 0L) stop("No abundance samples for group: ", tgt)
  means <- colMeans(sdf_all[idx, , drop = FALSE])
  data.frame(feature = names(means), target = tgt,
             target_mean_abundance = as.numeric(means))
}))
# Raw mean extrema and significant model directions must agree.
for (min_count_methods in c(2L, 3L)) {
  ver_count <- paste0(min_count_methods, "method")
  count_result <- ova_extreme_counts(enh, mean_abundance_df, group_lab_order,
                                     tested_targets, min_count_methods, alpha = DA_ALPHA)
  fwrite(count_result$counts, file.path(figures_dir,
         paste0("figureA_consensus_counts_", ver_count, ".tsv")), sep = "\t")
  fwrite(count_result$detail, file.path(figures_dir,
         paste0("figureA_classification_", ver_count, ".tsv")), sep = "\t")
  figA <- ova_count_plot(count_result$counts, tested_targets)
  count_stem <- paste0("figureA_consensus_counts_by_target_", ver_count)
  save_both(file.path(figures_dir, paste0(count_stem, ".pdf")), figA,
            70 / 25.4, 57 / 25.4)
  if (min_count_methods == 2L)
    save_both(file.path(figures_dir, "figureA_consensus_counts_by_target.pdf"), figA,
              70 / 25.4, 57 / 25.4)
}
writeLines(c("Counts: distinct DA genera supported by >=2 (or >=3) methods in the target OVA comparison.",
             "Enriched: unique highest arithmetic mean relative abundance across plotted groups (zeros included).",
             "Depleted: unique lowest mean. Tied extrema and intermediate means are not counted.",
             "Require significant MaAsLin2 and >=2 (or >=3) significant methods with the matching coefficient sign (FDR < 0.05).",
             "Counts include all eligible genera; Figure B shows up to top N enriched candidates per target using identical eligibility."),
           file.path(figures_dir, "figureA_count_definition.txt"))

real_lfc_df <- do.call(rbind, lapply(seq_len(nrow(comps_all)), function(i) {
  tgt <- comps_all$target[i]
  other_groups <- setdiff(valid_sites, tgt)
  idx <- which(site_vec_all %in% c(tgt, other_groups))
  if (length(idx) == 0) return(NULL)
  grp <- ifelse(site_vec_all[idx] == tgt, tgt, "Other")
  grp <- factor(grp, levels = c("Other", tgt))
  lfc <- real_lfc_matrix(sdf_all[idx, , drop = FALSE], grp, ref = "Other", tgt = tgt)
  data.frame(feature = names(lfc), target = tgt, real_lfc = as.numeric(lfc))
}))

sym_limits <- function(x) { s <- max(abs(x), na.rm = TRUE); c(-s, s) }
BUBBLE_SIZE_CAP <- 50

abundance_plots <- function(feats, min_m, real_lfc_df, hm_mode = "zscore") {
  rel <- read_relabund_rows(feats)
  if (nrow(rel) < length(feats)) {
    warning("Only ", nrow(rel), "/", length(feats), " features found in the relabund table")
  }
  samples <- intersect(colnames(rel), rownames(meta))
  sdf <- as.data.frame(t(rel[, samples]))
  colnames(sdf) <- rel[[1]]
  sdf$group_lab <- meta[rownames(sdf), "group_lab"]

  group_n <- sdf %>% group_by(group_lab) %>% summarise(n_total = n(), .groups = "drop")

  mean_by_group <- sdf %>%
    pivot_longer(-group_lab, names_to = "feature", values_to = "abundance") %>%
    group_by(group_lab, feature) %>%
    summarise(mean_abundance = mean(abundance),
              n_present = sum(abundance > 0), .groups = "drop") %>%
    left_join(group_n, by = "group_lab")

  mean_by_group <- mean_by_group %>%
    group_by(feature) %>%
    mutate(z_abun = { x <- mean_abundance; s <- sd(x);
      if (is.finite(s) && s > 0) (x - mean(x)) / s else NA_real_ },
      raw_abun = ifelse(mean_abundance > 0, log10(mean_abundance), NA_real_)) %>%
    ungroup()

  fill_col  <- if (hm_mode == "raw") "raw_abun" else "z_abun"
  fill_name <- if (hm_mode == "raw") "mean relative abundance\n(log10)" else
                                       "mean relative abundance\n(z-score)"
  fill_lims <- if (hm_mode == "raw") {
    r <- mean_by_group$raw_abun
    c(floor(min(r, na.rm = TRUE)), ceiling(max(r, na.rm = TRUE)))
  } else sym_limits(mean_by_group$z_abun)

  feats_keep <- feats[feats %in% rel[[1]]]
  lab_map <- setNames(enh$genus_lab[match(feats_keep, enh$feature)], feats_keep)
  row_levels <- rev(unname(lab_map))
  block_of <- attr(feats, "block")
  row_y <- setNames(seq_along(row_levels), row_levels)
  if (select_mode == "target_topn" && length(block_of) > 0) {
    bottom_blocks <- rev(block_of[match(feats_keep, feats)])
    bottom_block_index <- cumsum(c(TRUE, tail(bottom_blocks, -1L) !=
                                    head(bottom_blocks, -1L)))
    row_y <- row_y + (bottom_block_index - 1L) * OVA_BLOCK_GAP
  }
  y_limits <- range(row_y) + c(-0.5, 0.5)

  mean_by_group$genus_lab <- factor(lab_map[mean_by_group$feature], levels = row_levels)
  mean_by_group$group_lab <- factor(mean_by_group$group_lab, levels = group_lab_order)
  mean_by_group$y_plot <- unname(row_y[as.character(mean_by_group$genus_lab)])
  mean_by_group$x_plot <- as.numeric(mean_by_group$group_lab)

  mean_heatmap <- ggplot(mean_by_group,
                         aes(x = x_plot, y = y_plot, fill = .data[[fill_col]])) +
    geom_tile(color = NA, width = 0.9, height = 1) +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = fill_lims,
                         na.value = "grey85", name = fill_name) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8, face = "italic", family = OVA_GENUS_FONT, hjust = 1),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1, face = "bold"),
      axis.title = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right", legend.direction = "vertical",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.spacing = unit(0, "cm"), legend.box.spacing = unit(0, "cm"),
      legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.25, "cm"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
      plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2)
    ) +
    scale_x_continuous(breaks = seq_along(group_lab_order), labels = group_lab_order,
                       limits = c(0.55, length(group_lab_order) + 0.45),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = unname(row_y), labels = row_levels,
                       limits = y_limits, expand = c(0, 0)) +
    coord_cartesian(clip = "off")

  # Single prevalence panel, matching the age OVA figure: small seeded
  # vertical jitter and transparency reduce overlap without changing x values.
  prev_df <- mean_by_group %>%
    filter(n_total > 0) %>%
    mutate(prevalence = 100 * n_present / n_total,
           genus_lab = factor(genus_lab, levels = row_levels),
           group_lab = factor(group_lab, levels = levels(group_lab))) %>%
    arrange(genus_lab, group_lab)
  set.seed(42)
  prev_df$y_jit <- prev_df$y_plot +
    runif(nrow(prev_df), -0.15, 0.15)

  prevalence_plot <- ggplot(prev_df, aes(x = prevalence, y = y_jit)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45",
               linewidth = 0.4) +
    geom_point(aes(fill = group_lab), shape = 21, size = 1.6, stroke = 0.15,
               alpha = 0.6) +
    labs(x = "Prevalence (%)", y = NULL) +
    scale_fill_manual(values = site_colors, name = "Site",
                      breaks = levels(prev_df$group_lab),
                      guide = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    scale_x_continuous(limits = c(-8, 108), breaks = c(0, 50, 100),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits,
                       breaks = unname(row_y), labels = row_levels,
                       expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8, face = "italic", family = OVA_GENUS_FONT, color = "black", hjust = 1),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 8, face = "bold", margin = margin(t = 3)),
      axis.text.x = element_text(size = 6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
      legend.position = "right", legend.direction = "vertical",
      legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8),
      legend.spacing = unit(0, "cm"), legend.box.spacing = unit(0, "cm"),
      legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.25, "cm"),
      plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2)
    ) + coord_cartesian(clip = "off")

  # DA bubble
  sub <- enh[enh$feature %in% feats, , drop = FALSE]
  diff_bubble <- NULL
  if (nrow(sub) > 0) {
    sub <- sub[!is.na(sub$target), , drop = FALSE]
    if (nrow(sub) > 0) {
      sub$maaslin2_coef <- sub$maaslin2_lfc
      sub$pair_lab <- factor(short_target(sub$target),
                             levels = short_target(group_order))
      sub$genus_lab <- factor(enh$genus_lab[match(sub$feature, enh$feature)],
                              levels = row_levels)
      sub$y_plot <- unname(row_y[as.character(sub$genus_lab)])
      sig_mat <- as.matrix(sub[, c("in_linda", "in_maaslin2", "in_ancombc2")])
      sub$maa_q <- sub$maaslin2_padj
      sub <- sub[is.finite(sub$maa_q) & is.finite(sub$maaslin2_coef) & rowSums(sig_mat) >= min_m, , drop = FALSE]
      if (nrow(sub) > 0) {
        sub$neg_log10_maaq <- -log10(pmax(sub$maa_q, .Machine$double.xmin))
        diff_bubble <- ggplot(sub, aes(x = pair_lab, y = y_plot,
                                       size = neg_log10_maaq, fill = maaslin2_coef)) +
          geom_point(shape = 21, alpha = 0.9, color = "black", stroke = 0.4) +
          scale_size_continuous(range = c(0.8, 3.2), name = "-log10(FDR)") +
          scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                               midpoint = 0, na.value = "grey80",
                               name = "MaAsLin2 coefficient\n(Target vs. Rest)") +
          scale_x_discrete(expand = expansion(add = 0.55)) +
          scale_y_continuous(breaks = unname(row_y), labels = row_levels,
                             limits = y_limits, expand = c(0, 0)) +
          labs(x = NULL, y = NULL) +
          theme_minimal() +
          theme(
            axis.text.y = element_text(size = 8, face = "italic", family = OVA_GENUS_FONT, color = "black", hjust = 1),
            axis.text.x = element_text(size = 8, angle = 45, hjust = 1, face = "bold"),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
            legend.position = "right", legend.direction = "vertical",
            legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8),
            legend.spacing = unit(0, "cm"), legend.box.spacing = unit(0, "cm"),
            legend.key.width = unit(0.2, "cm"), legend.key.height = unit(0.25, "cm"),
            legend.box = "vertical",
            plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2)
          ) + coord_cartesian(clip = "on")
      }
    }
  }

  prev_panel <- prevalence_plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  # Keep the title next to its own ticks despite longer labels in other panels.
  prev_panel <- patchwork::free(prev_panel, type = "label", side = "b")
  combined <- if (is.null(diff_bubble))
    (mean_heatmap | prev_panel) + patchwork::plot_layout(widths = c(6, 3), guides = "collect")
  else
    (diff_bubble | mean_heatmap + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) |
       prev_panel) + patchwork::plot_layout(widths = c(6, 6, 3), guides = "collect")
  list(heatmap = mean_heatmap, prevalence_plot = prevalence_plot,
       diff_bubble = diff_bubble, combined = combined,
       heatmap_data = mean_by_group, bubble_data = sub,
       prevalence_data = prev_df, row_y = row_y)
}

# ==========================================================================
# 4. Figure C: top-5 consensus genera ridgelines per site
# ==========================================================================
top5_abundance_plot <- function(feats) {
  rel <- read_relabund_rows(feats)
  samples <- intersect(colnames(rel), rownames(meta))
  taxa_means <- data.frame(feature = rel[[1]],
                           mean_abundance = rowMeans(rel[, samples])) %>%
    arrange(desc(mean_abundance))
  top5 <- head(taxa_means$feature, 5)
  lab_map <- setNames(enh$genus_lab[match(taxa_means$feature, enh$feature)], taxa_means$feature)

  sdf <- as.data.frame(t(rel[, samples]))
  colnames(sdf) <- rel[[1]]
  sdf$group_lab <- meta[rownames(sdf), "group_lab"]

  groups_present <- group_lab_order[group_lab_order %in% unique(sdf$group_lab)]
  genus_long <- sdf %>%
    pivot_longer(-group_lab, names_to = "feature", values_to = "abundance") %>%
    filter(feature %in% top5) %>%
    mutate(abundance = ifelse(abundance == 0, 1e-10, abundance),
           genus_lab = factor(lab_map[feature], levels = lab_map[top5]),
           group_lab = factor(group_lab, levels = rev(groups_present)))
  grp_cols <- unname(site_colors[groups_present])

  ggplot(genus_long) +
    geom_density_ridges(aes(x = log10(abundance), y = group_lab, fill = group_lab),
                        quantile_lines = TRUE, quantiles = 2, alpha = 0.7) +
    scale_fill_manual(values = setNames(grp_cols, groups_present)) +
    facet_grid(~ genus_lab) +
    labs(x = "log10(Relative Abundance)", y = NULL) +
    theme_minimal() +
    theme(
      legend.position = "none", axis.text.y = element_text(size = 8, face = "bold"),
      axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 14, face = "bold"),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold.italic", family = OVA_GENUS_FONT),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
      plot.margin = margin(10, 10, 10, 10)
    )
}

for (ver in "2method") {
  min_m <- if (ver == "2method") 2L else 3L
  if (select_mode == "recurrence") {
    feats <- select_recurrence_features(enh, min_m, min_targets, top_n)
    cat("    ", ver, ": ", length(feats), " genera robust in >= ", min_targets,
        " targets\n", sep = "")
  } else if (select_mode == "target_topn") {
    feats <- select_target_topn(enh, min_m, mean_abundance_df, each_n,
                                group_levels = group_lab_order, alpha = DA_ALPHA)
    cat("    ", ver, ": ", length(feats), " target-enriched genera (top ",
        each_n, " per global mean-abundance peak; selected by Maaslin2 coef)\n",
        sep = "")
    selection <- attr(feats, "selection")
    fwrite(selection, file.path(figures_dir, paste0("figureB_selected_", ver, ".tsv")),
           sep = "\t")
    fwrite(attr(feats, "candidates"),
           file.path(figures_dir, paste0("figureB_candidates_", ver, ".tsv")), sep = "\t")
    selection_counts <- data.frame(target = group_order) %>%
      left_join(attr(feats, "candidates") %>% filter(passes_peak_filter) %>%
                  count(target, name = "eligible_genera"), by = "target") %>%
      left_join(selection %>% count(target, name = "selected_genera"), by = "target") %>%
      mutate(across(c(eligible_genera, selected_genera), ~replace_na(.x, 0L)))
    fwrite(selection_counts,
           file.path(figures_dir, paste0("figureB_selection_counts_", ver, ".tsv")), sep = "\t")
    checked_counts <- fread(file.path(figures_dir,
      paste0("figureA_consensus_counts_", ver, ".tsv")), data.table = FALSE)
    enriched_counts <- checked_counts %>% filter(direction == "Enriched")
    stopifnot(all(selection_counts$eligible_genera ==
      enriched_counts$n_genera[match(selection_counts$target, enriched_counts$target)] |
      (!selection_counts$target %in% tested_targets & selection_counts$eligible_genera == 0)))
    print(table(factor(selection$target, levels = group_order)))
  } else {
    feats <- select_heatmap_features(enh, min_m, real_lfc_df)
    cat("    ", ver, ": ", length(feats), " genera in the 5+5 union\n", sep = "")
  }
  if (length(feats) == 0) next
  adaptive <- select_mode != "block55"
  gap_height <- if (select_mode == "target_topn")
    max(0, length(unique(attr(feats, "block"))) - 1L) * OVA_BLOCK_GAP else 0
  display_rows <- length(feats) + gap_height
  h_hm <- if (adaptive) max(6, 0.13 * display_rows + 1.8) else 8
  h_cmb <- if (adaptive) max(5.3333, 0.13 * display_rows + 1.5) else 5.3333
  ab_zscore <- abundance_plots(feats, min_m, real_lfc_df, hm_mode = "zscore")
  ab_raw    <- abundance_plots(feats, min_m, real_lfc_df, hm_mode = "raw")
  fwrite(ab_zscore$heatmap_data,
         file.path(figures_dir, paste0("figureB_heatmap_data_", ver, ".tsv")), sep = "\t")
  fwrite(ab_zscore$bubble_data,
         file.path(figures_dir, paste0("figureB_bubble_data_", ver, ".tsv")), sep = "\t")
  fwrite(ab_zscore$prevalence_data,
         file.path(figures_dir, paste0("figureB_prevalence_data_", ver, ".tsv")), sep = "\t")
  save_both(file.path(figures_dir, paste0("figureB_abundance_heatmap_", ver, "_zscore.pdf")),
            ab_zscore$heatmap, 4.5, h_hm)
  save_both(file.path(figures_dir, paste0("figureB_abundance_heatmap_", ver, "_raw.pdf")),
            ab_raw$heatmap, 4.5, h_hm)
  save_both(file.path(figures_dir, paste0("figureB_abundance_heatmap_", ver, ".pdf")),
            ab_zscore$heatmap, 4.5, h_hm)
  save_both(file.path(figures_dir, paste0("figureB_prevalence_", ver, ".pdf")),
            ab_zscore$prevalence_plot, 6, h_hm)
  if (!is.null(ab_zscore$diff_bubble)) {
    save_both(file.path(figures_dir, paste0("figureB_da_bubble_", ver, ".pdf")),
              ab_zscore$diff_bubble, 8, h_hm)
  }
  save_both(file.path(figures_dir, paste0("figureB_abundance_", ver, "_combined_zscore.pdf")),
            ab_zscore$combined, 7.7, h_cmb)
  save_both(file.path(figures_dir, paste0("figureB_abundance_", ver, "_combined_raw.pdf")),
            ab_raw$combined, 7.7, h_cmb)
  save_both(file.path(figures_dir, paste0("figureB_abundance_", ver, "_combined.pdf")),
            ab_zscore$combined, 7.7, h_cmb)
  top5 <- top5_abundance_plot(feats)
  save_both(file.path(figures_dir, paste0("figureC_top5_abundance_", ver, ".pdf")),
            top5, 7, 6.6667)
}

cat("==========================================\n")
if (select_mode == "target_topn") {
  source("differential_abundance/ova_publication_helpers.R")
  save_ova_publication_figures(figures_dir, group_order, "Site", versions = "2method")
}
cat("Done. Figures in:", figures_dir, "\n")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")

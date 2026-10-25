# ==============================================================================
# 03-plot-saliva-ova-bubble.R  (age OVA da_v2; saliva-only)
# Consensus DA figures for the saliva-only OVA age pipeline. Reads the
# consensus tables produced by stage 02 (no re-run of the DA pipeline).
# Comparisons are OVA: each valid age group (target) vs the pooled "Other"
# of the other valid age groups, name "<site>_enrich_<target>", so
#
#   Figure A: target bubble grid of consensus counts
#             (size = #genera; two panels: >=2 methods robust / 3 of 3)
#   Figure B: DA bubble (rows = genera, cols = age targets;
#             size = -log10(Maaslin2 q), color = Maaslin2 coefficient
#             (adjusted log2-scale contrast, target vs Other)),
#             abundance heatmap (cols = saliva x age group) + prevalence
#             panel, combined
#   Figure C: top-5 consensus genera, log10(relabund) ridgelines per
#             saliva x age group
#
# Input:  <out>/consensus_all.tsv / consensus_all_enhanced.tsv
#         data/metadata-and-biom-v2-qc/genus_table_relabund.tsv
# Default directory: <out>/figures_saliva_target_top10_globalpeak_maaslin2coef/
# Legacy --select=block55 output:
#         <out>/figures_saliva/figureA_consensus_counts_saliva_by_target.{pdf,png}
#         <out>/figures_saliva/figureB_da_bubble_saliva_{2method,3method}.{pdf,png}
#         <out>/figures_saliva/figureB_abundance_heatmap_saliva_{2method,3method}.{pdf,png}
#         <out>/figures_saliva/figureB_prevalence_saliva_{2method,3method}.{pdf,png}
#         <out>/figures_saliva/figureB_abundance_saliva_{2method,3method}_combined.{pdf,png}
#         <out>/figures_saliva/figureC_top5_abundance_saliva_{2method,3method}.{pdf,png}
#
# Usage (from project root):
#   Rscript differential_abundance/age_differential_ova/03-plot-saliva-ova-bubble.R
#   Rscript differential_abundance/age_differential_ova/03-plot-saliva-ova-bubble.R --scope=hc
#   DA_RESULTS_ROOT=differential_abundance/result Rscript differential_abundance/age_differential_ova/03-plot-saliva-ova-bubble.R --scope=hc --select=target_topn --each-n=10
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(data.table)
  library(patchwork); library(viridis); library(ggridges); library(grid)
})

source("differential_abundance/age_differential_ova/da_helpers.R")   # rd(), SITE_ORDER, parse_scope_arg, parse_ova_comparison, valid_age_groups
source("differential_abundance/real_lfc_helper.R")               # real_lfc_matrix(), REAL_LFC_PS

source("differential_abundance/ova_plot_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
scope <- parse_scope_arg(args)

arg_val <- function(args, prefix, default) {
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) default else sub(prefix, "", hit[length(hit)])
}
# --select=block55 (legacy: per-comparison 5+5)
#           | recurrence: rows = genera robust in >= min_targets distinct
#             age targets (ranked by #targets, then best Maaslin2 q)
#           | target_topn (default): require a unique mean-abundance peak across all age
#             groups and significant positive Maaslin2 + other-method support;
#             select by Maaslin2 coef; display by group mean abundance
select_mode <- arg_val(args, "--select=", "target_topn")
min_targets <- as.integer(arg_val(args, "--min-targets=", "3"))
top_n       <- as.integer(arg_val(args, "--top-n=", "-1"))
each_n      <- as.integer(arg_val(args, "--each-n=", "10"))
# --heatmap=zscore (default) | raw  : heatmap fill. zscore = per-genus z-score
# of mean relabund across age groups (row-standardized, NA where the genus is
# absent in a group). raw = mean relabund on log10 scale (the genus's absolute
# abundance; also NA where absent). Both versions are emitted every run so the
# two can be compared side by side.
heatmap_mode <- arg_val(args, "--heatmap=", "zscore")
combined_width <- 7.7 * 3 / 4  # inches; retain the original height and text sizes

# The OVA pipeline is saliva-only; the site filter is fixed.
site_filter <- "saliva"

root <- Sys.getenv("DA_RESULTS_ROOT"); if (!nzchar(root)) root <- "differential_abundance/result"
out_dir <- file.path(age_root(scope), "age_da_consensus")
scope_tag <- if (scope == "hc") "age OVA (HC only)" else "age OVA (all)"
site_tag <- site_filter   # used in filenames + titles

group_order  <- SITE_ORDER
age_levels   <- AGE_ORDER
# Keep metadata group keys intact; use an ASCII > for reliable figure labels.
age_display  <- setNames(c("0-6 months", "6-36 months", "3-12 years",
                          "12-18 years", "18-60 years", ">60 year"), age_levels)
age_short    <- c("0-6m", "6-36m", "3-12y", "12-18y", "18-60y", ">60y")
names(age_short) <- age_levels
# Extra blank space between age blocks, measured in genus-row heights.
AGE_BLOCK_GAP <- OVA_BLOCK_GAP
mode_tag <- switch(select_mode,
  recurrence = paste0("_recurrence_min", min_targets,
                      if (top_n > 0) paste0("_top", top_n) else ""),
  target_topn = paste0("_target_top", each_n, "_globalpeak_maaslin2coef"),
  "")
figures_dir <- file.path(out_dir, paste0("figures_", site_tag, mode_tag))
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# OVA comparison "<site>_enrich_<target>" -> column label "<age_short> vs rest"
short_target <- function(g) paste0(age_short[g], " vs rest")

save_both <- function(path_pdf, plot, width, height, dpi = 300) {
  ggsave(path_pdf, plot, device = grDevices::cairo_pdf,
         width = width, height = height, limitsize = FALSE)
  png_path <- sub("\\.pdf$", ".png", path_pdf)
  ggsave(png_path, plot, device = ragg::agg_png,
         width = width, height = height, dpi = dpi, bg = "white",
         limitsize = FALSE)
}

cat("==========================================\n")
cat("Consensus DA plots - age OVA (each age group vs Other, saliva only)\n")
cat("Scope: ", scope_tag, "\n", sep = "")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ==========================================================================
# 1. Load consensus tables
# ==========================================================================
cons_file <- file.path(out_dir, "consensus_all.tsv")
if (!file.exists(cons_file)) stop("Missing: ", cons_file)
cons <- fread(cons_file, data.table = FALSE)
cat("    Consensus rows (all sites):", nrow(cons), "\n")

enh_file <- file.path(out_dir, "consensus_all_enhanced.tsv")
if (!file.exists(enh_file)) stop("Missing: ", enh_file)
enh <- fread(enh_file, data.table = FALSE)
cat("    Enhanced rows (all sites):", nrow(enh), "\n")

# OVA comparison -> site + target (no second age group; "Other" is implicit)
parsed <- lapply(enh$comparison, parse_ova_comparison)
enh$site      <- vapply(parsed, `[[`, character(1), "site")
enh$target    <- vapply(parsed, `[[`, character(1), "target")
enh$pair_lab  <- ifelse(!is.na(enh$site) & !is.na(enh$target),
                        short_target(enh$target), enh$comparison)

parsed_c <- lapply(cons$comparison, parse_ova_comparison)
cons$site      <- vapply(parsed_c, `[[`, character(1), "site")
cons$target    <- vapply(parsed_c, `[[`, character(1), "target")
cons$pair_lab  <- ifelse(!is.na(cons$site) & !is.na(cons$target),
                         short_target(cons$target), cons$comparison)

# ---- subset to saliva only (the only OVA site) ---------------------------
cons <- cons[!is.na(cons$site) & cons$site == site_filter, , drop = FALSE]
enh <- enh[!is.na(enh$site) & enh$site == site_filter, , drop = FALSE]
cat("    Consensus rows (", site_tag, "):", nrow(cons), "\n", sep = "")
cat("    Enhanced rows (", site_tag, "):", nrow(enh), "\n", sep = "")
if (nrow(cons) == 0) stop("No consensus rows for site = ", site_filter,
                          "; check --site / --scope.")

# ==========================================================================
# 2. Figure A: target bubble grid of consensus counts (saliva, one row)
# ==========================================================================
pair_counts <- cons %>%
  filter(!is.na(site), !is.na(target)) %>%
  group_by(site, target) %>%
  summarise(n2 = sum(n_methods >= 2),
            n3 = sum(n_methods == 3), .groups = "drop") %>%
  mutate(pair_lab = factor(short_target(target),
                           levels = short_target(age_levels)))

sites <- group_order[group_order %in% pair_counts$site]
pair_counts$site <- factor(pair_counts$site, levels = sites)
pair_counts_long <- pair_counts %>%
  select(site, pair_lab, n2, n3) %>%
  pivot_longer(c(n2, n3), names_to = "ver", values_to = "n") %>%
  mutate(ver = recode(ver, n2 = ">=2 methods (robust)", n3 = "3/3 methods"))

figA <- ggplot(pair_counts_long,
               aes(x = pair_lab, y = site, size = n, color = n)) +
  geom_point() +
  geom_text(aes(label = n), size = 3, color = "black", vjust = -1.4, show.legend = FALSE) +
  scale_size_continuous(range = c(1.5, 8), name = "# genera") +
  scale_color_viridis_c(option = "magma", direction = -1, name = "# genera") +
  facet_wrap(~ ver, ncol = 1, scales = "free_x") +
  labs(x = NULL, y = NULL,
       title = paste0("Age DA OVA consensus per target age group - ",
                      scope_tag, " (", site_tag, ")")) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
save_both(file.path(figures_dir,
                    paste0("figureA_consensus_counts_", site_tag, "_by_target.pdf")),
          figA, 8, 6)

cat("    Sites:", paste(sites, collapse = ", "), "\n")
cat("    Targets with any consensus:", nrow(pair_counts), "\n\n")

# ==========================================================================
# 3. Figure B: DA bubble + abundance heatmap + barplot  (saliva only)
#    Rows = consensus genera; bubble columns = age pairs (saliva only);
#    heatmap columns = saliva x age group (ordered by AGE_ORDER).
# ==========================================================================

# genus label: genus_short, fallback to last non-empty taxonomy level
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

# features per version: per-comparison "5+5" selection. Candidates are the
# robust genera (>= min_methods methods with padj < 0.05 = the q filter),
# then within each direction (real relabund LFC > 0 = TARGET-enriched,
# < 0 = Other-enriched) the 5 with the LARGEST |real LFC| are picked (top
# effect size instead of top q, so near-noise high-significance genera do
# not fill the slots). Union across comparisons with dedupe (first-come:
# a genus is kept only in the first block that claims it). The returned
# order IS the figure row order: blocks follow the comparison columns
# left-to-right (age order), target-enriched genera first within each block.
select_heatmap_features <- function(enh, min_methods, real_lfc_df,
                                    n_each = 5L) {
  keep <- enh[enh$n_methods >= min_methods & !is.na(enh$maaslin2_padj) &
              !is.na(enh$maaslin2_lfc), , drop = FALSE]
  keep$real_lfc <- real_lfc_df$real_lfc[match(
    paste(keep$feature, keep$site, keep$target),
    paste(real_lfc_df$feature, real_lfc_df$site, real_lfc_df$target))]
  keep <- keep[!is.na(keep$real_lfc), , drop = FALSE]

  comps <- unique(keep[, c("site", "target")])
  comps$tgt_f <- factor(comps$target, levels = age_levels)
  comps <- comps[order(comps$tgt_f), , drop = FALSE]

  ordered_rows <- character(0)
  block_of <- character(0)
  for (i in seq_len(nrow(comps))) {
    cc <- keep[keep$site == comps$site[i] & keep$target == comps$target[i],
               , drop = FALSE]
    up <- cc[which(cc$real_lfc > 0), , drop = FALSE]
    dn <- cc[which(cc$real_lfc < 0), , drop = FALSE]
    # largest |real LFC| first (top effect size, not top q)
    up <- up[order(-up$real_lfc), , drop = FALSE]
    dn <- dn[order(dn$real_lfc), , drop = FALSE]
    take <- c(head(up$feature, n_each), head(dn$feature, n_each))
    take <- take[!take %in% ordered_rows]
    if (length(take) > 0) {
      # claim order = figure row order; the block label is kept as an
      # attribute so the bubble panel can draw dashed block separators
      ordered_rows <- c(ordered_rows, take)
      block_of <- c(block_of, rep(comps$target[i], length(take)))
    }
  }
  attr(ordered_rows, "block") <- block_of
  ordered_rows
}

# Alternate row selection (draft): recurrence across targets. A genus is kept
# when it is robust (>= min_methods methods with padj < 0.05) in AT LEAST
# min_targets distinct OVA age targets (each target = "this age vs the pooled
# Other"). Ranked by #targets (desc), then best (min) Maaslin2 q across those
# targets; optional top_n cap. NOTE: multi-target significance is NOT
# independent replication (the "Other" pools overlap across targets); a genus
# significant in several targets with sign reversal reads as a non-monotonic
# (peak-shaped) age trajectory rather than a replicated same-direction effect.
select_recurrence_features <- function(enh, min_methods, min_targets,
                                       top_n = -1L) {
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

# ---- load metadata + relative abundance (only the needed feature rows) ----
meta <- load_age_metadata(hc_only = (scope == "hc"))
meta <- meta[meta[[COL_SITE]] %in% group_order &
             meta[[COL_AGE]] %in% age_levels, , drop = FALSE]
# restrict metadata to the selected site so heatmap columns are saliva-only
meta <- meta[meta[[COL_SITE]] == site_filter, , drop = FALSE]
meta$group_lab <- paste(meta[[COL_SITE]], meta[[COL_AGE]])
group_lab_order <- as.vector(t(outer(sites, age_levels, paste)))

# read only the selected feature rows from the ~2.4 GB relabund tsv
# (grep -F with a trailing tab anchors each pattern to the feature column)
read_relabund_rows <- function(feats) {
  pats <- paste0(feats, "\t")
  cmd <- paste0("head -1 ", RELABUND_TSV, "; grep -F -e '",
                paste(pats, collapse = "' -e '"), "' ", RELABUND_TSV)
  fread(cmd = cmd, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
}

# age-ordered palette (earlier -> older along viridis)
age_colors <- setNames(viridis(length(age_levels)), age_levels)

# ---- real relabund LFC for every robust feature (top-level, shared by the
# 5+5 selection; the DA bubble uses Maaslin2 coefficient) ----
# The "Other" pool is rebuilt with the SAME validity rule as the DA driver
# (valid_age_groups = age groups with >= MIN_SAMPLES_PER_GROUP samples in the
# site), so the log2 group-mean ratios are computed on exactly the samples
# each model saw.
valid_groups <- valid_age_groups(meta, site_filter, MIN_SAMPLES_PER_GROUP)
robust_feats <- unique(enh$feature[enh$n_methods >= 2 &
                                   !is.na(enh$maaslin2_padj)])
rel_all <- read_relabund_rows(robust_feats)
samples_all <- intersect(colnames(rel_all), rownames(meta))
sdf_all <- as.data.frame(t(rel_all[, samples_all]))
colnames(sdf_all) <- rel_all[[1]]
site_vec_all <- meta[rownames(sdf_all), COL_SITE]
age_vec_all  <- meta[rownames(sdf_all), COL_AGE]
comps_all <- unique(enh[, c("site", "target")])
comps_all <- comps_all[!is.na(comps_all$target), , drop = FALSE]
comps_all$tgt_f <- factor(comps_all$target, levels = age_levels)
comps_all <- comps_all[order(comps_all$tgt_f), , drop = FALSE]
# The same HC saliva samples used for the heatmap, including zero abundances.
mean_abundance_df <- do.call(rbind, lapply(age_levels, function(tgt) {
  idx <- which(site_vec_all == site_filter & age_vec_all == tgt)
  if (length(idx) == 0L) stop("No abundance samples for age group: ", tgt)
  means <- colMeans(sdf_all[idx, , drop = FALSE])
  data.frame(feature = names(means), site = site_filter, target = tgt,
             target_mean_abundance = as.numeric(means))
}))
real_lfc_df <- do.call(rbind, lapply(seq_len(nrow(comps_all)), function(i) {
  st <- comps_all$site[i]; tgt <- comps_all$target[i]
  other_groups <- setdiff(valid_groups, tgt)
  # samples of the target group + the Other pool (only valid groups),
  # recoded to a two-level factor {Other, target}
  idx <- which(site_vec_all == st & age_vec_all %in% c(tgt, other_groups))
  if (length(idx) == 0) return(NULL)
  grp <- ifelse(age_vec_all[idx] == tgt, tgt, "Other")
  grp <- factor(grp, levels = c("Other", tgt))
  lfc <- real_lfc_matrix(sdf_all[idx, , drop = FALSE],
                         grp, ref = "Other", tgt = tgt)
  data.frame(feature = names(lfc), site = st, target = tgt,
             real_lfc = as.numeric(lfc))
}))

# symmetric fill limits around 0 for the z-score heatmap
sym_limits <- function(x) {
  s <- max(abs(x), na.rm = TRUE)
  c(-s, s)
}

# bubble size: -log10(Maaslin2 q) capped here (see diff-bubble panel)
BUBBLE_SIZE_CAP <- 50

abundance_plots <- function(feats, min_m, real_lfc_df, hm_mode = "zscore") {
  rel <- read_relabund_rows(feats)
  if (nrow(rel) < length(feats)) {
    warning("Only ", nrow(rel), "/", length(feats),
            " features found in the relabund table")
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

  # fill = per-genus z-score of mean relative abundance across the age
  # groups (row-standardized so within-genus age patterns show regardless of
  # the genus's absolute abundance level; >0 = higher than that genus's own
  # cross-age mean). No log transform — z-score already removes the magnitude
  # axis, and an extra log would over-compress the high-abundance tiles and
  # over-amplify the low ones. Tiles where the genus is absent in a group
  # (mean_abundance == 0) are kept as NA and rendered grey by na.value, so
  # "genuinely not detected" reads distinctly from "detected but below the
  # genus's own mean".
  mean_by_group <- mean_by_group %>%
    group_by(feature) %>%
    mutate(z_abun = {
      x <- mean_abundance
      s <- sd(x)
      if (is.finite(s) && s > 0) (x - mean(x)) / s else NA_real_
    },
    raw_abun = ifelse(mean_abundance > 0, log10(mean_abundance), NA_real_)) %>%
    ungroup()

  fill_col   <- if (hm_mode == "raw") "raw_abun" else "z_abun"
  fill_name  <- if (hm_mode == "raw") "mean relative abundance\n(log10)" else
                                       "mean relative abundance\n(z-score)"
  fill_lims  <- if (hm_mode == "raw") {
    r <- mean_by_group$raw_abun
    c(floor(min(r, na.rm = TRUE)), ceiling(max(r, na.rm = TRUE)))
  } else sym_limits(mean_by_group$z_abun)

  taxa_means <- data.frame(
    feature = rel[[1]],
    mean_abundance = rowMeans(rel[, samples])
  ) %>% arrange(desc(mean_abundance))

  # feats is already in display order: target blocks follow age order.
  feats_keep <- feats[feats %in% rel[[1]]]
  lab_map <- setNames(enh$genus_lab[match(feats_keep, enh$feature)],
                      feats_keep)
  row_levels <- rev(unname(lab_map))
  block_of <- attr(feats, "block")
  row_y <- setNames(seq_along(row_levels), row_levels)
  if (select_mode == "target_topn" && length(block_of) > 0) {
    blocks_kept <- block_of[match(feats_keep, feats)]
    bottom_blocks <- rev(blocks_kept)
    bottom_block_index <- cumsum(c(TRUE, tail(bottom_blocks, -1L) !=
                                    head(bottom_blocks, -1L)))
    row_y <- row_y + (bottom_block_index - 1L) * AGE_BLOCK_GAP
  }
  y_limits <- range(row_y) + c(-0.5, 0.5)

  mean_by_group$genus_lab <- factor(lab_map[mean_by_group$feature], levels = row_levels)
  mean_by_group$y_plot <- unname(row_y[as.character(mean_by_group$genus_lab)])
  mean_by_group$group_lab <- factor(mean_by_group$group_lab, levels = group_lab_order)
  # x-axis label drops the redundant site prefix; only age group shown
  mean_by_group$age_lab <- factor(sub("^\\S+ ", "", mean_by_group$group_lab),
                                  levels = sub("^\\S+ ", "", levels(mean_by_group$group_lab)))
  mean_by_group$x_plot <- as.numeric(mean_by_group$age_lab)

  mean_heatmap <- ggplot(mean_by_group,
                         aes(x = x_plot, y = y_plot,
                             fill = .data[[fill_col]])) +
    # Solid cells within blocks; white gutters separate age columns and row blocks.
    geom_tile(color = NA, width = 0.9, height = 1) +
    # viridis with symmetric limits so 0 (the genus's own cross-age mean)
    # sits at the palette midpoint, dark = below own mean, bright = above.
    # na.value = grey for groups where the genus is absent (mean == 0).
    scale_fill_viridis_c(option = "mako", direction = -1,
                         limits = fill_lims,
                         na.value = "grey85",
                         name = fill_name) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8, face = "italic", family = OVA_GENUS_FONT, hjust = 1),
      axis.text.x = element_text(size = 8, angle = 60, hjust = 1, face = "bold"),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.spacing = unit(0, "cm"),
      legend.box.spacing = unit(0, "cm"),
      legend.key.width = unit(0.2, "cm"),
      legend.key.height = unit(0.25, "cm"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
      plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2)
    ) +
    scale_x_continuous(breaks = seq_along(age_levels),
                       labels = unname(age_display[age_levels]),
                       limits = c(0.55, length(age_levels) + 0.45),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = unname(row_y), labels = row_levels,
                       limits = y_limits, expand = c(0, 0)) +
    coord_cartesian(clip = "off")

  # ---- prevalence point panel (right of the combined figure) ---------------
  # For each genus x age group: prevalence = n_present / n_total (from the
  # same per-group counts already used for the heatmap tiles). Single panel
  # (one column, no age facets): every age group's point overlays the same
  # genus row, jittered vertically within the row so per-group values stay
  # distinguishable; translucent points reveal overlap. Fill color = age
  # group (legend on the right, collected by patchwork). Seeded vertical
  # jitter keeps the figure reproducible without changing prevalence values.
  prev_df <- mean_by_group %>%
    filter(n_total > 0) %>%
    mutate(prevalence = 100 * n_present / n_total,
           genus_lab = factor(genus_lab, levels = row_levels),
           age_lab  = factor(age_lab, levels = levels(age_lab))) %>%
    arrange(genus_lab, age_lab)
  set.seed(42)
  prev_df$y_jit <- prev_df$y_plot +
    runif(nrow(prev_df), -0.15, 0.15)
  prev_ranges <- prev_df %>%
    group_by(genus_lab, y_plot) %>%
    summarise(prev_min = min(prevalence), prev_max = max(prevalence),
              .groups = "drop")

  prevalence_plot <- ggplot(prev_df,
                            aes(x = prevalence, y = y_jit)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45",
               linewidth = 0.4) +
    geom_segment(data = prev_ranges,
                  aes(x = prev_min, xend = prev_max, y = y_plot, yend = y_plot),
                  inherit.aes = FALSE, colour = "grey75", linewidth = 0.3) +
    geom_point(aes(fill = age_lab), shape = 21, size = 1.6, stroke = 0.15,
               alpha = 0.6) +
    labs(x = "Prevalence (%)", y = NULL) +
    scale_fill_manual(values = age_colors, labels = age_short,
                      name = "Age group",
                      guide = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    scale_x_continuous(limits = c(-8, 108),
                       breaks = c(0, 50, 100),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits,
                       breaks = unname(row_y), labels = row_levels,
                       expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 8, face = "italic", family = OVA_GENUS_FONT, color = "black", hjust = 1),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 7, face = "bold", margin = margin(t = 3)),
      axis.text.x  = element_text(size = 6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.spacing = unit(0, "cm"),
      legend.box.spacing = unit(0, "cm"),
      legend.key.width = unit(0.2, "cm"),
      legend.key.height = unit(0.25, "cm"),
      plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2)
    ) +
    coord_cartesian(clip = "off")

  # ---- differential-abundance bubble panel (left of the combined figure) ----
  # rows = same genera, columns = age targets, saliva only. Point size =
  # -log10(Maaslin2 q) (q small = large point); color = native Maaslin2 coef
  # (adjusted log2-scale contrast; >0 = positive association with the target).
  sub <- enh[enh$feature %in% feats, , drop = FALSE]
  diff_bubble <- NULL
  if (nrow(sub) > 0) {
    sub <- sub[!is.na(sub$site) & !is.na(sub$target), , drop = FALSE]
    if (nrow(sub) > 0) {
      sub$maaslin2_coef <- sub$maaslin2_lfc

      sub$site     <- factor(sub$site, levels = sites)
      sub$pair_lab <- factor(short_target(sub$target),
                             levels = levels(pair_counts$pair_lab))
      sub$genus_lab <- factor(enh$genus_lab[match(sub$feature, enh$feature)],
                              levels = row_levels)
      sub$y_plot <- unname(row_y[as.character(sub$genus_lab)])

      sig_mat  <- as.matrix(sub[, c("in_linda", "in_maaslin2", "in_ancombc2")])
      # Show comparisons supported by >= min_m significant methods, with
      # finite Maaslin2 coefficient/q available for the two visual mappings.
      sub$maa_q <- sub$maaslin2_padj
      sub <- sub[is.finite(sub$maa_q) & is.finite(sub$maaslin2_coef) &
                 rowSums(sig_mat) >= min_m, , drop = FALSE]

      if (nrow(sub) > 0) {
        sub$neg_log10_maaq <- -log10(pmax(sub$maa_q, .Machine$double.xmin))
        diff_bubble <- ggplot(sub, aes(x = pair_lab, y = y_plot,
                                       size = neg_log10_maaq, fill = maaslin2_coef)) +
          # shape 21: filled circles with a black outline so pale mid-gradient
          # bubbles stay visible on the white background
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
            axis.text.x = element_text(size = 8, angle = 60, hjust = 1, face = "bold"),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
            legend.position = "right",
            legend.direction = "vertical",
            legend.title = element_text(size = 9, face = "bold"),
            legend.text = element_text(size = 8),
            legend.spacing = unit(0, "cm"),
            legend.box.spacing = unit(0, "cm"),
            legend.key.width = unit(0.2, "cm"),
            legend.key.height = unit(0.25, "cm"),
            legend.box = "vertical",
            plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2)
          ) +
          coord_cartesian(clip = "on")   # clip so large bubbles never cross the panel border
        }
      }
    }

  # All panels share row_y/y_limits, including blank spaces between age blocks.
  prev_panel <- prevalence_plot + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.x = element_text(hjust = 0))
  # Keep the title next to its own ticks despite longer labels in other panels.
  prev_panel <- patchwork::free(prev_panel, type = "label", side = "b")
  # widths: bubble : heatmap : prevalence. Legends collected to a single
  # shared area on the right via guides="collect", including the age-group
  # legend for the single prevalence panel.
  combined <- if (is.null(diff_bubble))
    (mean_heatmap | prev_panel) +
    patchwork::plot_layout(widths = c(6, 3), guides = "collect")
  else
    (diff_bubble | mean_heatmap + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank()) |
       prev_panel) +
    patchwork::plot_layout(widths = c(6, 6, 3), guides = "collect")
  list(heatmap = mean_heatmap, prevalence_plot = prevalence_plot,
       diff_bubble = diff_bubble, combined = combined,
       heatmap_data = mean_by_group, bubble_data = sub,
       prevalence_data = prev_df, row_y = row_y)
}

# ==========================================================================
# 4. Top 5 Genera Abundance Distribution (Figure C, density ridgeline)
#    Among the consensus genera (same selection as Figure B), take the top 5
#    by global mean abundance and draw log10(relative abundance) ridgelines
#    per site x age group. Zero abundances are floored at 1e-10 before log10.
# ==========================================================================
top5_abundance_plot <- function(feats) {
  rel <- read_relabund_rows(feats)
  samples <- intersect(colnames(rel), rownames(meta))
  taxa_means <- data.frame(
    feature = rel[[1]],
    mean_abundance = rowMeans(rel[, samples])
  ) %>% arrange(desc(mean_abundance))
  top5 <- head(taxa_means$feature, 5)

  # genus labels for the top-5 features (same mapping as abundance_plots)
  lab_map <- setNames(enh$genus_lab[match(taxa_means$feature, enh$feature)],
                      taxa_means$feature)

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

  # fill color cycles through the age-ordered palette per site
  grp_site <- sub(" .*", "", groups_present)
  grp_age  <- sub("^\\S+ ", "", groups_present)
  grp_cols <- unname(age_colors[grp_age])
  grp_labels <- setNames(paste(grp_site, unname(age_display[grp_age])),
                         groups_present)

  ggplot(genus_long) +
    geom_density_ridges(aes(x = log10(abundance), y = group_lab, fill = group_lab),
                        quantile_lines = TRUE, quantiles = 2, alpha = 0.7) +
    scale_fill_manual(values = setNames(grp_cols, groups_present)) +
    scale_y_discrete(labels = grp_labels) +
    facet_grid(~ genus_lab) +
    labs(x = "log10(Relative Abundance)", y = NULL) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 8, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold.italic", family = OVA_GENUS_FONT),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
      plot.margin = margin(10, 10, 10, 10)
    )
}

for (ver in c("2method", "3method")) {
  min_m <- if (ver == "2method") 2L else 3L
  if (select_mode == "recurrence") {
    feats <- select_recurrence_features(enh, min_m, min_targets, top_n)
    cat("    ", ver, ": ", length(feats), " genera robust in >= ", min_targets,
        " targets\n", sep = "")
  } else if (select_mode == "target_topn") {
    feats <- select_target_topn(enh, min_m, mean_abundance_df, each_n,
                                group_levels = age_levels, alpha = DA_ALPHA)
    cat("    ", ver, ": ", length(feats), " target-enriched genera (top ",
        each_n, " per global mean-abundance peak; selected by Maaslin2 coef)\n",
        sep = "")
    selection <- attr(feats, "selection")
    fwrite(selection, file.path(figures_dir, paste0("figureB_selected_", ver, ".tsv")),
           sep = "\t")
    fwrite(attr(feats, "candidates"),
           file.path(figures_dir, paste0("figureB_candidates_", ver, ".tsv")), sep = "\t")
    selection_counts <- data.frame(target = age_levels) %>%
      left_join(attr(feats, "candidates") %>% filter(passes_peak_filter) %>%
                  count(target, name = "eligible_genera"), by = "target") %>%
      left_join(selection %>% count(target, name = "selected_genera"), by = "target") %>%
      mutate(across(c(eligible_genera, selected_genera), ~replace_na(.x, 0L)))
    fwrite(selection_counts,
           file.path(figures_dir, paste0("figureB_selection_counts_", ver, ".tsv")), sep = "\t")
    print(table(factor(selection$target, levels = age_levels)))
  } else {
    feats <- select_heatmap_features(enh, min_m, real_lfc_df)
    cat("    ", ver, ": ", length(feats), " genera in the 5+5 union\n", sep = "")
  }
  if (length(feats) == 0) next
  # scale the heights with the row count in the alternate modes (more rows
  # than the legacy 5+5 union); block55 keeps the original fixed sizes
  adaptive <- select_mode != "block55"
  gap_height <- if (select_mode == "target_topn")
    max(0, length(unique(attr(feats, "block"))) - 1L) * AGE_BLOCK_GAP else 0
  display_rows <- length(feats) + gap_height
  h_hm <- if (adaptive) max(6, 0.13 * display_rows + 1.8) else 8
  h_prev <- h_hm
  h_cmb <- if (adaptive) max(5.3333, 0.13 * display_rows + 1.5) else 5.3333
  # Both heatmap versions are emitted so z-score and raw mean relabund can be
  # compared side by side. prevalence / DA bubble / Figure C do not depend on
  # the heatmap fill, so they are taken from the zscore run (the default) and
  # saved once without a heatmap-mode suffix to preserve the legacy filenames.
  ab_zscore <- abundance_plots(feats, min_m, real_lfc_df, hm_mode = "zscore")
  ab_raw    <- abundance_plots(feats, min_m, real_lfc_df, hm_mode = "raw")
  fwrite(ab_zscore$heatmap_data,
         file.path(figures_dir, paste0("figureB_heatmap_data_", ver, ".tsv")), sep = "\t")
  fwrite(ab_zscore$bubble_data,
         file.path(figures_dir, paste0("figureB_bubble_data_", ver, ".tsv")), sep = "\t")
  fwrite(ab_zscore$prevalence_data,
         file.path(figures_dir, paste0("figureB_prevalence_data_", ver, ".tsv")), sep = "\t")
  save_both(file.path(figures_dir,
                      paste0("figureB_abundance_heatmap_", site_tag, "_", ver, "_zscore.pdf")),
            ab_zscore$heatmap, 4.5, h_hm)
  save_both(file.path(figures_dir,
                      paste0("figureB_abundance_heatmap_", site_tag, "_", ver, "_raw.pdf")),
            ab_raw$heatmap, 4.5, h_hm)
  # legacy unsuffixed heatmap = zscore (kept for backward compatibility)
  save_both(file.path(figures_dir,
                      paste0("figureB_abundance_heatmap_", site_tag, "_", ver, ".pdf")),
            ab_zscore$heatmap, 4.5, h_hm)
  save_both(file.path(figures_dir,
                      paste0("figureB_prevalence_", site_tag, "_", ver, ".pdf")),
            ab_zscore$prevalence_plot, 6, h_prev)
  if (!is.null(ab_zscore$diff_bubble)) {
    save_both(file.path(figures_dir,
                        paste0("figureB_da_bubble_", site_tag, "_", ver, ".pdf")),
              ab_zscore$diff_bubble, 8, h_hm)
  }
  # combined: zscore (default) and raw, both suffixed; plus a legacy unsuffixed
  # zscore version so old references still resolve.
  save_both(file.path(figures_dir,
                      paste0("figureB_abundance_", site_tag, "_", ver, "_combined_zscore.pdf")),
            ab_zscore$combined, combined_width, h_cmb)
  save_both(file.path(figures_dir,
                      paste0("figureB_abundance_", site_tag, "_", ver, "_combined_raw.pdf")),
            ab_raw$combined, combined_width, h_cmb)
  save_both(file.path(figures_dir,
                      paste0("figureB_abundance_", site_tag, "_", ver, "_combined.pdf")),
            ab_zscore$combined, combined_width, h_cmb)
  # Figure C: top 5 consensus genera by mean abundance — density ridgelines
  top5 <- top5_abundance_plot(feats)
  save_both(file.path(figures_dir,
                      paste0("figureC_top5_abundance_", site_tag, "_", ver, ".pdf")),
            top5, 7, 6.6667)
}

cat("==========================================\n")
cat("Done. Figures in:", figures_dir, "\n")
cat("  figureA_consensus_counts_", site_tag, "_by_target.{pdf,png}\n", sep = "")
cat("  figureB_da_bubble_", site_tag, "{2method,3method}.{pdf,png}\n", sep = "")
cat("  figureB_abundance_heatmap_", site_tag, "{2method,3method}{,_zscore,_raw}.{pdf,png}\n", sep = "")
cat("  figureB_prevalence_", site_tag, "{2method,3method}.{pdf,png}\n", sep = "")
cat("  figureB_abundance_", site_tag, "{2method,3method}_combined{,_zscore,_raw}.{pdf,png}\n", sep = "")
cat("  figureC_top5_abundance_", site_tag, "{2method,3method}.{pdf,png}\n", sep = "")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")

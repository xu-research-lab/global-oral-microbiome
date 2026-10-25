# ==============================================================================
# 03-plot-consensus-bubble.R  (sex-within-continent da_v2; per-continent F_vs_M,
#                             saliva — this directory)
# Consensus DA plots for the sex-within-continent analysis. Comparisons are
# "saliva_<continent>_F_vs_M" (one per valid continent), so Figure A is a
# per-continent consensus-count bar chart. Figures B/C follow the da_v2 style
# (rows = genera, columns = continents or continent x sex groups).
#
#   Figure A: per-continent bar chart of consensus counts (n>=2 / 3 of 3)
#   Figure B: DA bubble (rows = genera, cols = continents; size = -log10(max q),
#             color = real relabund log2FC (group-mean ratio, F = ref)),
#             abundance heatmap (cols = continent x sex) + barplot, combined
#   Figure C: top-5 consensus genera, log10(relabund) ridgelines per
#             continent x sex group
#
# Input:  <out>/consensus_all.tsv / consensus_all_enhanced.tsv
#         data/metadata-and-biom-v2-qc/genus_table_relabund.tsv
# Output: <out>/figures/figureA_consensus_counts_per_continent.{pdf,png}
#         <out>/figures/figureB_da_bubble_{2method,3method}.{pdf,png}
#         <out>/figures/figureB_abundance_heatmap_{2method,3method}.{pdf,png}
#         <out>/figures/figureB_abundance_barplot_{2method,3method}.{pdf,png}
#         <out>/figures/figureB_abundance_{2method,3method}_combined.{pdf,png}
#         <out>/figures/figureC_top5_abundance_{2method,3method}.{pdf,png}
#
# Usage (from project root):
#   Rscript differential_abundance/sex_differential_continent/03-plot-consensus-bubble.R --scope=hc
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(data.table)
  library(patchwork); library(viridis); library(ggridges)
})

source("differential_abundance/sex_differential_continent/da_helpers.R")   # rd(), CONTINENT_ORDER, parse_scope_arg
source("differential_abundance/real_lfc_helper.R")                         # real_lfc_matrix(), REAL_LFC_PS

args <- commandArgs(trailingOnly = TRUE)
scope <- parse_scope_arg(args)

root <- Sys.getenv("DA_RESULTS_ROOT"); if (!nzchar(root)) root <- "differential_abundance/result"
out_dir <- file.path(continent_root(scope), "sex_continent_da_consensus")
scope_tag <- if (scope == "hc") "sex within continent, saliva (HC only)" else "sex within continent, saliva (all)"

group_order  <- CONTINENT_ORDER
sex_colors <- c("F" = "#C51B7D", "M" = "#1F78B4")
figures_dir <- file.path(out_dir, "figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

save_both <- function(path_pdf, plot, width, height, dpi = 300) {
  ggsave(path_pdf, plot, width = width, height = height, limitsize = FALSE)
  png_path <- sub("\\.pdf$", ".png", path_pdf)
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white",
         limitsize = FALSE)
}

cat("==========================================\n")
cat("Consensus DA plots — sex (per-continent F vs M, saliva)\n")
cat("Scope: ", scope_tag, "\n", sep = "")
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

# comparison "saliva_<continent>_F_vs_M" -> continent.
# Strip "_F_vs_M", drop the "saliva_" prefix, then reverse sanitize
# ("North_America" -> "North America") and match directly to CONTINENT_ORDER
# (which uses the display names with spaces, e.g. "North America").
parse_continent <- function(comp) {
  s <- sub("_F_vs_M$", "", comp)            # "saliva_North_America"
  s <- sub("^saliva_", "", s)               # "North_America"
  gsub("_", " ", s)                         # "North America"
}
cons$continent <- parse_continent(cons$comparison)
enh$continent  <- parse_continent(enh$comparison)

# ==========================================================================
# 2. Figure A: consensus counts per continent
# ==========================================================================
pair_counts <- cons %>%
  group_by(continent) %>%
  summarise(n2 = sum(n_methods >= 2),
            n3 = sum(n_methods == 3), .groups = "drop")
continents <- group_order[group_order %in% pair_counts$continent]
pair_counts$continent <- factor(pair_counts$continent, levels = continents)
max2 <- max(pair_counts$n2, 0); max3 <- max(pair_counts$n3, 0)
ymax <- max(max2, max3, 1)

figA <- ggplot(pair_counts) +
  geom_col(aes(x = continent, y = n2, fill = ">=2 methods (robust)"),
           width = 0.7, color = "grey25", linewidth = 0.3) +
  geom_col(aes(x = continent, y = n3, fill = "3/3 methods"),
           width = 0.45, color = "grey15", linewidth = 0.3) +
  geom_text(aes(x = continent, y = n2, label = n2), vjust = -0.5, size = 4,
            fontface = "bold") +
  geom_text(aes(x = continent, y = n3, label = n3), vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = c(">=2 methods (robust)" = "grey70",
                               "3/3 methods" = "grey25"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "Number of differently abundant genera",
       title = paste0("Sex DA consensus (F vs M) per continent (saliva) — ", scope_tag)) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1, face = "bold"),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )
save_both(file.path(figures_dir, "figureA_consensus_counts_per_continent.pdf"),
          figA, 9, 6)

cat("    Continents:", paste(continents, collapse = ", "), "\n")
cat("    Counts (>=2 methods):", paste(pair_counts$n2, collapse = ", "), "\n")
cat("    Counts (3/3 methods):", paste(pair_counts$n3, collapse = ", "), "\n\n")

# ==========================================================================
# 3. Figure B: DA bubble + abundance heatmap + barplot
#    Rows = consensus genera (same selection as the reference da_v2 plots),
#    columns = sites (bubble) or site x sex (heatmap).
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

# consensus LFC = median over the significant methods (used for row ordering)
lfc_cols <- c("linda_lfc", "maaslin2_lfc", "ancombc2_lfc")
sig_cols <- c("in_linda", "in_maaslin2", "in_ancombc2")
lfc_mat <- as.matrix(enh[, lfc_cols]); lfc_mat[!as.matrix(enh[, sig_cols])] <- NA
enh$consensus_lfc <- apply(lfc_mat, 1, function(x)
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE))

# features per version: significant in >=2 comparisons, top 40 by
# (n_comparisons, max |consensus LFC|); fallback top 20 by max |LFC|
select_heatmap_features <- function(enh, min_methods) {
  keep <- enh[enh$n_methods >= min_methods, , drop = FALSE]
  feat_n <- keep %>%
    group_by(feature) %>%
    summarise(n_comp = n_distinct(comparison),
              m = max(abs(consensus_lfc), na.rm = TRUE), .groups = "drop")
  shared <- feat_n %>% filter(n_comp >= 2) %>% pull(feature)
  if (length(shared) == 0) {
    return(feat_n %>% arrange(desc(m)) %>% head(20) %>% pull(feature))
  }
  feat_n %>% filter(feature %in% shared) %>%
    arrange(desc(n_comp), desc(m)) %>% head(40) %>% pull(feature)
}

# ---- load metadata + relative abundance (only the needed feature rows) ----
meta <- load_sex_metadata(hc_only = (scope == "hc"))
meta <- meta[meta[[COL_SEX]] %in% c("F", "M") &
             meta[[COL_SITE]] == SITE_SALIVA &
             meta[[COL_CONTINENT]] %in% group_order, , drop = FALSE]
meta$group_lab <- paste(meta[[COL_CONTINENT]], meta[[COL_SEX]])
group_lab_order <- as.vector(t(outer(continents, c("F", "M"), paste)))

# read only the selected feature rows from the 2.4 GB relabund tsv
# (grep -F with a trailing tab anchors each pattern to the feature column)
read_relabund_rows <- function(feats) {
  pats <- paste0(feats, "\t")
  cmd <- paste0("head -1 ", RELABUND_TSV, "; grep -F -e '",
                paste(pats, collapse = "' -e '"), "' ", RELABUND_TSV)
  fread(cmd = cmd, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
}

abundance_plots <- function(feats, min_m) {
  rel <- read_relabund_rows(feats)
  if (nrow(rel) < length(feats)) {
    warning("Only ", nrow(rel), "/", length(feats),
            " features found in the relabund table")
  }
  samples <- intersect(colnames(rel), rownames(meta))
  sdf <- as.data.frame(t(rel[, samples]))
  colnames(sdf) <- rel[[1]]
  sdf$group_lab <- meta[rownames(sdf), "group_lab"]

  # real relabund log2FC per feature x continent (ratio of group arithmetic
  # means, F = ref, M = tgt; pseudocount 1e-6)
  sdf_cont <- sub(" [FM]$", "", sdf$group_lab)
  sdf_sex  <- sub("^.* ", "", sdf$group_lab)
  real_lfc_df <- do.call(rbind, lapply(continents, function(ct) {
    idx <- which(sdf_cont == ct)
    if (length(idx) == 0) return(NULL)
    lfc <- real_lfc_matrix(sdf[idx, setdiff(colnames(sdf), "group_lab"), drop = FALSE],
                           sdf_sex[idx], ref = "F", tgt = "M")
    data.frame(feature = names(lfc), continent = ct, real_lfc = as.numeric(lfc))
  }))

  mean_by_group <- sdf %>%
    pivot_longer(-group_lab, names_to = "feature", values_to = "abundance") %>%
    group_by(group_lab, feature) %>%
    summarise(mean_abundance = mean(abundance), .groups = "drop")

  taxa_means <- data.frame(
    feature = rel[[1]],
    mean_abundance = rowMeans(rel[, samples])
  ) %>% arrange(desc(mean_abundance))

  # row order: global mean abundance, most abundant at top
  lab_map <- setNames(enh$genus_lab[match(taxa_means$feature, enh$feature)],
                      taxa_means$feature)
  row_levels <- rev(unname(lab_map))

  mean_by_group$genus_lab <- factor(lab_map[mean_by_group$feature], levels = row_levels)
  mean_by_group$group_lab <- factor(mean_by_group$group_lab, levels = group_lab_order)

  mean_heatmap <- ggplot(mean_by_group,
                         aes(x = group_lab, y = genus_lab,
                             fill = log10(mean_abundance))) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_viridis(option = "mako", direction = -1, name = "log10(Mean\nAbundance)") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, face = "italic", hjust = 1),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "bold"),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))

  bar_df <- taxa_means
  bar_df$genus_lab <- factor(lab_map[bar_df$feature], levels = row_levels)
  bar_breaks <- pretty(c(0, max(bar_df$mean_abundance, na.rm = TRUE)), n = 4)

  mean_barplot <- ggplot(bar_df) +
    geom_col(aes(x = mean_abundance, y = genus_lab), fill = "steelblue", width = 0.8) +
    scale_x_continuous(trans = "sqrt", breaks = bar_breaks,
                       expand = expansion(mult = c(0, 0.05))) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Mean Relative\nAbundance (sqrt)") +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
      plot.margin = margin(t = 0, r = 20, b = 0, l = 0)
    )

  # ---- differential-abundance bubble panel (left of the combined figure) ----
  # rows = same genera, columns = continents (F vs M). Point size = -log10(max q)
  # where max q is the LARGEST padj among the significant methods for that
  # feature x continent (the most conservative significance); color = real
  # relabund log2FC (group-mean ratio, >0 = enriched in MALES).
  sub <- enh[enh$feature %in% feats, , drop = FALSE]
  sub$real_lfc <- real_lfc_df$real_lfc[match(
    paste(sub$feature, sub$continent),
    paste(real_lfc_df$feature, real_lfc_df$continent))]
  diff_bubble <- NULL
  if (nrow(sub) > 0) {
    sub$continent <- factor(sub$continent, levels = continents)
    sub$genus_lab <- factor(enh$genus_lab[match(sub$feature, enh$feature)],
                            levels = row_levels)

    sig_mat  <- as.matrix(sub[, c("in_linda", "in_maaslin2", "in_ancombc2")])
    q_mat    <- as.matrix(sub[, c("linda_padj", "maaslin2_padj", "ancombc2_padj")])
    lfc_mat2 <- as.matrix(sub[, c("linda_lfc", "maaslin2_lfc", "ancombc2_lfc")])
    q_sig    <- q_mat;   q_sig[!sig_mat]   <- NA
    lfc_sig  <- lfc_mat2; lfc_sig[!sig_mat] <- NA
    sub$maxq <- apply(q_sig, 1, function(x)
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    sub$median_lfc <- apply(lfc_sig, 1, function(x)
      if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE))
    sub <- sub[!is.na(sub$maxq) & !is.na(sub$median_lfc) &
               rowSums(sig_mat) >= min_m, , drop = FALSE]

    if (nrow(sub) > 0) {
      sub$neg_log10_maxq <- -log10(pmax(sub$maxq, .Machine$double.xmin))
      diff_bubble <- ggplot(sub, aes(x = continent, y = genus_lab,
                                     size = neg_log10_maxq, color = real_lfc)) +
        geom_point(alpha = 0.9) +
        scale_size_continuous(range = c(1, 6), name = "-log10(max q)") +
        scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                              midpoint = 0, na.value = "grey80",
                              name = "log2FC\n(rel. abundance,\n>0 = males)") +
        scale_x_discrete(expand = c(0.02, 0.02)) +
        scale_y_discrete(expand = expansion(add = 0.6)) +
        labs(x = NULL, y = NULL) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 10, face = "italic", hjust = 1),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1, face = "bold"),
          axis.ticks = element_blank(),
          panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
          legend.position = "top",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9)
        )
    }
  }

  combined <- if (is.null(diff_bubble)) mean_heatmap | mean_barplot else
    diff_bubble | mean_heatmap + theme(axis.text.y = element_blank(),
                                       axis.ticks.y = element_blank()) | mean_barplot
  combined <- combined +
    plot_layout(widths = if (is.null(diff_bubble)) c(6, 7) else c(7, 8, 7))
  list(heatmap = mean_heatmap, barplot = mean_barplot, diff_bubble = diff_bubble,
       combined = combined)
}

# ==========================================================================
# 4. Top 5 Genera Abundance Distribution (Figure C, density ridgeline)
#    Among the consensus genera (same selection as Figure B), take the top 5
#    by global mean abundance and draw log10(relative abundance) ridgelines
#    per site x sex group. Zero abundances are floored at 1e-10 before log10.
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

  ggplot(genus_long) +
    geom_density_ridges(aes(x = log10(abundance), y = group_lab, fill = group_lab),
                        quantile_lines = TRUE, quantiles = 2, alpha = 0.7) +
    scale_fill_manual(values = rep(unname(sex_colors), length.out = length(groups_present))) +
    facet_grid(~ genus_lab) +
    labs(x = "log10(Relative Abundance)", y = NULL) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.ticks.y = element_blank(),
      strip.text = element_text(size = 12, face = "bold.italic"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
      plot.margin = margin(10, 10, 10, 10)
    )
}

for (ver in c("2method", "3method")) {
  min_m <- if (ver == "2method") 2L else 3L
  feats <- select_heatmap_features(enh, min_m)
  cat("    ", ver, ": ", length(feats), " abundance features\n", sep = "")
  if (length(feats) == 0) next
  ab <- abundance_plots(feats, min_m)
  save_both(file.path(figures_dir, paste0("figureB_abundance_heatmap_", ver, ".pdf")),
            ab$heatmap, 9, 12)
  save_both(file.path(figures_dir, paste0("figureB_abundance_barplot_", ver, ".pdf")),
            ab$barplot, 6, 12)
  if (!is.null(ab$diff_bubble)) {
    save_both(file.path(figures_dir, paste0("figureB_da_bubble_", ver, ".pdf")),
              ab$diff_bubble, 9, 12)
  }
  save_both(file.path(figures_dir, paste0("figureB_abundance_", ver, "_combined.pdf")),
            ab$combined, 17, 12)
  # Figure C: top 5 consensus genera by mean abundance — density ridgelines
  top5 <- top5_abundance_plot(feats)
  save_both(file.path(figures_dir, paste0("figureC_top5_abundance_", ver, ".pdf")),
            top5, 14, 7)
}

cat("==========================================\n")
cat("Done. Figures in:", figures_dir, "\n")
cat("  figureA_consensus_counts_per_continent.{pdf,png}\n")
cat("  figureB_da_bubble_{2method,3method}.{pdf,png}\n")
cat("  figureB_abundance_heatmap_{2method,3method}.{pdf,png}\n")
cat("  figureB_abundance_barplot_{2method,3method}.{pdf,png}\n")
cat("  figureB_abundance_{2method,3method}_combined.{pdf,png}\n")
cat("  figureC_top5_abundance_{2method,3method}.{pdf,png}\n")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")
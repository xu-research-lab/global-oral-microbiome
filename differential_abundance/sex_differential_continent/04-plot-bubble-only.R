# ==============================================================================
# 04-plot-bubble-only.R  (sex-within-continent, saliva, HC) - shared genera
# Single bubble plot: columns = continents, rows = genera.
#   point COLOR  = real relabund log2FC (log2 ratio of group-mean relative
#                  abundance, M vs F; red = enriched in males M,
#                  blue = enriched in females F)
#   point SIZE   = -log10(max q) where max q is the LARGEST (most conservative)
#                  padj among the significant methods for that feature x continent
#
# Only robust features (significant in >=2 methods) are shown.
# Restricts x-axis to genera robust in >=2 continents (shared signal).
#
# Input:  <out>/consensus_all_enhanced.tsv
# Output: <out>/figures/figure_bubble_per_continent.{pdf,png}
#
# Usage (from project root):
#   export DA_RESULTS_ROOT=differential_abundance/result
#   Rscript differential_abundance/sex_differential_continent/04-plot-bubble-only.R --scope=hc
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

source("differential_abundance/sex_differential_continent/da_helpers.R")
source("differential_abundance/real_lfc_helper.R")   # real_lfc_matrix(), REAL_LFC_PS

args <- commandArgs(trailingOnly = TRUE)
scope <- parse_scope_arg(args)

root <- Sys.getenv("DA_RESULTS_ROOT"); if (!nzchar(root)) root <- "differential_abundance/result"
out_dir <- file.path(continent_root(scope), "sex_continent_da_consensus")
figures_dir <- file.path(out_dir, "figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

scope_tag <- if (scope == "hc") "HC" else "all"

cat("==========================================\n")
cat("Bubble plot — sex within continent (saliva, shared genera), scope: ", scope_tag, "\n", sep = "")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ==========================================================================
# 1. Load enhanced consensus
# ==========================================================================
enh_file <- file.path(out_dir, "consensus_all_enhanced.tsv")
if (!file.exists(enh_file)) stop("Missing: ", enh_file)
enh <- fread(enh_file, data.table = FALSE)
cat("    Enhanced rows: ", nrow(enh), "\n", sep = "")

# Parse continent from comparison name ("saliva_<continent>_F_vs_M")
parse_continent <- function(comp)
  gsub("_", " ", sub("^saliva_", "", sub("_F_vs_M$", "", comp)))
enh$continent <- parse_continent(enh$comparison)

# Genus label: genus_short, fallback to last non-empty taxonomy level
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

uni_feat <- unique(enh$feature)
lab_map <- vapply(uni_feat, genus_label, character(1))
lab_map <- make.unique(lab_map)
names(lab_map) <- uni_feat
enh$genus_lab <- lab_map[enh$feature]

# ==========================================================================
# 2. Restrict to robust features (>=2 methods significant)
# ==========================================================================
rob <- enh[enh$robust == TRUE, , drop = FALSE]
cat("    Robust features (>=2 methods): ", nrow(rob), "\n", sep = "")
if (nrow(rob) == 0) stop("No robust features to plot.")

# Per feature x continent: max q over SIGNIFICANT methods only (point size)
sig_cols <- c("in_linda", "in_maaslin2", "in_ancombc2")
q_cols   <- c("linda_padj", "maaslin2_padj", "ancombc2_padj")

sig_mat  <- as.matrix(rob[, sig_cols])
q_mat    <- as.matrix(rob[, q_cols])

# Mask to significant methods only
q_sig   <- q_mat;   q_sig[!sig_mat]   <- NA

rob$maxq <- apply(q_sig, 1, function(x)
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))

# Drop rows without a usable max q
rob <- rob[!is.na(rob$maxq), , drop = FALSE]
cat("    Plottable robust points: ", nrow(rob), "\n", sep = "")
if (nrow(rob) == 0) stop("No plottable points after stat computation.")

# keep only features that are robust in >= 2 continents (shared signal, not a
# single-continent hit); keyed on full `feature` so same-genus distinct features
# are not merged
n_continent <- tapply(rob$continent, rob$feature, function(x) length(unique(x)))
keep_feat <- names(n_continent)[n_continent >= 2]
rob <- rob[rob$feature %in% keep_feat, , drop = FALSE]
cat("    Features robust in >=2 continents: ", length(keep_feat),
    " (", nrow(rob), " points)\n", sep = "")
if (length(keep_feat) == 0) stop("No feature is robust in >=2 continents.")

# -log10(max q), floored to avoid log(0)
rob$neg_log10_maxq <- -log10(pmax(rob$maxq, .Machine$double.xmin))

# ==========================================================================
# 2b. Real relabund log2FC per feature x continent
#     (log2 ratio of group-mean relative abundance, F = ref, M = tgt,
#     pseudocount 1e-6), computed on the saliva samples of this scope
# ==========================================================================
continent_y_order <- c("Europe", "Asia", "North America", "Oceania", "South America")

meta <- load_sex_metadata(hc_only = (scope == "hc"))
meta <- meta[meta[[COL_SEX]] %in% c("F", "M") &
             meta[[COL_SITE]] == SITE_SALIVA &
             meta[[COL_CONTINENT]] %in% continent_y_order, , drop = FALSE]

# read only the kept feature rows from the 2.4 GB relabund tsv
# (grep -F with a trailing tab anchors each pattern to the feature column)
pats <- paste0(keep_feat, "\t")
cmd <- paste0("head -1 ", RELABUND_TSV, "; grep -F -e '",
              paste(pats, collapse = "' -e '"), "' ", RELABUND_TSV)
rel <- fread(cmd = cmd, sep = "\t", header = TRUE, check.names = FALSE,
             data.table = FALSE)
if (nrow(rel) < length(keep_feat))
  warning("Only ", nrow(rel), "/", length(keep_feat),
          " features found in the relabund table")

real_lfc_df <- do.call(rbind, lapply(continent_y_order, function(ct) {
  m_ct <- meta[meta[[COL_CONTINENT]] == ct, , drop = FALSE]
  samples <- intersect(colnames(rel), rownames(m_ct))
  if (length(samples) == 0) return(NULL)
  sdf <- as.data.frame(t(rel[, samples, drop = FALSE]))
  colnames(sdf) <- rel[[1]]
  lfc <- real_lfc_matrix(sdf, m_ct[rownames(sdf), COL_SEX], ref = "F", tgt = "M")
  data.frame(feature = names(lfc), continent = ct, real_lfc = as.numeric(lfc))
}))
rob$real_lfc <- real_lfc_df$real_lfc[match(
  paste(rob$feature, rob$continent),
  paste(real_lfc_df$feature, real_lfc_df$continent))]
dropped_na_lfc <- sum(is.na(rob$real_lfc))
rob <- rob[!is.na(rob$real_lfc), , drop = FALSE]
if (dropped_na_lfc > 0)
  cat("    Dropped ", dropped_na_lfc, " points without real relabund LFC\n", sep = "")

# ==========================================================================
# 3. Ordering
# ==========================================================================
# Keep continent ordering
continents <- continent_y_order[continent_y_order %in% unique(rob$continent)]
rob$continent <- factor(rob$continent, levels = continents)

# Representative LFC per genus = the continent row with max |real_lfc|
genus_stat <- rob %>%
  group_by(genus_lab) %>%
  summarise(rep_lfc = real_lfc[which.max(abs(real_lfc))], .groups = "drop")

# X-axis ordering: M enriched (LFC > 0) on left, F enriched (LFC < 0) on right,
# within each block by descending |LFC|
genus_ordered <- genus_stat %>%
  mutate(
    sign = ifelse(rep_lfc > 0, "M", "F"),
    abs_lfc = abs(rep_lfc)
  ) %>%
  arrange(desc(sign), desc(abs_lfc))

rob$genus_lab <- factor(rob$genus_lab, levels = genus_ordered$genus_lab)

# ==========================================================================
# 4. Plot (x = genera, y = continents)
# ==========================================================================
p <- ggplot(rob, aes(x = genus_lab, y = continent,
                     size = neg_log10_maxq, fill = real_lfc)) +
  # shape 21 = circle with a dark outline, so near-zero (white) points
  # remain visible against the white panel background
  geom_point(shape = 21, colour = "grey25", stroke = 0.4, alpha = 0.9) +
  scale_size_continuous(range = c(2, 5), name = "-log10(max q)") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, na.value = "grey80",
                       name = "log2FC") +
  labs(x = NULL, y = NULL,
       title = NULL,
       subtitle = NULL) +
  theme_minimal(base_size = 12) +
  scale_y_discrete(position = "right") +
  guides(
    fill = guide_colorbar(barheight = unit(0.9, "cm"), barwidth = unit(0.18, "cm"),
                          title.position = "top", title.hjust = 0),
    size  = guide_legend(keyheight = unit(0.18, "cm"), keywidth = unit(0.18, "cm"),
                         title.position = "top", title.hjust = 0,
                         label.position = "right")
  ) +
  theme(
    axis.text.x = element_text(size = 9, face = "italic", hjust = 1, vjust = 1,
                               angle = 45),
    axis.text.y = element_text(size = 11, face = "bold", hjust = 0,
                               margin = margin(l = 6)),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    axis.ticks.length = unit(2, "pt"),
    panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
    legend.position = "right",
    legend.box.just = "left",
    legend.key.size = unit(0.22, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -2),
    legend.spacing.y = unit(0.05, "cm"),
    legend.title = element_text(size = 9, face = "bold", hjust = 0, margin = margin(b = 1)),
    legend.text = element_text(size = 8),
    # left margin: y axis is on the right, so without extra left padding the
    # first rotated x label (Anaeroglobus_38550) gets clipped at the canvas edge
    plot.margin = margin(2, 2, 2, 60, unit = "pt")
  )



out_pdf <- file.path(figures_dir, "figure_bubble_per_continent.pdf")
out_png <- file.path(figures_dir, "figure_bubble_per_continent.png")
ggsave(out_pdf, p, width = 7.4, height = 2.5, limitsize = FALSE)
ggsave(out_png, p, width = 7.4, height = 2.5, dpi = 300, bg = "white", limitsize = FALSE)

cat("\n==========================================\n")
cat("Done.\n")
cat("  ", out_pdf, "\n", sep = "")
cat("End: ", format(Sys.time()), "\n")
cat("==========================================\n")
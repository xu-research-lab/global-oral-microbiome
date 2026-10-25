# ==============================================================================
# 03-plot-lefse-bars.R  (sex_differential_continent/pooled — POOLED saliva F_vs_M)
# LEfSe-style vertical bar plot of the POOLED consensus robust genera
# (>=2 methods, padj<0.05), ALL continents merged (continent as covariate).
#
#   Bar height = effect, sign = enrichment direction (M => up, F => down):
#     Figure A: real relative-abundance log2FC (group-mean ratio, ps=1e-6,
#               F = reference)
#     Figure B: Maaslin2 effect size (coef, TSS+LOG model, log2 scale)
#
# Input:  <out>/consensus_robust.tsv          (direction_consensus)
#         <out>/full_table.tsv               (maaslin2_lfc)
#         data/metadata-and-biom-v2-qc/genus_table_relabund.tsv (real LFC)
# Output: <out>/figures/figure_lefse_real_lfc.{pdf,png}
#         <out>/figures/figure_lefse_maaslin2_coef.{pdf,png}
#
# Usage (from project root):
#   Rscript differential_abundance/sex_differential_continent/pooled/03-plot-lefse-bars.R --scope=hc
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(data.table)
})

source("differential_abundance/sex_differential_continent/da_helpers.R")   # rd(), COL_*, load_sex_metadata
source("differential_abundance/real_lfc_helper.R")                         # real_lfc_matrix(), REAL_LFC_PS

args <- commandArgs(trailingOnly = TRUE)
scope <- parse_scope_arg(args)

root    <- Sys.getenv("DA_RESULTS_ROOT"); if (!nzchar(root)) root <- "differential_abundance/result"
out_dir <- if (scope == "hc") rd("sex_pooled_hc/sex_pooled_da_consensus") else
                              rd("sex_pooled/sex_pooled_da_consensus")
figures_dir <- file.path(out_dir, "figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

SEX_REF <- "F"; SEX_TGT <- "M"
sex_colors <- c("F" = "#C51B7D", "M" = "#1F78B4")

save_both <- function(path_pdf, plot, width, height, dpi = 300) {
  ggsave(path_pdf, plot, width = width, height = height, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", path_pdf), plot, width = width, height = height,
         dpi = dpi, bg = "white", limitsize = FALSE)
}

cat("==========================================\n")
cat("LEfSe-style bars — pooled saliva F_vs_M (", scope, ")\n", sep = "")
cat("Out: ", out_dir, "/figures\n", sep = "")
cat("Start: ", format(Sys.time()), "\n")
cat("==========================================\n\n")

# ==========================================================================
# 1. Robust features + directions + Maaslin2 coef
# ==========================================================================
robust_file <- file.path(out_dir, "consensus_robust.tsv")
full_file   <- file.path(out_dir, "full_table.tsv")
for (f in c(robust_file, full_file)) if (!file.exists(f)) stop("Missing: ", f)

rob <- fread(robust_file, data.table = FALSE)
full <- fread(full_file, data.table = FALSE)
rob <- merge(rob, full[, c("feature", "maaslin2_lfc", "maaslin2_padj")],
             by = "feature", all.x = TRUE, sort = FALSE)
cat("    Robust features:", nrow(rob), "\n")

# bar sign = consensus direction (agreed:F negative, agreed:M positive);
# guard against a method-specific coef that flips sign
rob$dir_sign <- ifelse(grepl("agreed:M", rob$direction_consensus, fixed = TRUE), 1L,
                ifelse(grepl("agreed:F", rob$direction_consensus, fixed = TRUE), -1L, NA_integer_))
if (any(is.na(rob$dir_sign))) stop("Feature without agreed direction in robust table!")
# all-three taxa get a "*" suffix on the label
rob$all3 <- rob$all_three == "TRUE"   # kept for console reporting only

# genus label: genus_short with fallback to the last non-empty taxon level
label_one <- function(feature) {
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
rob$genus_lab <- vapply(rob$feature, label_one, character(1), USE.NAMES = FALSE)

# disambiguate duplicates across taxonomic lines
d <- as.data.frame(table(rob$genus_lab))
dup_lab <- d$Var1[d$Freq > 1]
for (l in dup_lab) {
  idx <- which(rob$genus_lab == l)
  parts <- strsplit(rob$feature[idx], ";", fixed = TRUE)
  fam <- vapply(parts, function(p) {
    p <- p[nzchar(p)]; i <- grep("^f__", p)
    if (length(i)) sub("^f__", "", p[i]) else ""}, character(1))
  rob$genus_lab[idx] <- paste0(l, " (", fam, ")")
}
if (any(rob$all3)) cat("    All-3-methods genera:", sum(rob$all3), "\n")
rob$genus_lab <- factor(rob$genus_lab)   # levels set per figure below

# ==========================================================================
# 2. Real relative-abundance log2FC over the pooled analysis set
#    (saliva, HC by scope, F vs M across all continents)
# ==========================================================================
meta <- load_sex_metadata(hc_only = (scope == "hc"))
meta <- meta[meta[[COL_SEX]] %in% c(SEX_REF, SEX_TGT) &
             meta[[COL_SITE]] == SITE_SALIVA, , drop = FALSE]
cat("    Pooled analysis samples: F=", sum(meta[[COL_SEX]] == SEX_REF),
    " / M=", sum(meta[[COL_SEX]] == SEX_TGT), "\n", sep = "")

pats <- paste0(rob$feature, "\t")
cmd <- paste0("head -1 ", RELABUND_TSV, "; grep -F -e '",
              paste(pats, collapse = "' -e '"), "' ", RELABUND_TSV)
rel <- fread(cmd = cmd, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
if (nrow(rel) < nrow(rob)) stop("Only ", nrow(rel), "/", nrow(rob),
                                " robust features found in the relabund table")
feats <- rel[[1]]; rn <- colnames(rel)
if (anyDuplicated(feats)) stop("duplicated feature rows in relabund table")
sdf <- as.data.frame(t(rel[, -1]))
colnames(sdf) <- as.character(feats)
sdf <- sdf[rownames(sdf) %in% rownames(meta), , drop = FALSE]
sdf[[COL_SEX]] <- meta[rownames(sdf), COL_SEX]
cat("    Relabund matched samples: ", nrow(sdf), "\n", sep = "")
rob$real_lfc <- real_lfc_matrix(sdf[, rob$feature, drop = FALSE],
                                sdf[[COL_SEX]], ref = SEX_REF, tgt = SEX_TGT)[rob$feature]

# ==========================================================================
# 3. One LEfSe-style vertical bar figure per metric
# ==========================================================================
draw_lefse <- function(df, value_col, y_lab, top_n = TOP_N) {
  d <- df[order(df[[value_col]]), , drop = FALSE]          # 从负到正排列（左侧为F富集，右侧为M富集）
  # keep only the top-N each side (positive = M, negative = F)
  pos <- utils::tail(d[d[[value_col]] > 0, , drop = FALSE], top_n)
  neg <- utils::head(d[d[[value_col]] < 0, , drop = FALSE], top_n)
  d <- rbind(neg, pos)                                     # 左负(F)右正(M)
  if (nrow(d) == 0) stop("no bars to draw for ", value_col)
  d$genus_lab <- factor(d$genus_lab, levels = d$genus_lab)
  
  # side/color follow the SIGN OF THE PLOTTED METRIC
  d$grp <- ifelse(d[[value_col]] > 0, SEX_TGT, SEX_REF)
  flipped <- d$genus_lab[sign(d[[value_col]]) != d$dir_sign]
  if (length(flipped))
    cat("    [note] sign contradicts model consensus: ", length(flipped),
        " -> ", paste(flipped, collapse = ", "), "\n", sep = "")
  
  ggplot(d, aes(x = genus_lab, y = .data[[value_col]], fill = grp)) +
    geom_col(width = 0.66, color = "grey20", linewidth = 0.2) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.45) +
    scale_fill_manual(values = sex_colors,
                      labels = c("F" = "Enriched in F", "M" = "Enriched in M"),
                      guide = guide_legend(title = NULL, reverse = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0.13, 0.13))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = y_lab) +
    theme_bw() +
    theme(axis.text.x = element_text(face = "italic", size = 9, 
                                     angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 10),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          legend.margin = margin(0, 0, 0, 0),
          legend.spacing = unit(0, "pt"),
          legend.background = element_blank(),
          legend.text = element_text(size = 9),
          # Reserve space for the first right-aligned, 45-degree genus label.
          plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 25, unit = "pt"),
          plot.title = element_blank(),
          plot.subtitle = element_blank())
}

TOP_N <- 20L
# 动态计算宽度以适应横向排列的标签，高度固定
fig_w <- 6
fig_h <- 3

p1 <- draw_lefse(rob, "real_lfc",
                 "log2 fold change (real mean relative abundance, F -> M)")
save_both(file.path(figures_dir, "figure_lefse_real_lfc.pdf"), p1, fig_w, fig_h)

p2 <- draw_lefse(rob, "maaslin2_lfc", "Maaslin2 coefficient")
save_both(file.path(figures_dir, "figure_lefse_maaslin2_coef.pdf"), p2, fig_w, fig_h)

cat("\nDONE: figures written to ", figures_dir, "\n", sep = "")

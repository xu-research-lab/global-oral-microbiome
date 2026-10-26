suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(patchwork)
  library(viridis)
  library(grid)
})
sym_limits <- function(x) {
  m <- max(abs(x), na.rm = TRUE)
  c(-m, m)
}

# Read selected rows without putting taxonomy labels into a shell command.
read_relabund_rows <- function(feats) {
  if (!length(feats)) stop("No features to read")
  patterns <- tempfile()
  on.exit(unlink(patterns), add = TRUE)
  writeLines(paste0(feats, "\t"), patterns)
  cmd <- paste(
    "head -n 1", shQuote(RELABUND_TSV), "; LC_ALL=C rg -F -f",
    shQuote(patterns), "--", shQuote(RELABUND_TSV)
  )
  rel <- fread(cmd = cmd, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  if (anyDuplicated(rel[[1]]) || !all(feats %in% rel[[1]])) {
    stop("Missing or duplicate feature rows in the relative-abundance table")
  }
  rel
}

genus_label <- function(feature) {
  g <- make_genus_short(feature)
  if (nzchar(g)) {
    return(g)
  }
  parts <- strsplit(feature, ";", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  lab <- ""
  for (p in rev(parts)) {
    lab <- sub("^[a-z]__", "", p)
    if (nzchar(lab)) break
  }
  if (!nzchar(lab)) "unclassified" else lab
}

age_levels <- if (analysis == "age") AGE_ORDER else character()
if (analysis == "age") {
  age_display <- setNames(c("0-6 months", "6-36 months", "3-12 years", "12-18 years", "18-60 years", ">60 year"), age_levels)
  age_short <- setNames(c("0-6m", "6-36m", "3-12y", "12-18y", "18-60y", ">60y"), age_levels)
  age_colors <- setNames(viridis(length(age_levels)), age_levels)
}
short_target <- function(g) paste0(if (analysis == "age") age_short[g] else g, " vs rest")

save_final <- function(plot, width, height, cairo = TRUE) {
  dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(plot_output_dir, figure_stem)
  ggsave(paste0(path, ".pdf"), plot,
    width = width, height = height,
    device = if (cairo) grDevices::cairo_pdf else "pdf", limitsize = FALSE
  )
  # Match the device used by the selected figures.
  ggsave(paste0(path, ".png"), plot,
    width = width, height = height, dpi = 300,
    device = ragg::agg_png, bg = "white", limitsize = FALSE
  )
  message("Saved ", path, ".{png,pdf}")
}

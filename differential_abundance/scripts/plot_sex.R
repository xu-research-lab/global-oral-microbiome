# Pooled HC sex: MaAsLin2 coefficient, up to 20 genera per direction.
robust_file <- file.path(out_dir, "consensus_robust.tsv")
full_file <- file.path(out_dir, "full_table.tsv")
for (f in c(robust_file, full_file)) if (!file.exists(f)) stop("Missing: ", f)

rob <- fread(robust_file, data.table = FALSE)
full <- fread(full_file, data.table = FALSE)
rob <- merge(rob, full[, c("feature", "maaslin2_lfc", "maaslin2_padj")],
  by = "feature", all.x = TRUE, sort = FALSE
)
cat("    Robust features:", nrow(rob), "\n")

# bar sign = consensus direction (agreed:F negative, agreed:M positive);
# guard against a method-specific coef that flips sign
rob$dir_sign <- ifelse(grepl("agreed:M", rob$direction_consensus, fixed = TRUE), 1L,
  ifelse(grepl("agreed:F", rob$direction_consensus, fixed = TRUE), -1L, NA_integer_)
)
if (any(is.na(rob$dir_sign))) stop("Feature without agreed direction in robust table!")
# all-three taxa get a "*" suffix on the label
rob$all3 <- rob$all_three == "TRUE" # kept for console reporting only

# genus label: genus_short with fallback to the last non-empty taxon level
label_one <- function(feature) {
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
rob$genus_lab <- vapply(rob$feature, label_one, character(1), USE.NAMES = FALSE)

# disambiguate duplicates across taxonomic lines
d <- as.data.frame(table(rob$genus_lab))
dup_lab <- d$Var1[d$Freq > 1]
for (l in dup_lab) {
  idx <- which(rob$genus_lab == l)
  parts <- strsplit(rob$feature[idx], ";", fixed = TRUE)
  fam <- vapply(parts, function(p) {
    p <- p[nzchar(p)]
    i <- grep("^f__", p)
    if (length(i)) sub("^f__", "", p[i]) else ""
  }, character(1))
  rob$genus_lab[idx] <- paste0(l, " (", fam, ")")
}
if (any(rob$all3)) cat("    All-3-methods genera:", sum(rob$all3), "\n")
rob$genus_lab <- factor(rob$genus_lab) # levels set per figure below


draw_lefse <- function(df, value_col, y_lab, top_n = TOP_N) {
  d <- df[order(df[[value_col]]), , drop = FALSE] # 从负到正排列（左侧为F富集，右侧为M富集）
  # keep only the top-N each side (positive = M, negative = F)
  pos <- utils::tail(d[d[[value_col]] > 0, , drop = FALSE], top_n)
  neg <- utils::head(d[d[[value_col]] < 0, , drop = FALSE], top_n)
  d <- rbind(neg, pos) # 左负(F)右正(M)
  if (nrow(d) == 0) stop("no bars to draw for ", value_col)
  d$genus_lab <- factor(d$genus_lab, levels = d$genus_lab)

  # side/color follow the SIGN OF THE PLOTTED METRIC
  d$grp <- ifelse(d[[value_col]] > 0, SEX_TGT, SEX_REF)
  flipped <- d$genus_lab[sign(d[[value_col]]) != d$dir_sign]
  if (length(flipped)) {
    cat("    [note] sign contradicts model consensus: ", length(flipped),
      " -> ", paste(flipped, collapse = ", "), "\n",
      sep = ""
    )
  }

  ggplot(d, aes(x = genus_lab, y = .data[[value_col]], fill = grp)) +
    geom_col(width = 0.66, color = "grey20", linewidth = 0.2) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.45) +
    scale_fill_manual(
      values = sex_colors,
      labels = c("F" = "Enriched in F", "M" = "Enriched in M"),
      guide = guide_legend(title = NULL, reverse = TRUE)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.13, 0.13))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = y_lab) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        face = "italic", size = 9,
        angle = 45, hjust = 1, vjust = 1
      ),
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
      plot.subtitle = element_blank()
    )
}


TOP_N <- 20L
sex_colors <- c(F = "#C51B7D", M = "#1F78B4")
plot <- draw_lefse(rob, "maaslin2_lfc", "Maaslin2 coefficient")
save_final(plot, 6, 3, cairo = FALSE)

#!/usr/bin/env Rscript
# pipeline_v2: ggtree rectangular phylotype tree
# DEREP ASV tips colored by phylotype novelty class; GG2 refs in gray.
# Heatmap rings: subjects, BioProjects, continents.
suppressMessages({
  library(ggtree)
  library(ape)
  library(tidyverse)
  library(ggnewscale)
})

PIPE_V2 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
REPORT  <- file.path(PIPE_V2, "analysis_report")
FIGOUT  <- file.path(REPORT, "figs")
dir.create(FIGOUT, showWarnings = FALSE)

# ---- Colors (colorblind-safe) ----
COL_UNCAP  <- "#E15759"   # homd_uncaptured
COL_CAP    <- "#4E79A7"   # homd_captured
COL_GG2    <- "#D9D9D9"   # GG2 reference tips
INK        <- "#2B2B2B"

# ---- Load tree ----
tree <- read.tree(file.path(REPORT, "phylotype_subtree_top150.nwk"))
cat(sprintf("Tree tips: %d\n", length(tree$tip.label)))

# ---- Tip types ----
tip_labels <- tree$tip.label
is_derep <- grepl("^DEREP_", tip_labels)
cat(sprintf("  DEREP tips: %d | GG2 ref tips: %d\n",
            sum(is_derep), sum(!is_derep)))

# ---- Load DEREP → phylotype metadata ----
d2p <- read.delim(file.path(REPORT, "derep_to_phylotype.tsv"),
                  sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  d2p rows: %d\n", nrow(d2p)))

# ---- Load phylotype table for continent info ----
phy_all <- read.delim(file.path(PIPE_V2, "06_phylotypes_novelty_classified.tsv"),
                      sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
# Map phylotype -> continents (Python set literal -> parse)
parse_continents <- function(x) {
  if (is.na(x) || x == "") return("")
  y <- gsub("[{}' ]", "", x)
  if (y == "") return("")
  # keep as comma list
  gsub(";", ",", y)
}
phy_cont <- setNames(sapply(phy_all$continents, parse_continents), phy_all$phylotype_id)

# ---- Build tip annotation table ----
# For DEREP tips: join metadata; for GG2 refs: mark as reference
d2p$tip <- as.character(d2p$tip)
tip_annot <- data.frame(
  tip   = tree$tip.label,
  type  = ifelse(is_derep, "derep", "gg2"),
  stringsAsFactors = FALSE
)
tip_annot <- merge(tip_annot, d2p, by = "tip", all.x = TRUE)
rownames(tip_annot) <- tip_annot$tip

# novelty label (factor for legend) — binary (2026-08-13)
tip_annot$novelty_label <- "GG2 reference"
tip_annot$novelty_label[tip_annot$type == "derep" & tip_annot$novelty == "homd_captured"] <- "HOMD-captured"
tip_annot$novelty_label[tip_annot$type == "derep" & tip_annot$novelty == "homd_uncaptured"] <- "HOMD-uncaptured"
tip_annot$novelty_label <- factor(tip_annot$novelty_label,
                                  levels = c("HOMD-uncaptured", "HOMD-captured", "GG2 reference"))

# support values (defaults for GG2 refs)
tip_annot$subjects <- as.numeric(tip_annot$subjects)
tip_annot$bps     <- as.numeric(tip_annot$bps)
tip_annot$subjects[is.na(tip_annot$subjects)] <- 0
tip_annot$bps[is.na(tip_annot$bps)] <- 0

# continents per tip (from phylotype)
tip_annot$cont_list <- sapply(tip_annot$phylotype, function(p) {
  if (is.na(p) || !p %in% names(phy_cont)) return("")
  phy_cont[p]
})
n_cont <- function(cont_str) {
  if (is.na(cont_str) || cont_str == "") return(0)
  length(unique(trimws(unlist(strsplit(cont_str, ",")))))
}
tip_annot$n_continents <- sapply(tip_annot$cont_list, n_cont)

cat(sprintf("  DEREP tips annotated: %d\n", sum(tip_annot$type == "derep" & !is.na(tip_annot$phylotype))))
cat(sprintf("  Novelty distribution (DEREP):\n"))
print(table(tip_annot$novelty[tip_annot$type == "derep"]))

# ---- Order tree tips to match annotation ----
tree <- keep.tip(tree, tree$tip.label)  # no-op safeguard
p0 <- ggtree(tree, layout = "rectangular", size = 0.15) %<+% tip_annot

# ---- Base plot: tip color by novelty ----
p <- p0 +
  geom_tippoint(aes(color = novelty_label), size = 0.6) +
  scale_color_manual(name = "Novelty",
                     values = c("HOMD-uncaptured" = COL_UNCAP,
                                "HOMD-captured" = COL_CAP,
                                "GG2 reference" = COL_GG2),
                     breaks = c("HOMD-uncaptured", "HOMD-captured", "GG2 reference"))

# ---- Heatmap rings (gheatmap: rownames = tip labels, one value column) ----
# Log1p transform for reads-like distributions; subjects & continents as-is
ring_subj <- data.frame(Subjects = tip_annot$subjects, row.names = tip_annot$tip)
ring_bp   <- data.frame(BioProjects = tip_annot$bps, row.names = tip_annot$tip)
ring_cont <- data.frame(Continents = tip_annot$n_continents, row.names = tip_annot$tip)

p <- p + new_scale_fill()
p <- gheatmap(p, ring_subj, offset = 0.01, width = 0.06,
              colnames = TRUE, colnames_angle = 90, colnames_offset_y = 0.5,
              font.size = 2.5) +
  scale_fill_viridis_c(option = "C", name = "Subjects", na.value = "#F5F5F5")

p <- p + new_scale_fill()
p <- gheatmap(p, ring_bp, offset = 0.09, width = 0.06,
              colnames = TRUE, colnames_angle = 90, colnames_offset_y = 0.5,
              font.size = 2.5) +
  scale_fill_viridis_c(option = "D", name = "BioProjects", na.value = "#F5F5F5")

p <- p + new_scale_fill()
p <- gheatmap(p, ring_cont, offset = 0.17, width = 0.06,
              colnames = TRUE, colnames_angle = 90, colnames_offset_y = 0.5,
              font.size = 2.5) +
  scale_fill_viridis_c(option = "B", name = "Continents", na.value = "#F5F5F5")

# Save
ggsave(file.path(FIGOUT, "fig_phylotype_tree_rect.png"),
       p + theme_tree2() +
         ggtitle("Top-150 candidate novel oral phylotypes on GG2 tree") +
         theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
               legend.position = "right"),
       width = 11, height = 22, dpi = 300)

cat("Saved: fig_phylotype_tree_rect.png\n")
cat("Done!\n")

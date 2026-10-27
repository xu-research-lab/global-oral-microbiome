#!/usr/bin/env Rscript
# pipeline_v2: phylotype-circle tree (CIRCULAR) for the strict candidates.
# Tree pruned to ONLY the representative ASV tips (no GG2 reference tips).
# Each phylotype = ONE circle:
#   size  = number of member ASVs
#   fill  = novelty depth (red = novel lineage outside HOMD genera, blue = novel
#           species within a known oral genus) — from the oral_clade annotation
#   label = species name (GG2 taxonomy), fallback genus / phylotype id
# Rings (ggtreeExtra::geom_fruit): Prevalence, BioProjects, Continents, Regions.
suppressMessages({
  library(ggtree)
  library(ape)
  library(tidyverse)
  library(ggnewscale)
  library(ggtreeExtra)
})

PIPE_V2 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
REPORT  <- file.path(PIPE_V2, "analysis_report")
FIGOUT  <- file.path(REPORT, "figs")
dir.create(FIGOUT, showWarnings = FALSE)

INK <- "#2B2B2B"
COL_OUTSIDE <- "#E15759"   # outside HOMD genus (oral_clade = FALSE)
COL_IN_HOMD  <- "#4E79A7"   # within HOMD genus (oral_clade = TRUE)

# ---- Load subtree (strict-candidate tips incl. GG2 refs) ----
# suffix auto-derived from the strict-candidate count (09_build_strict_tree.py)
SUFFIX <- "strict71"
tree <- read.tree(file.path(REPORT, paste0("phylotype_subtree_", SUFFIX, ".nwk")))
cat(sprintf("Subtree tips: %d\n", length(tree$tip.label)))

# ---- Load anchors ----
anchors <- read.delim(file.path(PIPE_V2, "07_phylotype_anchors.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  anchors: %d\n", nrow(anchors)))

# ---- Prune to ONLY the rep tips (drop all GG2 refs) ----
keep <- as.character(anchors$anchor_tip)
missing <- setdiff(keep, tree$tip.label)
cat(sprintf("  anchors missing from tree: %d\n", length(missing)))
drop <- setdiff(tree$tip.label, keep)
tree <- ape::drop.tip(tree, drop)   # collapse.singles = TRUE (default)
cat(sprintf("  rep-only tree tips: %d\n", length(tree$tip.label)))

# ---- Tip annotation (all tips are phylotype reps) ----
tip_annot <- data.frame(tip = tree$tip.label, stringsAsFactors = FALSE)
tip_annot <- merge(tip_annot, anchors, by.x = "tip", by.y = "anchor_tip", all.x = TRUE)
rownames(tip_annot) <- tip_annot$tip
for (col in c("max_subject_prevalence", "n_bioprojects_total", "n_continents",
              "n_regions", "n_members", "n_flagged")) {
  tip_annot[[col]] <- as.numeric(tip_annot[[col]])
  tip_annot[[col]][is.na(tip_annot[[col]])] <- 0
}
# oral_clade arrives as "True"/"False" (Python bool str) → normalize for the
# manual fill scale; NA = no placement-clade genus
tip_annot$oral_clade <- toupper(as.character(tip_annot$oral_clade))

# ---- Circular base plot ----
p0 <- ggtree(tree, layout = "circular", size = 0.3) %<+% tip_annot
rep_data <- p0$data[p0$data$isTip, ]

p <- p0 +
  geom_tippoint(data = rep_data,
                aes(size = log1p(n_members), fill = oral_clade),
                shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(name = "Novelty depth",
                    values = c("TRUE" = COL_IN_HOMD, "FALSE" = COL_OUTSIDE),
                    breaks = c("TRUE", "FALSE"),
                    labels = c("novel species in known oral genus",
                               "novel lineage outside HOMD genera"),
                    na.value = "grey70") +
  scale_size_continuous(name = "member ASVs (log)", range = c(0.2, 3)) +
  theme(legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8))

# ---- Tree labels (novel-to-HOMD genera; set by 06d_apply_tree_labels.py) ----
# label only the defensible novel-to-oral-catalog phylotypes (top-30); the
# circle fill encodes the oral clade (red = outside HOMD genus, blue = in HOMD)
# read.delim reads empty cells as "" not NA — filter both
labelled_tips <- tip_annot$tip[!is.na(tip_annot$tree_label) & tip_annot$tree_label != ""]
label_data <- rep_data[rep_data$label %in% labelled_tips, ]
cat(sprintf("  labelled tips (novel genera): %d / %d\n",
            nrow(label_data), nrow(rep_data)))
p <- p + geom_tiplab2(data = label_data, aes(label = tree_label),
                      size = 2.0, offset = 0.04, color = INK)

# ---- Rings (ggtreeExtra; matches fruit data to tips by the `tip` column) ----
# note: `%<+%` drops the `tip` column from p0$data, so build from tip_annot
# subject-level prevalence (distinct (bioproject, subject_id) individuals)
ring_prev <- data.frame(tip = tip_annot$tip, value = tip_annot$max_subject_prevalence)
ring_bp   <- data.frame(tip = tip_annot$tip, value = tip_annot$n_bioprojects_total)
ring_cont <- data.frame(tip = tip_annot$tip, value = tip_annot$n_continents)

p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_prev, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.7,
                    offset = 0.05, pwidth = 0.06,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "C", name = "Prevalence", na.value = "grey92")

p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_bp, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.7,
                    offset = 0.05, pwidth = 0.06,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "D", name = "BioProjects", na.value = "grey92")

p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_cont, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.7,
                    offset = 0.05, pwidth = 0.06,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "B", name = "Continents", na.value = "grey92")

p <- p + ggtitle("71 candidate novel oral phylotypes — circle size = member ASVs") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

ggsave(file.path(FIGOUT, "fig_tree_phylotype_circles.png"), p,
       width = 14, height = 14, dpi = 300, limitsize = FALSE)
ggsave(file.path(FIGOUT, "fig_tree_phylotype_circles.pdf"), p,
       width = 14, height = 14, limitsize = FALSE)
cat("Saved: fig_tree_phylotype_circles.png / .pdf\n")
cat("Done!\n")

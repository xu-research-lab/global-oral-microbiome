#!/usr/bin/env Rscript
# pipeline_v2: CIRCULAR phylotype-circle tree for the SUPPORTED CATALOG (186).
# Adapted from 10_ggtree_phylotype_circles.R. Tree pruned to ONLY the rep tips.
# Each phylotype = ONE circle:
#   size = log member ASVs; fill = oral_clade (red outside / blue within HOMD genus)
# Rings: Prevalence, BioProjects, Continents, Regions.
suppressMessages({
  library(ggtree)
  library(ape)
  library(tidyverse)
  library(ggnewscale)
  library(ggtreeExtra)
})

PIPE_V2 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
REPORT  <- file.path(PIPE_V2, "analysis_report")
FIGOUT  <- file.path(REPORT, "figs_v2")
dir.create(FIGOUT, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  SUFFIX <- args[1]
} else {
  # auto-detect the latest phylotype_subtree_supported*.nwk (current = supported258)
  cand <- list.files(REPORT, pattern = "^phylotype_subtree_supported[0-9]+\\.nwk$")
  if (length(cand) >= 1) {
    nums <- as.integer(sub("^phylotype_subtree_supported([0-9]+)\\.nwk$", "\\1", cand))
    SUFFIX <- paste0("supported", max(nums))
  } else {
    SUFFIX <- "supported258"
  }
}
cat(sprintf("SUFFIX: %s\n", SUFFIX))

INK <- "#2B2B2B"
COL_OUTSIDE <- "#E15759"
COL_IN_HOMD  <- "#4E79A7"

tree <- read.tree(file.path(REPORT, paste0("phylotype_subtree_", SUFFIX, ".nwk")))
cat(sprintf("Subtree tips: %d\n", length(tree$tip.label)))

anchors <- read.delim(file.path(PIPE_V2, "12d_anchors_supported.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  anchors: %d\n", nrow(anchors)))

# phylum per phylotype from step 12g profile — use the TREE-BASED phylum
# (placement-edge descendant-tip majority, consistent with the tree);
# fallback to the vsearch phylum when tree phylum is missing.
prof <- read.delim(file.path(PIPE_V2, "12g_phylotype_profile.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)
prof$phylum <- ifelse(prof$tree_phylum != "", prof$tree_phylum, prof$vsearch_phylum)
phy_count <- sort(table(prof$phylum[prof$phylum != ""]), decreasing = TRUE)
prof$phylum_group <- ifelse(prof$phylum == "", "unknown", prof$phylum)
anchors <- merge(anchors, prof[, c("phylotype_id", "phylum_group")],
                 by = "phylotype_id", all.x = TRUE)
# annotation-status flag: is EVERY member ASV's GG2 identity < 97%?
anchors <- merge(anchors, prof[, c("phylotype_id", "phylo_max_gg2_id")],
                 by = "phylotype_id", all.x = TRUE)
anchors$all_members_lt97 <- is.na(anchors$phylo_max_gg2_id) |
                            anchors$phylo_max_gg2_id < 97
anchors$annot_status <- ifelse(anchors$all_members_lt97,
                               "all members <97%", "has >=97% member")

keep <- as.character(anchors$anchor_tip)
missing <- setdiff(keep, tree$tip.label)
cat(sprintf("  anchors missing from tree: %d\n", length(missing)))
drop <- setdiff(tree$tip.label, keep)
tree <- ape::drop.tip(tree, drop)
cat(sprintf("  rep-only tree tips: %d\n", length(tree$tip.label)))

tip_annot <- data.frame(tip = tree$tip.label, stringsAsFactors = FALSE)
tip_annot <- merge(tip_annot, anchors, by.x = "tip", by.y = "anchor_tip", all.x = TRUE)
rownames(tip_annot) <- tip_annot$tip
for (col in c("max_subject_prevalence", "n_bioprojects_total", "n_continents",
              "n_regions", "n_members")) {
  tip_annot[[col]] <- as.numeric(tip_annot[[col]])
  tip_annot[[col]][is.na(tip_annot[[col]])] <- 0
}
tip_annot$oral_clade <- toupper(as.character(tip_annot$oral_clade))
tip_annot$tree_label[is.na(tip_annot$tree_label) | tip_annot$tree_label == ""] <- NA
tip_annot$phylum_group[is.na(tip_annot$phylum_group) | tip_annot$phylum_group == ""] <- "unknown"
tip_annot$annot_status[is.na(tip_annot$annot_status)] <- "all members <97%"

p0 <- ggtree(tree, layout = "circular", size = 0.3) %<+% tip_annot
rep_data <- p0$data[p0$data$isTip, ]

p <- p0 +
  geom_tippoint(data = rep_data,
                aes(size = log2(n_members + 1), fill = oral_clade),
                shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(name = "Oral clade",
                    values = c("TRUE" = COL_IN_HOMD, "FALSE" = COL_OUTSIDE),
                    breaks = c("TRUE", "FALSE"),
                    labels = c("within HOMD genus", "outside HOMD genera"),
                    na.value = "grey70") +
  scale_size_continuous(name = "member ASVs (log2)", range = c(0.3, 3)) +
  theme(legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8))

# tip labels intentionally omitted (186 circles too dense to read)

ring_prev <- data.frame(tip = tip_annot$tip, value = tip_annot$max_subject_prevalence)
ring_bp   <- data.frame(tip = tip_annot$tip, value = tip_annot$n_bioprojects_total)
ring_cont <- data.frame(tip = tip_annot$tip, value = tip_annot$n_continents)
ring_stat <- data.frame(tip = tip_annot$tip, status = tip_annot$annot_status)
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_prev, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "C", name = "Prevalence", na.value = "grey92",
                       guide = guide_colorbar(order = 1))
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_bp, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "D", name = "BioProjects", na.value = "grey92",
                       guide = guide_colorbar(order = 2))
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_cont, geom = geom_tile,
                    mapping = aes(y = tip, fill = value), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_viridis_c(option = "B", name = "Continents", na.value = "grey92",
                       guide = guide_colorbar(order = 3))
# ---- Annotation status ring (all members <97% GG2?) ----
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_stat, geom = geom_tile,
                    mapping = aes(y = tip, fill = status), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_manual(values = c("all members <97%" = "#F2F2F2",
                               "has >=97% member" = "#5B6E8C"),
                    name = "GG2 anchor",
                    guide = guide_legend(order = 4, keyheight = 0.5, keywidth = 0.5,
                                         title.theme = element_text(size = 8),
                                         label.theme = element_text(size = 6)))
# ---- Phylum ring (categorical) ----
ring_phy <- data.frame(tip = tip_annot$tip, phylum = tip_annot$phylum_group)
# legend order = count descending (most frequent phylum first)
phy_levels <- names(sort(table(ring_phy$phylum), decreasing = TRUE))
ring_phy$phylum <- factor(ring_phy$phylum, levels = phy_levels)
phy_colors <- hcl.colors(max(3, length(phy_levels)), "Spectral")
names(phy_colors) <- phy_levels
phy_colors["unknown"] <- "white"
p <- p + new_scale_fill()
p <- p + geom_fruit(data = ring_phy, geom = geom_tile,
                    mapping = aes(y = tip, fill = phylum), width = 0.35,
                    offset = 0.05, pwidth = 0.05,
                    axis.params = list(axis = "none")) +
  scale_fill_manual(values = phy_colors,
                    breaks = setdiff(phy_levels, "unknown"),
                    name = "Phylum",
                    guide = guide_legend(ncol = 3, keyheight = 0.5, keywidth = 0.5,
                                         order = 5,
                                         title.theme = element_text(size = 8),
                                         label.theme = element_text(size = 6)))

p <- p + ggtitle(sprintf("Supported catalog (%d clean phylotypes) — circle size = member ASVs",
                         length(tree$tip.label))) +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

ggsave(file.path(FIGOUT, paste0("fig_tree_circles_", SUFFIX, ".png")), p,
       width = 14, height = 14, dpi = 300, limitsize = FALSE)
ggsave(file.path(FIGOUT, paste0("fig_tree_circles_", SUFFIX, ".pdf")), p,
       width = 14, height = 14, limitsize = FALSE)
cat("Saved:", paste0("fig_tree_circles_", SUFFIX, ".png / .pdf"), "\nDone!\n")

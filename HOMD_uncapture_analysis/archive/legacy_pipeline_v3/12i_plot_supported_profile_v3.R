#!/usr/bin/env Rscript
# pipeline_v3: SUPPORTED CATALOG profile (replaces fig_tree_circles_supported_v3).
# Non-tree layout, one row per supported cluster (486), sorted by:
#   1) phylum (count desc, "unknown" last)  2) max subject prevalence (desc)
# Panels: Phylum | No. derep ASVs | BioProjects | Continents | Regions |
#         nearest_ref_dist | HOMD identity | annotation in HOMD (species) |
#         HOMD capture | Region matrix (no x labels) | per-ASV prevalence box+jitter
suppressPackageStartupMessages({
  library(ggplot2)
  library(aplot)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

split_values <- function(x) {
  if (is.na(x) || x == "" || x == "NA" || x == "nan") return(character())
  vals <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
  vals[vals != "" & vals != "NA" & vals != "nan"]
}
clean <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA" | x == "nan" | x == "Unassigned"] <- NA_character_
  x
}

message("[12i_v3] loading supported-catalog anchors + master ASV annotation")
anchors <- read.delim(file.path(PIPE_V3, "12d_anchors_supported_v3.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
master  <- read.delim(file.path(PIPE_V3, "07_master_asv_table.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
prof    <- read.delim(file.path(PIPE_V3, "12d_headline_profile_v3.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# ---- per-ASV annotation columns (rep ASV = anchor_tip) ----
mcols <- c("derep_id", "nearest_ref_dist", "homd_identity", "annotation_in_homd",
           "annotation", "annotation_level", "annotation_lineage",
           "amplicon_regions", "continents", "subject_prevalence")
master <- master[, intersect(mcols, colnames(master))]
master$derep_id <- as.character(master$derep_id)
anchors$anchor_tip <- as.character(anchors$anchor_tip)

rep <- merge(anchors, master, by.x = "anchor_tip", by.y = "derep_id", all.x = TRUE)
stopifnot(nrow(rep) == nrow(anchors))

# cluster-level homd_capture from headline profile (one row per derep ASV)
prof$derep_id <- as.character(prof$derep_id)
cap <- unique(prof[, c("derep_id", "homd_capture", "cluster_regions", "cluster_continents")])
rep <- merge(rep, cap, by.x = "anchor_tip", by.y = "derep_id", all.x = TRUE)

# ---- numeric coercion ----
for (cc in c("n_members", "max_subject_prevalence", "n_bioprojects_total",
             "n_continents", "n_regions", "nearest_ref_dist", "homd_identity")) {
  rep[[cc]] <- suppressWarnings(as.numeric(rep[[cc]]))
}
rep$annotation <- clean(rep$annotation)
rep$annotation_level <- clean(rep$annotation_level)
# phylum taken from annotation_lineage (2nd rank), NOT gg2_phylum
rep$phylum <- vapply(strsplit(rep$annotation_lineage, ";", fixed = TRUE), function(x) {
  if (length(x) >= 2) trimws(x[2]) else ""
}, character(1))
rep$phylum[is.na(rep$phylum) | rep$phylum == ""] <- "unknown"
rep$homd_capture[is.na(rep$homd_capture)] <- "homd_uncaptured"
rep$annotation_in_homd[is.na(rep$annotation_in_homd)] <- FALSE
rep$label <- ifelse(is.na(rep$annotation),
                    paste0(rep$cluster_id, " | ", rep$tree_label),
                    paste0(rep$cluster_id, " | ", rep$annotation))

# ---- sorting: phylum -> nearest_ref_dist (asc) -> max subject prevalence (desc) ----
phy_tab <- sort(table(rep$phylum), decreasing = TRUE)
phy_levels <- names(phy_tab)
if ("unknown" %in% phy_levels) {
  phy_levels <- c(setdiff(phy_levels, "unknown"), "unknown")
}
rep$phylum <- factor(rep$phylum, levels = phy_levels)
rep <- rep[order(rep$phylum, rep$nearest_ref_dist, -rep$max_subject_prevalence,
                 rep$cluster_id), ]
rep$rank <- seq_len(nrow(rep))
rep$rank_f <- factor(rep$rank, levels = rev(rep$rank))

n <- nrow(rep)
message(sprintf("[12i_v3] %d supported clusters plotted (phylum groups: %d)",
                n, length(phy_levels)))

# ---- adaptive sizing (compressed to normal height) ----
base_sz <- 6
ylab_sz <- 2.2
xlab_sz <- 6
row_in  <- 0.010
fig_h   <- max(5.0, 1.2 + n * row_in)   # ~6.1 in for 486 rows
thin    <- ceiling(n / 60)

theme_h <- theme_minimal(base_size = base_sz) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = xlab_sz),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(2, 4, 2, 2))

mk_tile <- function(d, fill_var, scale_fun) {
  p <- ggplot(d, aes(y = rank_f, x = track, fill = !!sym(fill_var))) +
    geom_tile(width = 1, height = 1) +
    scale_fun +
    theme_h +
    theme(axis.text.y = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.title = element_text(size = 5.5),
          legend.text = element_text(size = 4.8),
          legend.key.size = unit(0.16, "cm"),
          legend.margin = margin(0, 0, 0, 0))
  p
}

# ---- panel data ----
phy_df <- data.frame(rank_f = rep$rank_f, track = "Phylum", phylum = rep$phylum)
asv_df <- data.frame(rank_f = rep$rank_f, track = "No. derep ASVs (log10)",
                      value = log10(pmax(rep$n_members, 1)))
bp_df  <- data.frame(rank_f = rep$rank_f, track = "No. bioprojects (log10)",
                      value = log10(pmax(rep$n_bioprojects_total, 1)))
cont_df <- data.frame(rank_f = rep$rank_f, track = "No. continents", value = rep$n_continents)
dist_df <- data.frame(rank_f = rep$rank_f, track = "nearest_ref_dist",
                      dist_cat = ifelse(!is.na(rep$nearest_ref_dist) & rep$nearest_ref_dist <= 0.03,
                                        "<= 0.03", "> 0.03"))
dist_df$dist_cat <- factor(dist_df$dist_cat, levels = c("<= 0.03", "> 0.03"))
homd_df <- data.frame(rank_f = rep$rank_f, track = "HOMD identity", value = rep$homd_identity)
inhomd_df <- data.frame(rank_f = rep$rank_f, track = "Annotation in HOMD",
                        in_homd = ifelse(rep$annotation_in_homd, "yes", "no"))
# sites matrix (cluster-level union of member ASV sampling sites; replaces continents matrix)
prof$sites <- as.character(prof$sites)
sites_by_cluster <- aggregate(sites ~ cluster_id, data = prof, FUN = function(x) {
  vals <- unique(unlist(lapply(x, split_values)))
  paste(vals, collapse = ";")
})
sites_map <- setNames(sites_by_cluster$sites, as.character(sites_by_cluster$cluster_id))
rep$sites <- sites_map[as.character(rep$cluster_id)]
all_sites <- sort(unique(unlist(lapply(rep$sites, split_values))))
sites_mat <- do.call(rbind, lapply(seq_len(nrow(rep)), function(i) {
  vals <- split_values(rep$sites[i])
  data.frame(rank_f = rep$rank_f[i], track = all_sites, detected = all_sites %in% vals)
}))
sites_mat$track <- factor(sites_mat$track, levels = all_sites)

# cluster-level MAX subject prevalence (bar chart)
rank_map <- setNames(rep$rank_f, rep$cluster_id)
prof$subject_prevalence <- suppressWarnings(as.numeric(prof$subject_prevalence))
prev_max <- aggregate(subject_prevalence ~ cluster_id, data = prof,
                      FUN = function(x) max(x, na.rm = TRUE))
prev_max$rank_f <- rank_map[match(prev_max$cluster_id, names(rank_map))]
prev_max <- prev_max[!is.na(prev_max$rank_f), ]
p_prev <- ggplot(prev_max, aes(y = rank_f, x = subject_prevalence)) +
  geom_col(fill = "#2F5D62", width = 0.8) +
  labs(y = NULL, x = "Max subject prevalence") +
  theme_classic(base_size = base_sz) +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size = xlab_sz),
        plot.margin = margin(2, 4, 2, 2))

# ---- assemble ----
phy_colors <- hcl.colors(max(3, length(phy_levels)), "Spectral")
names(phy_colors) <- phy_levels
if ("unknown" %in% phy_levels) {
  phy_colors["unknown"] <- "#D9D9D9"
}

p_phy <- mk_tile(phy_df, "phylum",
                 scale_fill_manual(values = phy_colors, name = "Phylum", drop = FALSE)) +
  labs(title = sprintf("Supported catalog (%d clusters) — sorted by phylum, then nearest_ref_dist, then subject prevalence",
                       n),
       subtitle = "one row per supported cluster | heatmaps: No. derep ASVs (log10), bioprojects (log10), continents, nearest_ref_dist (<=0.03 vs >0.03), HOMD identity, annotation-in-HOMD (species), sites matrix | right: max subject prevalence (bar)") +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5))
p_asv <- mk_tile(asv_df, "value",
                 scale_fill_gradient(low = "#F1F3F1", high = "#5B6E8C", name = "No. derep ASVs (log10)"))
p_bp <- mk_tile(bp_df, "value",
                scale_fill_gradient(low = "#F1F3F1", high = "#244C5A", name = "No. bioprojects (log10)"))
p_cont <- mk_tile(cont_df, "value",
                  scale_fill_gradient(low = "#F1F3F1", high = "#7A5C3E", name = "No. continents"))
p_dist <- mk_tile(dist_df, "dist_cat",
                  scale_fill_manual(values = c("<= 0.03" = "#F7F7F7", "> 0.03" = "#B2182B"),
                                    name = "nearest_ref_dist", drop = FALSE))
p_homd <- mk_tile(homd_df, "value",
                  scale_fill_gradient(low = "#F6F1E8", high = "#1A9850",
                                      name = "HOMD identity (%)"))
p_inhomd <- mk_tile(inhomd_df, "in_homd",
                    scale_fill_manual(values = c("yes" = "#8ECAE6", "no" = "#F3F3F3"),
                                      name = "Annotation in HOMD (species)"))
p_sites_mat <- mk_tile(sites_mat, "detected",
                        scale_fill_manual(values = c(`FALSE` = "#F3F3F3", `TRUE` = "#4A7C59"),
                                          name = "Sites"))

combined <- p_phy %>%
  insert_right(p_asv, width = 0.12) %>%
  insert_right(p_bp, width = 0.12) %>%
  insert_right(p_cont, width = 0.12) %>%
  insert_right(p_dist, width = 0.12) %>%
  insert_right(p_homd, width = 0.12) %>%
  insert_right(p_inhomd, width = 0.12) %>%
  insert_right(p_sites_mat, width = 0.84) %>%
  insert_right(p_prev, width = 0.5)
combined$width <- c(0.20, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.84, 0.5)

out_base <- file.path(FIGOUT, "fig_v3_supported_profile")
ggsave(paste0(out_base, ".pdf"), combined, width = 8.0, height = fig_h,
       device = cairo_pdf, limitsize = FALSE)
ggsave(paste0(out_base, ".png"), combined, width = 8.0, height = fig_h, dpi = 150,
       limitsize = FALSE)
cat(sprintf("Saved: %s.png / .pdf (height %.1f in, %d clusters)\nDone!\n",
            out_base, fig_h, n))

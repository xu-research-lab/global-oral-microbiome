#!/usr/bin/env Rscript
# pipeline_v2: HEADLINE candidates (33) prevalence profile — HORIZONTAL layout.
# Adapted from ASV_vs_HOMD_analysis/region_aware_pipeline/17_plot_bacterial_asv_prevalence_profile.R
# y-axis = phylotype rank by subject prevalence (species annotation labels).
# Panels (left -> right): subject prevalence (horizontal bars), No. ASVs,
# BioProjects, log10(reads), HOMD identity bins, Continents, Sites.
suppressPackageStartupMessages({
  library(ggplot2)
  library(aplot)
})

PIPE_V2 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v2"
FIGOUT  <- file.path(PIPE_V2, "analysis_report", "figs_v2")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

split_values <- function(x) {
  if (is.na(x) || x == "" || x == "NA" || x == "nan") return(character())
  vals <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
  vals[vals != "" & vals != "NA" & vals != "nan"]
}

clean <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "Unassigned"] <- NA_character_
  x
}

message("[12f] loading headline profile")
df <- read.delim(file.path(PIPE_V2, "12d_headline_profile.tsv"),
                 sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
message("[12f] retained ", nrow(df), " headline phylotypes")

df$species <- clean(df$gg2_species)
df$genus <- clean(df$gg2_genus)
df$axis_label <- ifelse(!is.na(df$species), df$species,
                 ifelse(!is.na(df$genus), paste0(df$genus, " sp."),
                        gsub("EDGE_", "E", df$phylotype_id)))

df <- df[order(-df$subject_prevalence, -df$n_bioprojects, -df$total_reads), ]
df$rank <- seq_len(nrow(df))
df$rank_f <- factor(df$rank, levels = rev(df$rank))   # highest prevalence on top
df$log10_total_reads <- log10(df$total_reads + 1)
df$homd_identity <- suppressWarnings(as.numeric(df$homd_identity))

identity_bins <- c("<90%", "90-95%", "95-97%", "97-99%", ">=99%")
df$homd_identity_bin <- factor(
  cut(df$homd_identity, breaks = c(-Inf, 90, 95, 97, 99, Inf),
      labels = identity_bins, right = FALSE), levels = identity_bins)

asv_df <- data.frame(rank_f = df$rank_f, track = "No. ASVs", value = df$n_asvs)
bioproject_df <- data.frame(rank_f = df$rank_f, track = "No. bioprojects", value = df$n_bioprojects)
reads_df <- data.frame(rank_f = df$rank_f, track = "log10 reads", value = df$log10_total_reads)
homd_df <- data.frame(rank_f = df$rank_f, track = "HOMD id", identity_bin = df$homd_identity_bin)

all_continents <- sort(unique(unlist(lapply(df$continents, split_values))))
continent_df <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
  vals <- split_values(df$continents[i])
  data.frame(rank_f = df$rank_f[i], track = all_continents, detected = all_continents %in% vals)
}))
continent_df$track <- factor(continent_df$track, levels = all_continents)

all_sites <- sort(unique(unlist(lapply(df$sites, split_values))))
site_df <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
  vals <- split_values(df$sites[i])
  data.frame(rank_f = df$rank_f[i], track = all_sites, detected = all_sites %in% vals)
}))
site_df$track <- factor(site_df$track, levels = all_sites)

theme_h <- theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 6, color = "#222222"),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(2, 6, 2, 2))

# ---- far-right panel: subject prevalence (horizontal bars) ----
p_prev <- ggplot(df, aes(y = rank_f, x = subject_prevalence)) +
  geom_col(fill = "#2F5D62", width = 0.78) +
  labs(y = NULL, x = "Subject prevalence") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(size = 7),
        plot.margin = margin(2, 6, 2, 2))

mk_tile <- function(d, fill_var, scale_fun, show_y = FALSE) {
  p <- ggplot(d, aes(y = rank_f, x = track, fill = !!sym(fill_var))) +
    geom_tile(width = 0.88, height = 0.78) +
    scale_fun +
    theme_h +
    theme(axis.text.y = if (show_y) element_text(size = 4.5, color = "#222222")
                         else element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 5.5),
          legend.key.size = unit(0.18, "cm"),
          legend.margin = margin(0, 0, 0, 0))
  p
}

p_asv <- mk_tile(asv_df, "value",
                 scale_fill_gradient(low = "#F1F3F1", high = "#5B6E8C", name = "No. ASVs"),
                 show_y = TRUE) +
  scale_y_discrete(labels = setNames(df$axis_label, df$rank)) +
  labs(title = "HEADLINE candidates (33): GG2>=99% & HOMD-uncaptured",
       subtitle = "ordered by subject-level prevalence") +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7))
p_bp <- mk_tile(bioproject_df, "value",
                scale_fill_gradient(low = "#F1F3F1", high = "#244C5A", name = "No. bioprojects"))
p_reads <- mk_tile(reads_df, "value",
                   scale_fill_gradient(low = "#F6F1E8", high = "#9A4F2D", name = "log10 reads"))
identity_palette <- c("<90%" = "#D73027", "90-95%" = "#FC8D59", "95-97%" = "#FEE08B",
                      "97-99%" = "#91BFDB", ">=99%" = "#1A9850")
p_homd <- mk_tile(homd_df, "identity_bin",
                  scale_fill_manual(values = identity_palette, drop = TRUE, name = "HOMD id"))
p_continent <- mk_tile(continent_df, "detected",
                       scale_fill_manual(values = c(`FALSE` = "#F3F3F3", `TRUE` = "#C96B3C"),
                                         name = "Continents"))
p_site <- mk_tile(site_df, "detected",
                  scale_fill_manual(values = c(`FALSE` = "#F3F3F3", `TRUE` = "#4A7C59"),
                                    name = "Site"))

# aplot chain; override the width slot so ALL panels have explicit equal widths
# (the base panel would otherwise always be width 1)
combined <- p_asv %>%
  insert_right(p_bp, width = 0.15) %>%
  insert_right(p_reads, width = 0.15) %>%
  insert_right(p_homd, width = 0.15) %>%
  insert_right(p_continent, width = 0.65) %>%
  insert_right(p_site, width = 0.90) %>%
  insert_right(p_prev, width = 0.75)
combined$width <- c(0.15, 0.15, 0.15, 0.15, 0.65, 0.90, 0.75)

out_base <- file.path(FIGOUT, "fig_v2_headline_prevalence_profile")
ggsave(paste0(out_base, ".pdf"), combined, width = 7.0, height = 4.5, device = cairo_pdf)
ggsave(paste0(out_base, ".png"), combined, width = 7.0, height = 4.5, dpi = 300)
cat("Saved:", paste0(out_base, ".png / .pdf"), "\nDone!\n")

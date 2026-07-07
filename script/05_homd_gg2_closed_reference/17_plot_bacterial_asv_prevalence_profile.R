#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(aplot)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[[1]] else "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
infile <- file.path(outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_only.tsv")
figdir <- file.path(outdir, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

split_values <- function(x) {
  if (is.na(x) || x == "" || x == "NA" || x == "nan") return(character())
  vals <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
  vals <- vals[vals != "" & vals != "NA" & vals != "nan"]
  unique(vals)
}

make_label <- function(species, genus, phylum, derep_id) {
  clean <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == "" | x == "Unassigned"] <- NA_character_
    x
  }
  species <- clean(species)
  genus <- clean(genus)
  phylum <- clean(phylum)
  taxon <- ifelse(!is.na(species), species, ifelse(!is.na(genus), paste0(genus, " sp."), phylum))
  taxon[is.na(taxon)] <- "Unassigned"
  paste0(derep_id, "\n", taxon)
}

message("[fig10] loading bacterial ASV table")
df <- read.delim(infile, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("subject_prevalence" %in% names(df))

message("[fig10] filtering to GG2 identity=100% and query_cov=100%")
df <- df[df$gg2_identity == 100 & df$gg2_qcov == 100, ]
message("[fig10] retained ", nrow(df), " ASVs with perfect GG2 match")

df <- df[order(-df$subject_prevalence, -df$n_bioprojects, -df$total_reads), ]
df$rank <- seq_len(nrow(df))
df$rank_f <- factor(df$rank, levels = df$rank)
df$log10_total_reads <- log10(df$total_reads + 1)
df$homd_best_identity <- suppressWarnings(as.numeric(df$identity))
df$gg2_best_identity_for_plot <- suppressWarnings(as.numeric(df$gg2_identity))
identity_bins <- c("<90%", "90-95%", "95-97%", "97-99%", ">=99%")
bin_identity <- function(x) {
  out <- cut(
    x,
    breaks = c(-Inf, 90, 95, 97, 99, Inf),
    labels = identity_bins,
    right = FALSE
  )
  factor(out, levels = identity_bins)
}
df$homd_identity_bin <- bin_identity(df$homd_best_identity)
df$gg2_identity_bin <- bin_identity(df$gg2_best_identity_for_plot)
df$top10_label <- make_label(df$gg2_species, df$gg2_genus, df$gg2_phylum, df$derep_id)

top10 <- df[seq_len(min(10, nrow(df))), ]
x_limits <- c(0.5, nrow(df) + 0.5)

message("[fig10] preparing heatmap data")
bioproject_df <- data.frame(
  rank = df$rank,
  rank_f = df$rank_f,
  track = "No. bioprojects",
  value = df$n_bioprojects,
  stringsAsFactors = FALSE
)
reads_df <- data.frame(
  rank = df$rank,
  rank_f = df$rank_f,
  track = "log10(total reads + 1)",
  value = df$log10_total_reads,
  stringsAsFactors = FALSE
)
homd_identity_df <- data.frame(
  rank = df$rank,
  rank_f = df$rank_f,
  track = "Best HOMD identity",
  identity_bin = df$homd_identity_bin,
  stringsAsFactors = FALSE
)
gg2_identity_df <- data.frame(
  rank = df$rank,
  rank_f = df$rank_f,
  track = "Best GG2 identity",
  identity_bin = df$gg2_identity_bin,
  stringsAsFactors = FALSE
)

all_regions <- sort(unique(unlist(lapply(df$subject_geo_regions, split_values))))
region_df <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
  vals <- split_values(df$subject_geo_regions[i])
  data.frame(rank = df$rank[i], rank_f = df$rank_f[i], track = all_regions, detected = all_regions %in% vals, stringsAsFactors = FALSE)
}))
region_df$track <- factor(region_df$track, levels = rev(all_regions))

all_sites <- sort(unique(unlist(lapply(df$sites, split_values))))
site_df <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
  vals <- split_values(df$sites[i])
  data.frame(rank = df$rank[i], rank_f = df$rank_f[i], track = all_sites, detected = all_sites %in% vals, stringsAsFactors = FALSE)
}))
site_df$track <- factor(site_df$track, levels = rev(all_sites))

base_theme <- theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(2, 6, 2, 6)
  )

p_top <- ggplot(df, aes(x = rank, y = subject_prevalence)) +
  geom_line(color = "#2F5D62", linewidth = 0.45) +
  geom_point(color = "#2F5D62", size = 0.9, alpha = 0.9) +
  geom_point(data = top10, color = "#B33A3A", size = 1.7) +
  labs(
    x = NULL,
    y = "Subject prevalence",
    title = "Prevalent bacterial ASVs with perfect GG2 match (100% identity & coverage)",
    subtitle = paste0("n = ", nrow(df), " ASVs; ordered by subject-level prevalence; labels mark the top 10")
  ) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  base_theme +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    axis.text.y = element_text(color = "#222222")
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_top <- p_top + ggrepel::geom_text_repel(
    data = top10,
    aes(label = top10_label),
    size = 2.4,
    min.segment.length = 0,
    segment.color = "#777777",
    box.padding = 0.25,
    point.padding = 0.18,
    max.overlaps = Inf,
    seed = 17
  )
} else {
  p_top <- p_top + geom_text(
    data = top10,
    aes(label = top10_label),
    size = 2.2,
    angle = 35,
    hjust = 0,
    vjust = -0.25,
    check_overlap = TRUE
  )
}

p_bioproject <- ggplot(bioproject_df, aes(x = rank, y = track, fill = value)) +
  geom_tile(width = 1, height = 0.88) +
  scale_fill_gradient(low = "#F1F3F1", high = "#244C5A", name = "No. bioprojects") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "#222222"),
    legend.position = "right",
    plot.margin = margin(0, 6, 0, 6)
  )

p_reads <- ggplot(reads_df, aes(x = rank, y = track, fill = value)) +
  geom_tile(width = 1, height = 0.88) +
  scale_fill_gradient(low = "#F6F1E8", high = "#9A4F2D", name = "log10(total reads + 1)") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "#222222"),
    legend.position = "right",
    plot.margin = margin(0, 6, 0, 6)
  )

identity_palette <- c(
  "<90%" = "#D73027",
  "90-95%" = "#FC8D59",
  "95-97%" = "#FEE08B",
  "97-99%" = "#91BFDB",
  ">=99%" = "#1A9850"
)

p_homd_identity <- ggplot(homd_identity_df, aes(x = rank, y = track, fill = identity_bin)) +
  geom_tile(width = 1, height = 0.88) +
  scale_fill_manual(values = identity_palette, drop = FALSE, name = "Best HOMD identity") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "#222222"),
    legend.position = "right",
    plot.margin = margin(0, 6, 0, 6)
  )

p_gg2_identity <- ggplot(gg2_identity_df, aes(x = rank, y = track, fill = identity_bin)) +
  geom_tile(width = 1, height = 0.88) +
  scale_fill_manual(values = identity_palette, drop = FALSE, name = "Best GG2 identity") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "#222222"),
    legend.position = "right",
    plot.margin = margin(0, 6, 0, 6)
  )

p_region <- ggplot(region_df, aes(x = rank, y = track, fill = detected)) +
  geom_tile(width = 1, height = 0.88, color = NA) +
  scale_fill_manual(values = c(`FALSE` = "#F3F3F3", `TRUE` = "#C96B3C"), name = "Region") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "#222222"),
    legend.position = "right",
    plot.margin = margin(0, 6, 0, 6)
  )

p_site <- ggplot(site_df, aes(x = rank, y = track, fill = detected)) +
  geom_tile(width = 1, height = 0.88, color = NA) +
  scale_fill_manual(values = c(`FALSE` = "#F3F3F3", `TRUE` = "#4A7C59"), name = "Site") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(x = "ASV rank by subject prevalence", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "#222222"),
    legend.position = "right",
    plot.margin = margin(0, 6, 2, 6)
  )

message("[fig10] combining panels with aplot")
combined <- p_top %>%
  insert_bottom(p_bioproject, height = 0.07) %>%
  insert_bottom(p_reads, height = 0.07) %>%
  insert_bottom(p_homd_identity, height = 0.07) %>%
  insert_bottom(p_gg2_identity, height = 0.07) %>%
  insert_bottom(p_region, height = 0.26) %>%
  insert_bottom(p_site, height = 0.32)

out_base <- file.path(figdir, "fig10_bacterial_perfect_gg2_match_asv_prevalence_profile")
ggsave(paste0(out_base, ".pdf"), combined, width = 13.8, height = 8.8, device = cairo_pdf)
ggsave(paste0(out_base, ".png"), combined, width = 13.8, height = 8.8, dpi = 300)

summary_path <- file.path(outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_only.perfect_gg2_match_top10.tsv")
write.table(top10[, c("derep_id", "subject_prevalence", "n_bioprojects", "total_reads", "subject_geo_regions", "sites", "gg2_identity", "gg2_phylum", "gg2_genus", "gg2_species")],
            file = summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("[fig10] wrote: ", paste0(out_base, ".pdf"))
message("[fig10] wrote: ", paste0(out_base, ".png"))
message("[fig10] wrote: ", summary_path)

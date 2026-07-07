#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(aplot)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[[1]] else "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/containment_derep/gg2_backbone_suspicious_filter"
infile <- file.path(outdir, "funnel_layer5_asvs.tsv")
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
  taxon <- ifelse(!is.na(species), species,
                  ifelse(!is.na(genus), paste0(genus, " sp."), phylum))
  taxon[is.na(taxon)] <- "Unassigned"
  paste0(derep_id, "\n", taxon)
}

message("[funnel L5] loading ASV table")
df <- read.delim(infile, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("subject_prevalence" %in% names(df))

df <- df[order(-df$subject_prevalence, -df$n_bioprojects, -df$total_reads), ]
df$rank <- seq_len(nrow(df))
df$rank_f <- factor(df$rank, levels = df$rank)
df$log10_total_reads <- log10(df$total_reads + 1)
df$top10_label <- make_label(df$gg2_species, df$gg2_genus, df$gg2_phylum, df$derep_id)

message("[funnel L5] n = ", nrow(df), " ASVs")

top10 <- df[seq_len(min(10, nrow(df))), ]
x_limits <- c(0.5, nrow(df) + 0.5)

message("[funnel L5] preparing heatmap data")
bioproject_df <- data.frame(
  rank = df$rank, rank_f = df$rank_f,
  track = "No. bioprojects", value = df$n_bioprojects,
  stringsAsFactors = FALSE
)
reads_df <- data.frame(
  rank = df$rank, rank_f = df$rank_f,
  track = "log10(total reads + 1)", value = df$log10_total_reads,
  stringsAsFactors = FALSE
)

all_regions <- sort(unique(unlist(lapply(df$subject_geo_regions, split_values))))
if (length(all_regions) > 0) {
  region_df <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    vals <- split_values(df$subject_geo_regions[i])
    if (length(vals) == 0) vals <- "NA"
    data.frame(rank = df$rank[i], rank_f = df$rank_f[i],
               track = all_regions, detected = all_regions %in% vals,
               stringsAsFactors = FALSE)
  }))
  region_df$track <- factor(region_df$track, levels = rev(all_regions))
} else {
  region_df <- data.frame(rank = integer(), rank_f = factor(),
                          track = character(), detected = logical())
}

all_sites <- sort(unique(unlist(lapply(df$sites, split_values))))
if (length(all_sites) > 0) {
  site_df <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    vals <- split_values(df$sites[i])
    if (length(vals) == 0) vals <- "NA"
    data.frame(rank = df$rank[i], rank_f = df$rank_f[i],
               track = all_sites, detected = all_sites %in% vals,
               stringsAsFactors = FALSE)
  }))
  site_df$track <- factor(site_df$track, levels = rev(all_sites))
} else {
  site_df <- data.frame(rank = integer(), rank_f = factor(),
                        track = character(), detected = logical())
}

base_theme <- theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(2, 6, 2, 6)
  )

p_top <- ggplot(df, aes(x = rank, y = subject_prevalence)) +
  geom_line(color = "#2F5D62", linewidth = 0.5) +
  geom_point(color = "#2F5D62", size = 1.1, alpha = 0.9) +
  geom_point(data = top10, color = "#B33A3A", size = 2.0) +
  labs(
    x = NULL, y = "Subject prevalence",
    title = paste0("Funnel Layer 5: High-Confidence Bacterial ASVs (n = ", nrow(df), ")"),
    subtitle = ">=3 BioProjects, >=50 subjects, >=100 reads, GG2 100% match, median RA >=0.1%"
  ) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  base_theme +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9),
    axis.text.y = element_text(color = "#222222")
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_top <- p_top + ggrepel::geom_text_repel(
    data = top10, aes(label = top10_label),
    size = 2.4, min.segment.length = 0, segment.color = "#777777",
    box.padding = 0.25, point.padding = 0.18, max.overlaps = Inf, seed = 17
  )
} else {
  p_top <- p_top + geom_text(
    data = top10, aes(label = top10_label),
    size = 2.2, angle = 35, hjust = 0, vjust = -0.25, check_overlap = TRUE
  )
}

p_bioproject <- ggplot(bioproject_df, aes(x = rank, y = track, fill = value)) +
  geom_tile(width = 1, height = 0.88) +
  scale_fill_gradient(low = "#F1F3F1", high = "#244C5A", name = "No. bioprojects") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.y = element_text(color = "#222222"),
        legend.position = "right", plot.margin = margin(0, 6, 0, 6))

p_reads <- ggplot(reads_df, aes(x = rank, y = track, fill = value)) +
  geom_tile(width = 1, height = 0.88) +
  scale_fill_gradient(low = "#F6F1E8", high = "#9A4F2D", name = "log10(total reads + 1)") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.y = element_text(color = "#222222"),
        legend.position = "right", plot.margin = margin(0, 6, 0, 6))

panels <- list(p_top)

if (nrow(bioproject_df) > 0) {
  panels <- c(panels, list(p_bioproject))
  panels <- c(panels, list(p_reads))
}

if (nrow(region_df) > 0) {
  p_region <- ggplot(region_df, aes(x = rank, y = track, fill = detected)) +
    geom_tile(width = 1, height = 0.88, color = NA) +
    scale_fill_manual(values = c(`FALSE` = "#F3F3F3", `TRUE` = "#C96B3C"), name = "Region") +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.y = element_text(color = "#222222"),
          legend.position = "right", plot.margin = margin(0, 6, 0, 6))
  panels <- c(panels, list(p_region))
}

if (nrow(site_df) > 0) {
  p_site <- ggplot(site_df, aes(x = rank, y = track, fill = detected)) +
    geom_tile(width = 1, height = 0.88, color = NA) +
    scale_fill_manual(values = c(`FALSE` = "#F3F3F3", `TRUE` = "#4A7C59"), name = "Site") +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    labs(x = "ASV rank by subject prevalence", y = NULL) +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.y = element_text(color = "#222222"),
          legend.position = "right", plot.margin = margin(0, 6, 2, 6))
  panels <- c(panels, list(p_site))
}

message("[funnel L5] combining ", length(panels), " panels with aplot")
# Build panel heights
n_panels <- length(panels)
heights <- c(0.5, rep(0.07, min(n_panels - 1, 2)),
             if (n_panels > 3) 0.15 else NULL,
             if (n_panels > 4) 0.25 else NULL)
heights <- heights[!sapply(heights, is.null)]
heights <- heights / sum(heights)

combined <- panels[[1]]
if (n_panels >= 2) combined <- combined %>% insert_bottom(panels[[2]], height = 0.07)
if (n_panels >= 3) combined <- combined %>% insert_bottom(panels[[3]], height = 0.07)
if (n_panels >= 4) combined <- combined %>% insert_bottom(panels[[4]], height = 0.18)
if (n_panels >= 5) combined <- combined %>% insert_bottom(panels[[5]], height = 0.28)

out_base <- file.path(figdir, "funnel_layer5_prevalence_profile")
ggsave(paste0(out_base, ".pdf"), combined, width = 13.8, height = 7.5, device = cairo_pdf)
ggsave(paste0(out_base, ".png"), combined, width = 13.8, height = 7.5, dpi = 300)

summary_path <- file.path(outdir, "funnel_layer5_top10.tsv")
write.table(top10[, c("derep_id", "subject_prevalence", "n_bioprojects", "total_reads",
                       "sites", "subject_geo_regions", "gg2_phylum", "gg2_genus", "gg2_species")],
            file = summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("[funnel L5] wrote: ", paste0(out_base, ".pdf"))
message("[funnel L5] wrote: ", paste0(out_base, ".png"))
message("[funnel L5] wrote: ", summary_path)

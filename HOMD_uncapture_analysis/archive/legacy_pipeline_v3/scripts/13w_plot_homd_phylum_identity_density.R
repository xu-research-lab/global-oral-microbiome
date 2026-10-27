#!/usr/bin/env Rscript
# pipeline_v3: alternative density version of the four identity panels.
#   x = best HOMD identity (%), y = density (per identity unit).
#   Panels (2x2 facet, All phyla combined first):
#     All phyla combined | ASV count
#     Reads | Subject prevalence
# Run 13v_homd_phylum_identity_density.py first.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(ggnewscale)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

cov <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE)
dens_file <- file.path(PIPE_V3, "13v_homd_phylum_identity_density.tsv")
stopifnot(file.exists(dens_file))
dens <- read.delim(dens_file, sep = "\t", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)

# legend order follows frac_ge97 descending (as in the other version)
meta <- cov %>% transmute(phylum = phylum,
                          frac_ge97 = n_asvs_ge97 / n_asvs_total) %>%
  arrange(desc(frac_ge97))
meta$label <- sprintf("%s (%.1f%%)", meta$phylum, 100 * meta$frac_ge97)
meta$label <- factor(meta$label, levels = meta$label)
PAL <- setNames(hue_pal()(15), levels(meta$label))

# ---- per-phylum panels ----
phy <- dens %>% filter(group != "All phyla combined") %>%
  left_join(meta, by = c("group" = "phylum")) %>%
  transmute(panel = metric, identity, density, label)

# ---- all-phyla-combined panel ----
tot <- dens %>% filter(group == "All phyla combined")
tot$metric <- factor(tot$metric,
                     levels = c("ASV count", "Reads",
                                "Subject prevalence"))
tot$panel <- "All phyla combined"

phy$panel <- factor(phy$panel,
                    levels = c("ASV count", "Reads",
                               "Subject prevalence"))
panel_levels <- c("All phyla combined", "ASV count", "Reads",
                  "Subject prevalence")

base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.box = "vertical",
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 9),
        legend.key.size = unit(0.34, "cm"),
        legend.text = element_text(size = 6.3))

p <- ggplot() +
  geom_vline(xintercept = 97, linetype = "dashed", colour = "grey40",
             linewidth = 0.35) +
  geom_line(data = phy, aes(x = identity, y = density, colour = label),
            linewidth = 0.7) +
  scale_colour_manual(values = PAL, name = NULL, drop = FALSE,
                      guide = guide_legend(order = 2)) +
  new_scale_colour() +
  geom_line(data = tot, aes(x = identity, y = density, colour = metric),
            linewidth = 0.9) +
  scale_colour_manual(values = c("ASV count" = "grey40",
                                 "Reads" = "#1f77b4",
                                 "Subject prevalence" = "#d62728"),
                      name = "Metric",
                      guide = guide_legend(order = 1,
                                           direction = "vertical")) +
  facet_wrap(~ factor(panel, levels = panel_levels), ncol = 2) +
  scale_x_continuous(breaks = seq(70, 100, 5), limits = c(70, 100),
                     name = "Best identity to HOMD reference (%)") +
  scale_y_continuous(limits = c(0, NA),
                     name = "Density") +
  base_theme

fig <- file.path(FIGOUT, "fig_v3_homd_phylum_identity_density")
ggsave(paste0(fig, ".png"), p, width = 8, height = 6, dpi = 300)
ggsave(paste0(fig, ".pdf"), p, width = 8, height = 6)
message("[13w] density alternative written")

#!/usr/bin/env Rscript
# pipeline_v3: per-phylum IDENTITY DISTRIBUTION curves.
# x = cumulative fraction of phylum ASVs (from low to high identity),
# y = best identity to HOMD reference. One monotonic curve per phylum.
# Curves sitting high/right = better HOMD coverage.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

d   <- read.delim(file.path(PIPE_V3, "13f_homd_phylum_ecdf.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cov <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# --- rotate: x = cumulative fraction of ASVs with identity <= x (i.e. 1 - frac_ge) ---
d <- d %>% filter(identity >= 70) %>%
  mutate(cumfrac = 1 - frac_ge)   # fraction of phylum ASVs below this identity

# legend order by frac_ge97 (high -> low)
meta <- cov %>% transmute(phylum = phylum, frac_ge97 = n_asvs_ge97 / n_asvs_total) %>%
  arrange(desc(frac_ge97))
meta$label <- sprintf("%s (%.1f%%)", meta$phylum, 100 * meta$frac_ge97)
meta$label <- factor(meta$label, levels = meta$label)
d <- d %>% left_join(meta, by = "phylum")

pal <- hue_pal()(15)

p <- ggplot(d, aes(x = 100 * cumfrac, y = identity, colour = label)) +
  geom_hline(yintercept = 97, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  annotate("text", x = 99, y = 97.35, label = "97%", size = 3, colour = "grey25", hjust = 1) +
  scale_x_continuous(name = "Cumulative fraction of phylum ASVs (%)",
                     breaks = seq(0, 100, 20), limits = c(0, 100)) +
  scale_y_continuous(name = "Best identity to HOMD reference (%)",
                     breaks = seq(70, 100, 3), limits = c(70, 100.5)) +
  scale_colour_manual(values = pal, name = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 7.5), panel.grid.minor = element_blank())

fig <- file.path(FIGOUT, "fig_v3_homd_phylum_identity_distribution")
ggsave(paste0(fig, ".png"), p, width = 8.5, height = 5.5, dpi = 300)
ggsave(paste0(fig, ".pdf"), p, width = 8.5, height = 5.5)
message("[13m] identity distribution curves written")

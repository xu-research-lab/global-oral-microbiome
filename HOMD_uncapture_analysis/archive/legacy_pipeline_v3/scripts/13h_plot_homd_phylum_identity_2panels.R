#!/usr/bin/env Rscript
# pipeline_v3: per-phylum cumulative identity curves, two windows:
#   Fig 1: identity 70-97%, denominator = ALL phylum ASVs
#   Fig 2: identity 97-100%, denominator = only ASVs with identity >= 97%
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

d <- read.delim(file.path(PIPE_V3, "13f_homd_phylum_ecdf.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cov <- read.delim(file.path(PIPE_V3, "13d_homd_phylum_coverage.tsv"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# per-phylum metadata: total ASVs, >=97% ASVs, label
meta <- data.frame(
  phylum     = cov$phylum,
  n_total    = cov$n_asvs_total,
  n_ge97     = cov$n_asvs_ge97,
  stringsAsFactors = FALSE
)
meta$label <- with(meta, sprintf("%s  (n=%s; >=97%%: %s)",
                                 phylum,
                                 format(round(n_total/1000, 0), trim = TRUE),
                                 format(round(n_ge97/1000, 1), trim = TRUE)))
meta$label <- with(meta, sprintf("%s (%.0fk; >=97%%: %.0fk)", phylum, n_total/1000, n_ge97/1000))
meta$label <- factor(meta$label, levels = meta$label)
meta$low_ge97 <- meta$n_ge97 < 500

d <- d %>% left_join(meta, by = "phylum")
d$label <- factor(d$label, levels = levels(meta$label))

# --- Fig 1: 70-97%, all ASVs as denominator ----
d1 <- d %>% filter(identity >= 70 & identity <= 97)
p1 <- ggplot(d1, aes(x = identity, y = frac_ge, colour = label)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = c(75, 80, 85, 90, 95), linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  annotate("text", x = c(75.3, 80.3, 85.3, 90.3, 95.3, 97.3), y = 1.04,
           label = c("phylum~75", "class~80", "order~85", "family~90", "genus~95", "97%"),
           size = 2.3, colour = "grey30", hjust = 0, vjust = 1) +
  scale_x_continuous(breaks = seq(70, 97, 3), limits = c(70, 97.6)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.08),
                     name = "Cumulative fraction of all phylum ASVs (identity >= x)") +
  scale_colour_manual(values = hue_pal()(15), name = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 7), panel.grid.minor = element_blank())

fig1 <- file.path(FIGOUT, "fig_v3_homd_phylum_identity_curve_70_97")
ggsave(paste0(fig1, ".png"), p1, width = 9.5, height = 5.5, dpi = 300)
ggsave(paste0(fig1, ".pdf"), p1, width = 9.5, height = 5.5)
message("[13h] fig1 (70-97) written")

# --- Fig 2: 97-100%, denominator = ASVs >=97% ----
# frac_ge at x=97 equals n_ge97/n_total (the max within x>=97 window);
# normalizing by it makes every curve start at 100%.
d2 <- d %>% filter(identity >= 97)   # inherits label / low_ge97 from meta join
d2 <- d2 %>% group_by(phylum) %>%
  mutate(frac2 = frac_ge / max(frac_ge)) %>% ungroup()

p2 <- ggplot(d2, aes(x = identity, y = frac2, colour = label)) +
  geom_line(data = d2 %>% filter(!low_ge97), linewidth = 0.9) +
  geom_line(data = d2 %>% filter(low_ge97), linewidth = 0.5, linetype = "dashed", alpha = 0.55) +
  geom_vline(xintercept = 99, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 97, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  annotate("text", x = 99.1, y = 1.02, label = "99%", size = 2.5, colour = "grey30", hjust = 0) +
  scale_x_continuous(breaks = seq(97, 100, 0.5), limits = c(97, 100)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.05),
                     name = "Cumulative fraction of ASVs >=97% (identity >= x)") +
  scale_colour_manual(values = hue_pal()(15), name = NULL) +
  labs(caption = "Dashed grey lines: phyla with <500 ASVs >=97% identity (unstable curves)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 7), panel.grid.minor = element_blank())

fig2 <- file.path(FIGOUT, "fig_v3_homd_phylum_identity_curve_97_100")
ggsave(paste0(fig2, ".png"), p2, width = 9.5, height = 5.5, dpi = 300)
ggsave(paste0(fig2, ".pdf"), p2, width = 9.5, height = 5.5)
message("[13h] fig2 (97-100) written")

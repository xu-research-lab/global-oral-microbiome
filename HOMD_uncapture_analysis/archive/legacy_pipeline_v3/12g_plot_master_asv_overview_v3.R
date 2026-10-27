#!/usr/bin/env Rscript
# pipeline_v3: Overview of 07_master_asv_table.tsv
# Panels:
#   A. ASV subject prevalence distribution (log10 x), split by contamination flag
#   B. Cluster size distribution (log10 x)
#   C. ASV counts by contamination flag type (incl. organelle)
#   D. Number of clusters by annotation level
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

PIPE_V3 <- "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3"
FIGOUT  <- file.path(PIPE_V3, "analysis_report", "figs")
dir.create(FIGOUT, showWarnings = FALSE, recursive = TRUE)

message("[12g_v3] loading 07_master_asv_table.tsv")
d <- read.delim(file.path(PIPE_V3, "07_master_asv_table.tsv"),
                sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))
d$subject_prevalence <- num(d$subject_prevalence)
d$cluster_size      <- num(d$cluster_size)
d$cluster_n_asvs    <- num(d$cluster_n_asvs)
d$contamination_flag <- ifelse(is.na(d$contamination_flag), "False", d$contamination_flag)
d$organelle_flag     <- ifelse(is.na(d$organelle_flag), "False", d$organelle_flag)

# ---- cluster-level table ----
cl <- d %>% distinct(cluster_id, .keep_all = TRUE)
cl$annotation_level[is.na(cl$annotation_level) | cl$annotation_level == ""] <- "unannotated"

# ---- Panel A: ASV subject prevalence ----
dA <- d %>% filter(!is.na(subject_prevalence))
pA <- ggplot(dA, aes(x = subject_prevalence + 0.5, fill = contamination_flag)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.65, colour = NA) +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE, big.mark = ",")) +
  scale_fill_manual(values = c("False" = "#4DBBD5", "True" = "#E64B35"),
                    name = "contamination flag") +
  labs(x = "subject prevalence (log10, +0.5)", y = "No. ASVs",
       title = "A. ASV subject prevalence") +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.75, 0.8), legend.background = element_rect(fill = "white", colour = "grey80"))

# ---- Panel B: cluster size ----
dB <- cl %>% filter(!is.na(cluster_size))
pB <- ggplot(dB, aes(x = cluster_size + 0.5)) +
  geom_histogram(bins = 50, fill = "#3C5488", colour = "white", alpha = 0.9) +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE, big.mark = ",")) +
  labs(x = "cluster size: No. derep ASVs (log10, +0.5)", y = "No. clusters",
       title = "B. Cluster size distribution") +
  theme_bw(base_size = 11)

# ---- Panel C: contamination flag composition (ASV level) ----
flag_order <- c("clean", "cluster97", "genus_list", "both",
                "clean;organelle", "cluster97;organelle", "genus_list;organelle", "both;organelle")
dC <- d %>%
  mutate(flag_type = ifelse(is.na(flag_type) | flag_type == "", "clean", flag_type)) %>%
  count(flag_type) %>%
  mutate(flag_type = factor(flag_type, levels = flag_order)) %>%
  filter(!is.na(flag_type))
dC$label <- paste0(dC$n, " (", round(100 * dC$n / sum(dC$n), 1), "%)")
pC <- ggplot(dC, aes(x = flag_type, y = n, fill = flag_type)) +
  geom_col(colour = "white", width = 0.7) +
  geom_text(aes(label = label), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("clean" = "#4DBBD5", "cluster97" = "#E64B35",
                               "genus_list" = "#F39B7F", "both" = "#B09C85",
                               "clean;organelle" = "#8491B4", "cluster97;organelle" = "#7E6148",
                               "genus_list;organelle" = "#91D1C2", "both;organelle" = "#00A087"),
                    guide = "none") +
  labs(x = NULL, y = "No. ASVs", title = "C. ASVs by contamination flag type") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# ---- Panel D: clusters by annotation level ----
lv_order <- c("species", "genus", "family", "order", "class", "phylum", "domain", "unannotated")
dD <- cl %>% count(annotation_level) %>%
  mutate(annotation_level = factor(annotation_level, levels = lv_order)) %>%
  filter(!is.na(annotation_level))
dD$label <- paste0(dD$n, " (", round(100 * dD$n / sum(dD$n), 1), "%)")
pD <- ggplot(dD, aes(x = annotation_level, y = n, fill = annotation_level)) +
  geom_col(colour = "white", width = 0.7) +
  geom_text(aes(label = label), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("species" = "#E64B35", "genus" = "#4DBBD5",
                               "family" = "#00A087", "order" = "#3C5488",
                               "class" = "#F39B7F", "phylum" = "#8491B4",
                               "domain" = "#B09C85", "unannotated" = "grey70"),
                    guide = "none") +
  labs(x = NULL, y = "No. clusters", title = "D. Clusters by annotation level") +
  theme_bw(base_size = 11)

# ---- combine ----
fig <- (pA | pB) / (pC | pD) + plot_layout(widths = c(1, 1), heights = c(1, 1))
ggsave(file.path(FIGOUT, "fig_v3_master_asv_overview.png"), fig, width = 12, height = 10, dpi = 300)
ggsave(file.path(FIGOUT, "fig_v3_master_asv_overview.pdf"), fig, width = 12, height = 10)
message("[12g_v3] saved fig_v3_master_asv_overview.png/pdf")

# ---- text summary for the record ----
cat("\n===== summary =====\n")
cat("ASVs:", nrow(d), "| clusters:", nrow(cl), "\n")
cat("contamination_flag: True =", sum(d$contamination_flag == "True"),
    "| False =", sum(d$contamination_flag == "False"), "\n")
cat("organelle_flag: True =", sum(d$organelle_flag == "True"),
    "| False =", sum(d$organelle_flag == "False"), "\n")
cat("cluster annotation levels:\n")
print(dD %>% arrange(annotation_level))

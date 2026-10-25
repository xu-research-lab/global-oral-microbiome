# Sample distributions for the 10 selected 6-36-month genera in the HC saliva audit.
# Run from the project root; no statistical inference or DA refitting is performed.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})
source("differential_abundance/age_differential_ova/da_helpers.R")
source("differential_abundance/ova_plot_helpers.R")

base <- "differential_abundance/result/age_hc_ova/age_da_consensus"
fig_dir <- file.path(base, "figures_saliva_target_top10_globalpeak_maaslin2coef")
out_dir <- file.path(base, "diagnostics_infant_age_project")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
selected <- fread(file.path(fig_dir, "figureB_selected_2method.tsv"), data.table = FALSE)
selected <- selected[selected$target == "6-36 months", ]
features <- selected$feature
genera <- sub(".*g__", "", features)
stopifnot(length(features) == 10L, !anyDuplicated(features), !anyDuplicated(genera))

meta <- load_age_metadata(hc_only = TRUE)
count_ids <- strsplit(readLines(COUNTS_TSV, n = 1), "\t", fixed = TRUE)[[1]][-1]
meta <- meta[meta$sample_id %in% count_ids & meta$age_cleaned_group %in% AGE_ORDER[1:2], ]
stopifnot(!anyDuplicated(meta$sample_id))
pattern_file <- tempfile()
rows_file <- tempfile()
writeLines(paste0(features, "\t"), pattern_file)
status <- system2("rg", c("-F", "-f", shQuote(pattern_file), shQuote(RELABUND_TSV)),
                  stdout = rows_file)
stopifnot(status == 0L)
rel <- fread(rows_file, header = FALSE, data.table = FALSE)
names(rel) <- c("feature", strsplit(readLines(RELABUND_TSV, n = 1), "\t", fixed = TRUE)[[1]][-1])
unlink(c(pattern_file, rows_file))
stopifnot(setequal(rel$feature, features), all(meta$sample_id %in% names(rel)))
x <- t(as.matrix(rel[, meta$sample_id, drop = FALSE]))
colnames(x) <- rel$feature
x <- x[, features, drop = FALSE]
stopifnot(identical(rownames(x), meta$sample_id), all(is.finite(x)), all(x >= 0 & x <= 1))

samples <- bind_rows(lapply(seq_along(features), function(i) {
  data.frame(sample_id = meta$sample_id, bioproject_id = meta$bioproject_id,
             age_group = meta$age_cleaned_group, genus = genera[i],
             abundance_pct = 100 * x[, i], stringsAsFactors = FALSE)
}))
plot_data <- samples
plot_data$age_group <- factor(plot_data$age_group, levels = AGE_ORDER[1:2])
plot_data$genus <- factor(plot_data$genus, levels = rev(genera))
plot_data$row_y <- as.numeric(plot_data$genus)
plot_data$box_y <- plot_data$row_y + ifelse(plot_data$age_group == AGE_ORDER[1], 0.18, -0.18)
project_counts <- meta %>% count(bioproject_id, name = "n_samples") %>%
  arrange(desc(n_samples), bioproject_id)
project_order <- project_counts$bioproject_id
plot_data$bioproject_id <- factor(plot_data$bioproject_id, levels = project_order)

# Compute quartiles and Tukey whiskers on untransformed relative abundance.
# Transform coordinates only after computing the summaries, preserving the
# standard 1.5-IQR boxplot definition despite the nonlinear displayed axis.
box_stats <- plot_data %>% group_by(genus, age_group, box_y) %>%
  group_modify(~ {
    v <- .x$abundance_pct
    q <- quantile(v, c(0.25, 0.5, 0.75), names = FALSE, type = 7)
    iqr <- q[3] - q[1]
    data.frame(n_samples = length(v), n_zero = sum(v == 0),
               prevalence_pct = 100 * mean(v > 0), mean_pct = mean(v),
               q1 = q[1], median = q[2], q3 = q[3],
               whisker_low = min(v[v >= q[1] - 1.5 * iqr]),
               whisker_high = max(v[v <= q[3] + 1.5 * iqr]))
  }) %>% ungroup()
expected_n <- c(786, 543)
got_n <- box_stats %>% distinct(age_group, n_samples) %>% arrange(age_group)
stopifnot(identical(as.integer(got_n$n_samples), as.integer(expected_n)),
          nrow(box_stats) == 20L, nrow(plot_data) == 13290L, length(project_order) == 12L,
          !anyDuplicated(plot_data[c("sample_id", "genus")]),
          sum(plot_data$abundance_pct == 0) == sum(box_stats$n_zero))
# Independent check against the previously calculated pooled age-group means.
audit_path <- file.path(out_dir, "descriptive_exclusion_sensitivity.tsv")
if (file.exists(audit_path)) {
  audit <- fread(audit_path, data.table = FALSE)
  for (age in AGE_ORDER[1:2]) {
    ref_col <- if (age == AGE_ORDER[1]) "mean_0_6_pct" else "mean_6_36_pct"
    obs <- box_stats %>% filter(age_group == age)
    ref <- audit[audit$scenario == "all_samples", ]
    stopifnot(max(abs(obs$mean_pct - ref[[ref_col]][match(obs$genus, ref$genus)])) < 1e-10)
  }
}
fwrite(box_stats, file.path(out_dir, "infant_project_boxplot_statistics.tsv"), sep = "\t")
fwrite(plot_data[, c("genus", "sample_id", "bioproject_id", "age_group", "abundance_pct")],
       file.path(out_dir, "infant_project_boxplot_samples.tsv"), sep = "\t")

project_colors <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                             "#E69F00", "#56B4E9", "#8C564B", "#6A3D9A",
                             "#A6A326", "#E377C2", "#17BECF", "#555555"), project_order)
age_fills <- c("0-6 months" = "#FFFFFF", "6-36 months" = "#BDBDBD")
p <- ggplot() +
  geom_hline(yintercept = seq(0.5, 10.5, 1), colour = "grey90", linewidth = 0.2) +
  geom_boxplot(data = box_stats,
               aes(y = box_y, xmin = whisker_low, xlower = q1, xmiddle = median,
                   xupper = q3, xmax = whisker_high, fill = age_group,
                   group = interaction(genus, age_group)),
               stat = "identity", orientation = "y", width = 0.26,
               colour = "grey15", linewidth = 0.3,
               staplewidth = 0.5, outlier.shape = NA) +
  geom_point(data = plot_data,
             aes(x = abundance_pct, y = box_y, colour = bioproject_id),
             position = position_jitter(width = 0, height = 0.085, seed = 20261006),
             size = 0.55, alpha = 0.5, shape = 16) +
  geom_segment(data = box_stats,
                aes(x = median, xend = median, y = box_y - 0.13, yend = box_y + 0.13),
                colour = "grey15", linewidth = 0.45) +
  scale_y_continuous(breaks = seq_along(genera), labels = rev(genera),
                     limits = c(0.5, 10.5), expand = c(0, 0)) +
  scale_x_continuous(transform = scales::pseudo_log_trans(sigma = 0.001, base = 10),
                     breaks = c(0, 0.01, 0.1, 1, 10, 100),
                     labels = c("0", "0.01", "0.1", "1", "10", "100"),
                     limits = c(0, 100), expand = expansion(mult = c(0.045, 0.025))) +
  scale_fill_manual(values = age_fills, name = "Recorded age group",
                     guide = guide_legend(order = 1)) +
  scale_colour_manual(values = project_colors, name = "BioProject",
                       guide = guide_legend(order = 2, ncol = 1,
                                            override.aes = list(size = 2, alpha = 1))) +
  labs(x = "Relative abundance (%)", y = NULL,
       caption = paste(
         "Boxes: median and interquartile range; whiskers: 1.5 x IQR. Points: individual samples; zero values retained.",
         "Upper / lower box per genus: recorded 0-6 months (n = 786) / 6-36 months (n = 543). Axis uses a pseudo-log scale.",
         "PRJNA997934 records integer years (0Y cannot establish 0-6 months). PRJEB35824 includes repeated visits.",
         sep = "\n")) +
  theme_bw(base_family = OVA_GENUS_FONT, base_size = 8) +
  theme(axis.text.y = element_text(face = "italic", size = 8, colour = "black"),
        axis.text.x = element_text(size = 7, colour = "black"),
        axis.ticks.y = element_blank(), axis.title.x = element_text(size = 9),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.3),
        legend.position = "right", legend.direction = "vertical", legend.box = "vertical",
        legend.box.just = "left", legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(3.5, "mm"),
        plot.caption.position = "plot",
        plot.caption = element_text(size = 7, hjust = 0, margin = margin(t = 7)),
        plot.margin = margin(5, 5, 5, 5))
stem <- file.path(out_dir, "figure_infant_age_project_boxplots")
ggsave(paste0(stem, ".pdf"), p, device = cairo_pdf,
       width = 200, height = 165, units = "mm", bg = "white")
ggsave(paste0(stem, ".png"), p, device = ragg::agg_png,
       width = 200, height = 165, units = "mm", dpi = 300, bg = "white")
cat("Saved one-panel boxplots for", length(features), "genera with", length(project_order),
    "project colors; all", nrow(plot_data),
    "sample observations retained, including", sum(plot_data$abundance_pct == 0), "zeros.\n")
cat("Sample counts and project means verified; output:", stem, "\n")

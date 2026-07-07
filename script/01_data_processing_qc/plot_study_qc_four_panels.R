suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(aplot)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default) {
  key <- paste0("--", name, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^", key), "", hit[[1]])
}

paths <- list(
  asv = get_arg("asv", "/beegfs/dongbiao/oral/alpha_diversity_plots/tables/observed_asv_sample_values.tsv"),
  chimera = get_arg("chimera", "QC/oral_chim_perc_YL_260519.csv"),
  keeping = get_arg("keeping", "/beegfs/dongbiao/oral/study_close_reference_keeping_rate_cleaned_all.csv"),
  summary = get_arg("summary", "/beegfs/dongbiao/oral/merge_closed_reference_tables_cleaned_all_summary.csv"),
  depth = get_arg("depth", "QC/depth_df.csv"),
  outdir = get_arg("outdir", "QC"),
  prefix = get_arg("prefix", "study_qc_four_panels_asv_cleaned_all")
)

asv_threshold <- as.numeric(get_arg("asv-threshold", "10"))
low_asv_ratio_threshold <- as.numeric(get_arg("low-asv-ratio-threshold", "0.25"))
chimera_threshold <- as.numeric(get_arg("chimera-threshold", "25"))
study_chimera_ratio_threshold <- as.numeric(get_arg("study-chimera-ratio-threshold", "0.30"))
read_retention_threshold <- as.numeric(get_arg("read-retention-threshold", "0.80"))

dir.create(paths$outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading study summary...")
study_summary <- fread(paths$summary)
setnames(
  study_summary,
  old = c("study_id", "sample_count_in_closed_reference_table", "sample_count_in_metadata"),
  new = c("study_id", "closed_reference_sample_count", "metadata_sample_count")
)

message("Reading closed-reference keeping rate...")
keeping <- fread(paths$keeping)
setnames(keeping, "ratio", "closed_reference_reads_rate")

message("Reading raw read depth table...")
depth <- fread(paths$depth)
stopifnot(all(c("project_id", "reads_depth") %in% names(depth)))
setnames(depth, old = c("project_id", "reads_depth"), new = c("study_id", "raw_reads"))
depth[, study_id := as.character(study_id)]
depth[, raw_reads := as.numeric(raw_reads)]
depth_qc <- depth[, .(
  read_depth_sample_count = .N,
  mean_raw_reads = mean(raw_reads, na.rm = TRUE),
  median_raw_reads = median(raw_reads, na.rm = TRUE)
), by = study_id]

message("Reading observed ASV table...")
asv <- fread(paths$asv)
stopifnot(all(c("study_id", "sample_id", "observed_asv") %in% names(asv)))
asv[, study_id := as.character(study_id)]
asv[, sample_id := as.character(sample_id)]
asv[, observed_asv := as.numeric(observed_asv)]
asv[, low_asv_sample := observed_asv <= asv_threshold]

asv_qc <- asv[, .(
  asv_sample_count = .N,
  low_asv_sample_count = sum(low_asv_sample, na.rm = TRUE),
  low_asv_sample_ratio = sum(low_asv_sample, na.rm = TRUE) / .N,
  observed_asv_median = as.numeric(median(observed_asv, na.rm = TRUE))
), by = study_id]

fwrite(
  asv[, .(study_id, sample_id, observed_asv, low_asv_sample)],
  file.path(paths$outdir, paste0(paths$prefix, "_sample_asv_values.csv"))
)

message("Reading chimera table...")
chimera <- fread(paths$chimera)
stopifnot(all(c("Project", "SampleID", "chim_perc") %in% names(chimera)))
chimera[, study_id := as.character(Project)]
chimera[, chim_perc := as.numeric(chim_perc)]
chimera_qc <- chimera[, .(
  chimera_sample_count = .N,
  high_chimera_sample_count = sum(chim_perc > chimera_threshold, na.rm = TRUE),
  high_chimera_sample_ratio = sum(chim_perc > chimera_threshold, na.rm = TRUE) / .N
), by = study_id]
chimera_qc[, high_chimera_flag := high_chimera_sample_ratio > study_chimera_ratio_threshold]

metrics <- Reduce(
  function(x, y) merge(x, y, by = "study_id", all = TRUE),
  list(study_summary, keeping, depth_qc, asv_qc, chimera_qc)
)

metrics[, sample_retention_ratio := closed_reference_sample_count / metadata_sample_count]
metrics[, study_id := as.character(study_id)]
metrics[, low_asv_flag := low_asv_sample_ratio > low_asv_ratio_threshold]
metrics[, low_reads_flag := closed_reference_reads_rate < read_retention_threshold]
metrics[, sort_sample_count := ifelse(is.na(closed_reference_sample_count), -1, closed_reference_sample_count)]

study_levels <- metrics[order(-sort_sample_count, study_id), study_id]
metrics[, study_id := factor(study_id, levels = study_levels)]

metrics_out <- copy(metrics)
metrics_out[, study_id := as.character(study_id)]
metrics_out[, sort_sample_count := NULL]
fwrite(metrics_out, file.path(paths$outdir, paste0(paths$prefix, "_study_metrics.csv")))
fwrite(
  metrics_out[low_asv_flag == TRUE],
  file.path(paths$outdir, paste0(paths$prefix, "_low_asv_flagged_studies.csv"))
)
fwrite(
  metrics_out[low_reads_flag == TRUE],
  file.path(paths$outdir, paste0(paths$prefix, "_low_reads_flagged_studies.csv"))
)
fwrite(
  metrics_out[high_chimera_flag == TRUE],
  file.path(paths$outdir, paste0(paths$prefix, "_high_chimera_flagged_studies.csv"))
)

theme_qc <- function(show_x = FALSE) {
  theme_bw(base_size = 18) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = if (show_x) element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10) else element_blank(),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.title.x = if (show_x) element_text(size = 18, face = "bold") else element_blank(),
      plot.title = element_text(face = "bold", size = 20),
      plot.margin = margin(5, 10, 5, 10)
    )
}

p0 <- ggplot(metrics, aes(study_id, mean_raw_reads)) +
  geom_col(fill = "#6B6ECF", width = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "#D62728", linewidth = 0.9) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Mean raw reads per study (log10 scale; red line: 1,000 reads)",
    x = NULL,
    y = "Mean reads"
  ) +
  theme_qc(FALSE)

p1 <- ggplot(metrics, aes(study_id, low_asv_sample_ratio, fill = low_asv_flag)) +
  geom_col(width = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = low_asv_ratio_threshold, linetype = "dashed", color = "#D62728", linewidth = 0.9) +
  scale_fill_manual(values = c("FALSE" = "#4C78A8", "TRUE" = "#D62728"), na.value = "grey80", guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = paste0("Samples with observed ASV <= ", asv_threshold, " (red: study ratio > ", low_asv_ratio_threshold * 100, "%)"),
    x = NULL,
    y = "Sample ratio"
  ) +
  theme_qc(FALSE)

p2 <- ggplot(metrics, aes(study_id, closed_reference_reads_rate, fill = low_reads_flag)) +
  geom_col(width = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = 0.70, linetype = "dotted", color = "#D62728", linewidth = 0.9) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "#D62728", linewidth = 0.9) +
  scale_fill_manual(values = c("FALSE" = "#59A14F", "TRUE" = "#D62728"), na.value = "grey80", guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Closed-reference retained reads rate (red bars: <80%; dotted/dashed lines: 70%/80%)",
    x = NULL,
    y = "Read ratio"
  ) +
  theme_qc(FALSE)

p3 <- ggplot(metrics, aes(study_id, high_chimera_sample_ratio, fill = high_chimera_flag)) +
  geom_col(width = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = study_chimera_ratio_threshold, linetype = "dashed", color = "#D62728", linewidth = 0.9) +
  scale_fill_manual(values = c("FALSE" = "#9C755F", "TRUE" = "#D62728"), na.value = "grey80", guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = paste0("Samples with chimera percent > ", chimera_threshold, "% (red: study ratio > ", study_chimera_ratio_threshold * 100, "%)"),
    x = NULL,
    y = "Sample ratio"
  ) +
  theme_qc(FALSE)

p4 <- ggplot(metrics, aes(study_id, closed_reference_sample_count)) +
  geom_col(fill = "#F28E2B", width = 0.85, na.rm = TRUE) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Retained samples in closed-reference table (log10 scale; studies sorted high to low)",
    x = "Study",
    y = "Samples (log10)"
  ) +
  theme_qc(TRUE)

combined <- aplot::plot_list(
  gglist = list(p0, p1, p2, p3, p4),
  ncol = 1,
  heights = c(1, 1, 1, 1, 1.25),
  output = "patchwork"
)

pdf_file <- file.path(paths$outdir, paste0(paths$prefix, ".pdf"))
png_file <- file.path(paths$outdir, paste0(paths$prefix, ".png"))

ggplot2::ggsave(pdf_file, combined, width = 30, height = 14, limitsize = FALSE)
ggplot2::ggsave(png_file, combined, width = 30, height = 14, dpi = 300, limitsize = FALSE)

message("Done.")
message("Observed ASV samples: ", nrow(asv))
message("Study metrics: ", file.path(paths$outdir, paste0(paths$prefix, "_study_metrics.csv")))
message("Flagged low ASV studies: ", file.path(paths$outdir, paste0(paths$prefix, "_low_asv_flagged_studies.csv")))
message("Flagged low reads studies: ", file.path(paths$outdir, paste0(paths$prefix, "_low_reads_flagged_studies.csv")))
message("Flagged chimera studies: ", file.path(paths$outdir, paste0(paths$prefix, "_high_chimera_flagged_studies.csv")))
message("Plot PDF: ", pdf_file)
message("Plot PNG: ", png_file)

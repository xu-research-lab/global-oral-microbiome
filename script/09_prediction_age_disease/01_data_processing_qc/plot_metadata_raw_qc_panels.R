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
  input = get_arg("input", "metadata_all_for_qc_raw.tsv"),
  outdir = get_arg("outdir", "QC"),
  prefix = get_arg("prefix", "metadata_raw_qc_panels")
)

raw_reads_threshold <- as.numeric(get_arg("raw-reads-threshold", "1000"))
cr_mapping_threshold <- as.numeric(get_arg("cr-mapping-threshold", "0.30"))
chimera_threshold <- as.numeric(get_arg("chimera-threshold", "25"))
asv_threshold <- as.numeric(get_arg("asv-threshold", "30"))
low_asv_ratio_threshold <- as.numeric(get_arg("low-asv-ratio-threshold", "0.50"))
study_chimera_ratio_threshold <- as.numeric(get_arg("study-chimera-ratio-threshold", "0.30"))

dir.create(paths$outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading metadata input...")
dt <- fread(paths$input)
required_cols <- c("bioproject_id", "sample_id", "total_reads", "chim_perc", "asv_count", "Observed", "subject_id", "study_cr_mapping_rate", "sample_cr_mapping_rate")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing columns in input table: ", paste(missing_cols, collapse = ", "))
}

dt[, study_id := as.character(bioproject_id)]
dt[, sample_id := as.character(sample_id)]
dt[, raw_reads := suppressWarnings(as.numeric(total_reads))]
dt[, chim_perc_num := suppressWarnings(as.numeric(chim_perc))]
dt[, asv_count_num := suppressWarnings(as.numeric(asv_count))]
dt[, observed_num := suppressWarnings(as.numeric(Observed))]
dt[, study_cr_mapping_rate_num := suppressWarnings(as.numeric(study_cr_mapping_rate))]
dt[, sample_cr_mapping_rate_num := suppressWarnings(as.numeric(sample_cr_mapping_rate))]
dt[, dechimera_reads := raw_reads * (1 - chim_perc_num / 100)]
dt[is.na(raw_reads) | is.na(chim_perc_num), dechimera_reads := NA_real_]
dt[, low_asv_sample := asv_count_num <= asv_threshold]
dt[, high_chimera_sample := chim_perc_num > chimera_threshold]
dt[, low_cr_sample := sample_cr_mapping_rate_num < cr_mapping_threshold]
dt[, observed_present := !is.na(observed_num)]

dt[, is_oral := site_group == "oral"]

# Overall study metrics (all samples, for display)
study_metrics <- dt[, .(
  metadata_sample_count = .N,
  oral_sample_count = sum(is_oral, na.rm = TRUE),
  unique_subject_count = uniqueN(fifelse(is.na(subject_id) | subject_id == "" | subject_id == "NA", NA_character_, as.character(subject_id)), na.rm = TRUE),
  mean_raw_reads = mean(raw_reads, na.rm = TRUE),
  median_raw_reads = median(raw_reads, na.rm = TRUE),
  observed_present_count = sum(observed_present, na.rm = TRUE),
  study_cr_mapping_rate = suppressWarnings(first(na.omit(study_cr_mapping_rate_num)))
), by = study_id]

# Oral-only metrics for QC flag computation
oral_metrics <- dt[is_oral == TRUE, .(
  oral_low_asv_count = sum(low_asv_sample, na.rm = TRUE),
  oral_high_chimera_count = sum(high_chimera_sample, na.rm = TRUE),
  oral_low_cr_count = sum(low_cr_sample, na.rm = TRUE)
), by = study_id]

study_metrics <- merge(study_metrics, oral_metrics, by = "study_id", all.x = TRUE)
study_metrics[is.na(oral_low_asv_count), oral_low_asv_count := 0]
study_metrics[is.na(oral_high_chimera_count), oral_high_chimera_count := 0]
study_metrics[is.na(oral_low_cr_count), oral_low_cr_count := 0]
study_metrics[is.infinite(study_cr_mapping_rate), study_cr_mapping_rate := NA_real_]
study_metrics[, additional_sample_count := pmax(metadata_sample_count - unique_subject_count, 0)]

# Ratios based on ORAL samples only (for QC flagging)
study_metrics[, low_asv_sample_ratio := fifelse(oral_sample_count > 0, oral_low_asv_count / oral_sample_count, 0)]
study_metrics[, high_chimera_sample_ratio := fifelse(oral_sample_count > 0, oral_high_chimera_count / oral_sample_count, 0)]
study_metrics[, observed_present_ratio := observed_present_count / metadata_sample_count]
study_metrics[, low_asv_flag := low_asv_sample_ratio > low_asv_ratio_threshold]
study_metrics[, high_chimera_flag := high_chimera_sample_ratio > study_chimera_ratio_threshold]
study_metrics[, low_cr_mapping_flag := study_cr_mapping_rate < cr_mapping_threshold]
study_metrics[, sort_sample_count := metadata_sample_count]
setorder(study_metrics, -sort_sample_count, study_id)
study_levels <- study_metrics$study_id
study_metrics[, study_id := factor(study_id, levels = study_levels)]

study_metrics_out <- copy(study_metrics)
study_metrics_out[, study_id := as.character(study_id)]
study_metrics_out[, sort_sample_count := NULL]
fwrite(study_metrics_out, file.path(paths$outdir, paste0(paths$prefix, "_study_metrics.csv")))
fwrite(study_metrics_out[low_asv_flag == TRUE], file.path(paths$outdir, paste0(paths$prefix, "_low_asv_flagged_studies.csv")))
fwrite(study_metrics_out[high_chimera_flag == TRUE], file.path(paths$outdir, paste0(paths$prefix, "_high_chimera_flagged_studies.csv")))
fwrite(study_metrics_out[low_cr_mapping_flag == TRUE], file.path(paths$outdir, paste0(paths$prefix, "_low_cr_mapping_studies.csv")))

theme_qc <- function(show_x = FALSE) {
  theme_bw(base_size = 18) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = if (show_x) element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8) else element_blank(),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.title.x = if (show_x) element_text(size = 18, face = "bold") else element_blank(),
      plot.title = element_text(face = "bold", size = 20),
      plot.margin = margin(5, 10, 5, 10)
    )
}

p0 <- ggplot(study_metrics, aes(study_id, mean_raw_reads)) +
  geom_col(fill = "#6B6ECF", width = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = raw_reads_threshold, linetype = "dashed", color = "#D62728", linewidth = 0.9) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Mean raw reads per study (log10 scale; red line: 1,000 reads)",
    x = NULL,
    y = "Mean reads"
  ) +
  theme_qc(FALSE)

p1 <- ggplot(study_metrics, aes(study_id, low_asv_sample_ratio, fill = low_asv_flag)) +
  geom_col(width = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = low_asv_ratio_threshold, linetype = "dashed", color = "#D62728", linewidth = 0.9) +
  scale_fill_manual(values = c("FALSE" = "#4C78A8", "TRUE" = "#D62728"), na.value = "grey80", guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Samples with ASV count <= 30 (red: study ratio > 50%)",
    x = NULL,
    y = "Sample ratio"
  ) +
  theme_qc(FALSE)

p2 <- ggplot(study_metrics, aes(study_id, study_cr_mapping_rate, fill = low_cr_mapping_flag)) +
  geom_col(width = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = cr_mapping_threshold, linetype = "dashed", color = "#7F7F7F", linewidth = 0.9) +
  scale_fill_manual(values = c("FALSE" = "#59A14F", "TRUE" = "#A0A0A0"), na.value = "grey80", guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = paste0("Closed-reference mapping rate (grey line: ", cr_mapping_threshold * 100, "%)"),
    x = NULL,
    y = "CR mapping rate"
  ) +
  theme_qc(FALSE)

p3 <- ggplot(study_metrics, aes(study_id, high_chimera_sample_ratio, fill = high_chimera_flag)) +
  geom_col(width = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = study_chimera_ratio_threshold, linetype = "dashed", color = "#D62728", linewidth = 0.9) +
  scale_fill_manual(values = c("FALSE" = "#9C755F", "TRUE" = "#D62728"), na.value = "grey80", guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Samples with chimera percent > 25% (red: study ratio > 30%)",
    x = NULL,
    y = "Sample ratio"
  ) +
  theme_qc(FALSE)

p4 <- ggplot(study_metrics, aes(study_id)) +
  geom_col(aes(y = metadata_sample_count, fill = "Total samples"), width = 0.85, na.rm = TRUE) +
  geom_col(aes(y = unique_subject_count, fill = "Unique subjects"), width = 0.45, na.rm = TRUE) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = c("Total samples" = "#F4B183", "Unique subjects" = "#C55A11")) +
  labs(
    title = "Samples in metadata_all_for_qc_raw.tsv (total samples with unique subjects inset)",
    x = "Study",
    y = "Samples",
    fill = NULL
  ) +
  theme_qc(TRUE) +
  theme(
    legend.position = c(0.87, 0.86),
    legend.background = element_rect(fill = alpha("white", 0.85), color = "grey70"),
    legend.key.size = unit(0.5, "cm")
  )

combined_study <- aplot::plot_list(
  gglist = list(p0, p1, p2, p3, p4),
  ncol = 1,
  heights = c(1, 1, 1, 1, 1.25),
  output = "patchwork"
)

study_pdf <- file.path(paths$outdir, paste0(paths$prefix, "_study_level.pdf"))
study_png <- file.path(paths$outdir, paste0(paths$prefix, "_study_level.png"))
ggsave(study_pdf, combined_study, width = 30, height = 16, limitsize = FALSE)
ggsave(study_png, combined_study, width = 30, height = 16, dpi = 300, limitsize = FALSE)

study_order <- study_metrics_out[, .(study_id, metadata_sample_count)]
setorder(study_order, -metadata_sample_count, study_id)
study_order[, study_rank := .I]
palette_cycle <- c(
  "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#B07AA1",
  "#EDC948", "#9C755F", "#BAB0AC", "#1F77B4", "#FF7F0E", "#2CA02C"
)
study_order[, study_color := palette_cycle[((study_rank - 1) %% length(palette_cycle)) + 1]]
fwrite(study_order, file.path(paths$outdir, paste0(paths$prefix, "_sample_study_order.csv")))

sample_dt <- merge(dt, study_order, by = "study_id", all.x = TRUE)
setorder(sample_dt, study_rank, -raw_reads, sample_id)
sample_dt[, sample_rank := .I]

study_boundaries <- sample_dt[, .(
  start_rank = min(sample_rank),
  end_rank = max(sample_rank)
), by = .(study_id, study_rank)][order(study_rank)]
study_boundaries[, boundary_rank := end_rank + 0.5]

fwrite(
  sample_dt[, .(
    study_id, sample_id, sample_rank, raw_reads, dechimera_reads,
    chim_perc = chim_perc_num, asv_count = asv_count_num, Observed = observed_num,
    metadata_sample_count
  )],
  file.path(paths$outdir, paste0(paths$prefix, "_sample_ordered_values.csv"))
)

threshold_summary <- data.table(
  metric = c("raw_reads", "dechimera_reads", "chim_perc", "asv_count", "Observed_na"),
  filter_rule = c(
    paste0("raw_reads < ", raw_reads_threshold),
    paste0("dechimera_reads < ", raw_reads_threshold),
    paste0("chim_perc > ", chimera_threshold),
    paste0("asv_count <= ", asv_threshold),
    "Observed is NA"
  ),
  threshold = c(raw_reads_threshold, raw_reads_threshold, chimera_threshold, asv_threshold, NA_real_),
  failed_samples = c(
    sum(sample_dt$raw_reads < raw_reads_threshold, na.rm = TRUE),
    sum(sample_dt$dechimera_reads < raw_reads_threshold, na.rm = TRUE),
    sum(sample_dt$chim_perc_num > chimera_threshold, na.rm = TRUE),
    sum(sample_dt$asv_count_num <= asv_threshold, na.rm = TRUE),
    sum(is.na(sample_dt$observed_num))
  ),
  non_missing_samples = c(
    sum(!is.na(sample_dt$raw_reads)),
    sum(!is.na(sample_dt$dechimera_reads)),
    sum(!is.na(sample_dt$chim_perc_num)),
    sum(!is.na(sample_dt$asv_count_num)),
    nrow(sample_dt)
  )
)
threshold_summary[, failed_ratio := failed_samples / non_missing_samples]
fwrite(threshold_summary, file.path(paths$outdir, paste0(paths$prefix, "_sample_threshold_summary.csv")))

theme_sample_qc <- function() {
  theme_bw(base_size = 18) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 17),
      axis.title.y = element_text(size = 19, face = "bold"),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold", size = 22),
      plot.margin = margin(5, 10, 5, 10),
      legend.position = "none"
    )
}

p_reads <- ggplot(sample_dt, aes(sample_rank, raw_reads, fill = study_color)) +
  geom_col(width = 1, na.rm = TRUE) +
  geom_hline(yintercept = raw_reads_threshold, linetype = "dashed", color = "#D62728", linewidth = 1.1) +
  geom_vline(
    data = study_boundaries[boundary_rank < max(sample_dt$sample_rank)],
    aes(xintercept = boundary_rank),
    inherit.aes = FALSE,
    color = "grey35",
    linewidth = 0.15,
    alpha = 0.35
  ) +
  scale_fill_identity() +
  scale_y_log10(labels = comma) +
  labs(
    title = "Sample raw reads by study/sample order",
    y = "Reads (log10)"
  ) +
  theme_sample_qc()

p_chimera <- ggplot(sample_dt, aes(sample_rank, chim_perc_num, fill = study_color)) +
  geom_col(width = 1, na.rm = TRUE) +
  geom_hline(yintercept = chimera_threshold, linetype = "dashed", color = "#D62728", linewidth = 1.1) +
  geom_vline(
    data = study_boundaries[boundary_rank < max(sample_dt$sample_rank)],
    aes(xintercept = boundary_rank),
    inherit.aes = FALSE,
    color = "grey35",
    linewidth = 0.15,
    alpha = 0.35
  ) +
  scale_fill_identity() +
  scale_y_continuous(labels = percent_format(scale = 1, accuracy = 1), limits = c(0, 100)) +
  labs(
    title = "Sample chimera percentage (red line: 25%)",
    y = "Chimera %"
  ) +
  theme_sample_qc()

p_asv <- ggplot(sample_dt, aes(sample_rank, asv_count_num, fill = study_color)) +
  geom_col(width = 1, na.rm = TRUE) +
  geom_hline(yintercept = asv_threshold, linetype = "dashed", color = "#D62728", linewidth = 1.1) +
  geom_vline(
    data = study_boundaries[boundary_rank < max(sample_dt$sample_rank)],
    aes(xintercept = boundary_rank),
    inherit.aes = FALSE,
    color = "grey35",
    linewidth = 0.15,
    alpha = 0.35
  ) +
  scale_fill_identity() +
  scale_y_log10(labels = comma, limits = c(1, NA)) +
  labs(
    title = "Sample ASV count (red line: 30)",
    y = "ASV count"
  ) +
  theme_sample_qc()

combined_sample <- aplot::plot_list(
  gglist = list(p_reads, p_chimera, p_asv),
  ncol = 1,
  heights = c(1, 1, 1),
  output = "patchwork"
)

sample_pdf <- file.path(paths$outdir, paste0(paths$prefix, "_sample_level.pdf"))
sample_png <- file.path(paths$outdir, paste0(paths$prefix, "_sample_level.png"))
ggsave(sample_pdf, combined_sample, width = 28, height = 10.5, limitsize = FALSE)
ggsave(sample_png, combined_sample, width = 28, height = 10.5, dpi = 300, limitsize = FALSE)

message("Done.")
message("Study metrics: ", file.path(paths$outdir, paste0(paths$prefix, "_study_metrics.csv")))
message("Study plot PDF: ", study_pdf)
message("Study plot PNG: ", study_png)
message("Sample plot PDF: ", sample_pdf)
message("Sample plot PNG: ", sample_png)
message("Threshold summary: ", file.path(paths$outdir, paste0(paths$prefix, "_sample_threshold_summary.csv")))

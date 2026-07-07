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
  input = get_arg("input", "QC/sample_filtering_input_table_after_study_filtering.csv"),
  outdir = get_arg("outdir", "QC"),
  prefix = get_arg("prefix", "sample_level_qc_panels_after_study_filtering")
)

raw_reads_threshold <- as.numeric(get_arg("raw-reads-threshold", "1000"))
dechimera_reads_threshold <- as.numeric(get_arg("dechimera-reads-threshold", "1000"))
chimera_threshold <- as.numeric(get_arg("chimera-threshold", "25"))
asv_threshold <- as.numeric(get_arg("asv-threshold", "10"))

dir.create(paths$outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading sample-level QC input...")
dt <- fread(paths$input)
required_cols <- c(
  "study_id", "normalized_sample_id", "raw_reads", "dechimera_reads",
  "chim_perc", "observed_asv"
)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing columns in input table: ", paste(missing_cols, collapse = ", "))
}

dt[, study_id := as.character(study_id)]
dt[, normalized_sample_id := as.character(normalized_sample_id)]
dt[, raw_reads := as.numeric(raw_reads)]
dt[, dechimera_reads := as.numeric(dechimera_reads)]
dt[, chim_perc := as.numeric(chim_perc)]
dt[, observed_asv := as.numeric(observed_asv)]

study_order <- dt[, .(
  study_sample_count = .N,
  median_raw_reads = median(raw_reads, na.rm = TRUE)
), by = study_id][order(-study_sample_count, study_id)]
study_order[, study_rank := .I]

dt <- merge(dt, study_order[, .(study_id, study_sample_count, study_rank)], by = "study_id", all.x = TRUE)
setorder(dt, study_rank, -raw_reads, normalized_sample_id)
dt[, sample_rank := .I]

palette_cycle <- c(
  "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#B07AA1",
  "#EDC948", "#9C755F", "#BAB0AC", "#1F77B4", "#FF7F0E", "#2CA02C"
)
study_order[, study_color := palette_cycle[((study_rank - 1) %% length(palette_cycle)) + 1]]
dt <- merge(dt, study_order[, .(study_id, study_color)], by = "study_id", all.x = TRUE)
setorder(dt, sample_rank)

study_boundaries <- dt[, .(
  start_rank = min(sample_rank),
  end_rank = max(sample_rank),
  center_rank = mean(range(sample_rank)),
  study_sample_count = .N
), by = .(study_id, study_rank)][order(study_rank)]
study_boundaries[, boundary_rank := end_rank + 0.5]

fwrite(
  study_order,
  file.path(paths$outdir, paste0(paths$prefix, "_study_order.csv"))
)
fwrite(
  dt[, .(
    study_id, normalized_sample_id, sample_rank, raw_reads, dechimera_reads,
    chim_perc, observed_asv, study_sample_count
  )],
  file.path(paths$outdir, paste0(paths$prefix, "_sample_ordered_values.csv"))
)

threshold_summary <- data.table(
  metric = c("raw_reads", "dechimera_reads", "chim_perc", "observed_asv"),
  filter_rule = c(
    paste0("raw_reads < ", raw_reads_threshold),
    paste0("dechimera_reads < ", dechimera_reads_threshold),
    paste0("chim_perc > ", chimera_threshold),
    paste0("observed_asv <= ", asv_threshold)
  ),
  threshold = c(raw_reads_threshold, dechimera_reads_threshold, chimera_threshold, asv_threshold),
  failed_samples = c(
    sum(dt$raw_reads < raw_reads_threshold, na.rm = TRUE),
    sum(dt$dechimera_reads < dechimera_reads_threshold, na.rm = TRUE),
    sum(dt$chim_perc > chimera_threshold, na.rm = TRUE),
    sum(dt$observed_asv <= asv_threshold, na.rm = TRUE)
  ),
  non_missing_samples = c(
    sum(!is.na(dt$raw_reads)),
    sum(!is.na(dt$dechimera_reads)),
    sum(!is.na(dt$chim_perc)),
    sum(!is.na(dt$observed_asv))
  )
)
threshold_summary[, failed_ratio := failed_samples / non_missing_samples]
fwrite(
  threshold_summary,
  file.path(paths$outdir, paste0(paths$prefix, "_threshold_summary.csv"))
)

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

bar_width <- 1
threshold_line_width <- 1.1
boundary_line_width <- 0.15

p_reads <- ggplot(dt, aes(sample_rank, raw_reads, fill = study_color)) +
  geom_col(width = bar_width, na.rm = TRUE) +
  geom_line(aes(y = dechimera_reads), color = "#2F2F2F", linewidth = 0.22, na.rm = TRUE) +
  geom_hline(
    yintercept = raw_reads_threshold,
    linetype = "dashed",
    color = "#D62728",
    linewidth = threshold_line_width
  ) +
  geom_hline(
    yintercept = dechimera_reads_threshold,
    linetype = "dotted",
    color = "#D62728",
    linewidth = threshold_line_width
  ) +
  geom_vline(
    data = study_boundaries[boundary_rank < max(dt$sample_rank)],
    aes(xintercept = boundary_rank),
    inherit.aes = FALSE,
    color = "grey35",
    linewidth = boundary_line_width,
    alpha = 0.35
  ) +
  scale_fill_identity() +
  scale_y_log10(labels = comma, limits = c(1, NA), oob = scales::squish) +
  labs(
    title = paste0(
      "Sample raw reads by study/sample order (bars: raw reads; black curve: reads after removing chimeras; red lines: ",
      raw_reads_threshold, " reads)"
    ),
    y = "Reads (log10)"
  ) +
  theme_sample_qc()

p_chimera <- ggplot(dt, aes(sample_rank, chim_perc, fill = study_color)) +
  geom_col(width = bar_width, na.rm = TRUE) +
  geom_hline(
    yintercept = chimera_threshold,
    linetype = "dashed",
    color = "#D62728",
    linewidth = threshold_line_width
  ) +
  geom_vline(
    data = study_boundaries[boundary_rank < max(dt$sample_rank)],
    aes(xintercept = boundary_rank),
    inherit.aes = FALSE,
    color = "grey35",
    linewidth = boundary_line_width,
    alpha = 0.35
  ) +
  scale_fill_identity() +
  scale_y_continuous(labels = percent_format(scale = 1, accuracy = 1), limits = c(0, 100)) +
  labs(
    title = paste0("Sample chimera percentage (red line: ", chimera_threshold, "%)"),
    y = "Chimera %"
  ) +
  theme_sample_qc()

p_asv <- ggplot(dt, aes(sample_rank, observed_asv, fill = study_color)) +
  geom_col(width = bar_width, na.rm = TRUE) +
  geom_hline(
    yintercept = asv_threshold,
    linetype = "dashed",
    color = "#D62728",
    linewidth = threshold_line_width
  ) +
  geom_vline(
    data = study_boundaries[boundary_rank < max(dt$sample_rank)],
    aes(xintercept = boundary_rank),
    inherit.aes = FALSE,
    color = "grey35",
    linewidth = boundary_line_width,
    alpha = 0.35
  ) +
  scale_fill_identity() +
  scale_y_log10(labels = comma, limits = c(1, NA), oob = scales::squish) +
  labs(
    title = paste0("Sample observed ASV count (red line: ASV = ", asv_threshold, ")"),
    y = "Observed ASV"
  ) +
  theme_sample_qc()

combined <- aplot::plot_list(
  gglist = list(p_reads, p_chimera, p_asv),
  ncol = 1,
  heights = c(1, 1, 1),
  output = "patchwork"
)

pdf_file <- file.path(paths$outdir, paste0(paths$prefix, ".pdf"))
png_file <- file.path(paths$outdir, paste0(paths$prefix, ".png"))

ggplot2::ggsave(pdf_file, combined, width = 28, height = 10.5, limitsize = FALSE)
ggplot2::ggsave(png_file, combined, width = 28, height = 10.5, dpi = 300, limitsize = FALSE)

message("Done.")
message("Samples: ", nrow(dt))
message("Studies: ", nrow(study_order))
message("Study order: ", file.path(paths$outdir, paste0(paths$prefix, "_study_order.csv")))
message("Ordered sample values: ", file.path(paths$outdir, paste0(paths$prefix, "_sample_ordered_values.csv")))
message("Threshold summary: ", file.path(paths$outdir, paste0(paths$prefix, "_threshold_summary.csv")))
message("Plot PDF: ", pdf_file)
message("Plot PNG: ", png_file)

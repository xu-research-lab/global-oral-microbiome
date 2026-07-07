# ============================================================
# analysis-oral-GAM-growth-curve-age-0610.R
#
# Adjusted oral microbiome diversity growth chart
#   - Same three-quantile presentation as original (5% / 50% / 95%)
#     fitted with qgam, but with mgcv-style covariate adjustment:
#       site, sequencing region, log10(reads), bioproject_id (random).
#   - Final curves are population-averaged across oral sites
#     (weighted by actual sample share) to give a single oral curve.
#
# Adjusted quantile formula (used at qu = 0.05, 0.50, 0.95):
#   Shannon ~ s(age_cleaned, k = 10, bs = "cs") +
#             site_cleaned +
#             sequencing_region_cleaned +
#             s(log10(total_reads), k = 5) +
#             s(bioproject_id, bs = "re")
#
# Literature support:
#   Pedersen et al. 2019 PeerJ        bs="re" usage + HGAMs
#   Wood 2017 GAM book 2nd ed.        mgcv / smooth-RE duality
#   Knight et al. Best practices       nested tests for batch effects
#   Goodrich et al. 2014 Cell          sequence count as covariate
# ============================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(qgam)
  library(ggplot2)
  library(readr)
  library(dplyr)
})

input_path <- "/home/yanglun/YSD/global_oral/global_oral_v2/metadata-v2/metadata_oral_qc_dedup.tsv"
output_dir <- "/home/yanglun/YSD/global_oral/global_oral_v2/scripts/sample-growth-curve-age/results/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_img              <- paste0(output_dir, "gam-age-adjusted-260610.pdf")
output_csv              <- paste0(output_dir, "merged_gam_age_adjusted_ggplot_pre.csv")
output_peak_driver_csv  <- paste0(output_dir, "gam-age-adjusted-260610_peak18_input_rows.csv")
output_peak_summary_csv <- paste0(output_dir, "gam-age-adjusted-260610_peak18_prediction_summary.csv")
output_model_summary    <- paste0(output_dir, "gam-age-adjusted-260610_model_summary.txt")
output_site_weights_csv <- paste0(output_dir, "gam-age-adjusted-260610_site_weights.csv")

# ============================================================
# 1. Read data + filter
# ============================================================
cat("Reading data...\n")
df <- read_tsv(input_path, show_col_types = FALSE)

required_cols <- c("Shannon", "age_cleaned", "site_cleaned",
                   "sequencing_region_cleaned", "total_reads", "bioproject_id")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Input file missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Diagnostic line marker (same convention as original script)
diagnosis_line_col <- "peak18_input_file_line"
if (diagnosis_line_col %in% names(df)) {
  stop(paste0("Input already contains column ", diagnosis_line_col,
              ", rename it before running."))
}
df[[diagnosis_line_col]] <- seq_len(nrow(df)) + 1L

# All model variables must be present; total_reads > 0 (for log10)
n_before <- nrow(df)
df <- df[
  !is.na(df$Shannon) &
    !is.na(df$age_cleaned) &
    !is.na(df$site_cleaned) & df$site_cleaned != "" &
    !is.na(df$sequencing_region_cleaned) & df$sequencing_region_cleaned != "" &
    !is.na(df$total_reads) & df$total_reads > 0 &
    !is.na(df$bioproject_id) & df$bioproject_id != "",
]
n_after <- nrow(df)
cat(sprintf("Before filter: %d samples; after filter: %d samples (dropped %d NA rows)\n",
            n_before, n_after, n_before - n_after))

# Factorise categorical predictors
df$site_cleaned              <- as.factor(df$site_cleaned)
df$sequencing_region_cleaned <- as.factor(df$sequencing_region_cleaned)
df$bioproject_id             <- as.factor(df$bioproject_id)

cat(sprintf("Age range: %.1f - %.1f years\n",
            min(df$age_cleaned), max(df$age_cleaned)))
cat(sprintf("Sites (n=%d): %s\n",
            nlevels(df$site_cleaned),
            paste(levels(df$site_cleaned), collapse = ", ")))
cat(sprintf("Sequencing regions (n=%d): %s\n",
            nlevels(df$sequencing_region_cleaned),
            paste(levels(df$sequencing_region_cleaned), collapse = ", ")))
cat(sprintf("BioProjects: %d\n\n", nlevels(df$bioproject_id)))

# ============================================================
# 2. Adjusted quantile-GAM fits (5%, 50%, 95%)
# ============================================================
adjusted_formula <- Shannon ~ s(age_cleaned, k = 10, bs = "cs") +
  site_cleaned +
  sequencing_region_cleaned +
  s(log10(total_reads), k = 5) +
  s(bioproject_id, bs = "re")

cat("Fitting adjusted qgam at qu=0.50 (this may take a few minutes)...\n")
fit_50 <- qgam(adjusted_formula, data = df, qu = 0.50, ncores = 4)
cat("  qu=0.50 done\n")
cat("Fitting adjusted qgam at qu=0.05...\n")
fit_05 <- qgam(adjusted_formula, data = df, qu = 0.05, ncores = 4)
cat("  qu=0.05 done\n")
cat("Fitting adjusted qgam at qu=0.95...\n")
fit_95 <- qgam(adjusted_formula, data = df, qu = 0.95, ncores = 4)
cat("  qu=0.95 done\n\n")

# Save full summary of the median model (per-term EDF / p / concurvity)
sink(output_model_summary)
cat("=== Adjusted qgam summary (qu = 0.50) ===\n\n")
print(summary(fit_50))
cat("\n\n=== Concurvity check (qu = 0.50) ===\n\n")
print(concurvity(fit_50))
sink()
cat("Median-model summary saved to:", output_model_summary, "\n\n")

# ============================================================
# 3. Prediction grid + population-averaged marginalisation
# ============================================================
age_grid  <- seq(min(df$age_cleaned), max(df$age_cleaned), length.out = 300)
all_sites <- levels(df$site_cleaned)

site_weights <- df %>%
  count(site_cleaned, name = "n") %>%
  mutate(w = n / sum(n))
cat("Site weights (sample-share based):\n")
print(site_weights)
write_csv(site_weights, output_site_weights_csv)

ref_region <- levels(df$sequencing_region_cleaned)[
  which.max(table(df$sequencing_region_cleaned))
]
ref_reads <- median(df$total_reads)
cat(sprintf("\nReference sequencing region: %s\n", ref_region))
cat(sprintf("Reference depth (median): %.0f reads\n\n", ref_reads))

# Grid: age x site, with other covariates fixed at reference levels
pred_grid <- expand.grid(
  age_cleaned  = age_grid,
  site_cleaned = all_sites,
  stringsAsFactors = FALSE
)
pred_grid$sequencing_region_cleaned <- ref_region
pred_grid$total_reads               <- ref_reads
pred_grid$bioproject_id             <- df$bioproject_id[1]  # placeholder; excluded below
pred_grid <- pred_grid %>%
  left_join(site_weights[, c("site_cleaned", "w")], by = "site_cleaned")

# Predict each quantile, then weight-average across sites at each age
predict_marginal <- function(fit) {
  pp <- as.numeric(predict(fit, newdata = pred_grid,
                           exclude = "s(bioproject_id)"))
  pred_grid$pred <- pp
  pred_grid %>%
    group_by(age_cleaned) %>%
    summarise(value = sum(pred * w) / sum(w), .groups = "drop") %>%
    arrange(age_cleaned)
}

cat("Generating marginal predictions for all three quantiles...\n")
pred_50 <- predict_marginal(fit_50)
pred_05 <- predict_marginal(fit_05)
pred_95 <- predict_marginal(fit_95)

pred_df <- data.frame(
  age_cleaned = pred_50$age_cleaned,
  Median      = pred_50$value,
  Lower       = pred_05$value,
  Upper       = pred_95$value
)
write_csv(pred_df, output_csv)
cat("Predicted curves exported to:", output_csv, "\n\n")

# ============================================================
# 4. 18-year-old peak diagnostics (mirrors original script)
# ============================================================
peak_age_center     <- 18
peak_age_half_width <- 2
peak_age_min        <- peak_age_center - peak_age_half_width
peak_age_max        <- peak_age_center + peak_age_half_width

peak_pred_df  <- pred_df[pred_df$age_cleaned >= peak_age_min &
                           pred_df$age_cleaned <= peak_age_max, ]
peak_input_df <- df[df$age_cleaned >= peak_age_min &
                      df$age_cleaned <= peak_age_max, ]

if (nrow(peak_pred_df) > 0) {
  i_med <- which.max(peak_pred_df$Median)
  i_low <- which.max(peak_pred_df$Lower)
  i_up  <- which.max(peak_pred_df$Upper)
  peak_median_age   <- peak_pred_df$age_cleaned[i_med]
  peak_median_value <- peak_pred_df$Median[i_med]
  peak_lower_age    <- peak_pred_df$age_cleaned[i_low]
  peak_lower_value  <- peak_pred_df$Lower[i_low]
  peak_upper_age    <- peak_pred_df$age_cleaned[i_up]
  peak_upper_value  <- peak_pred_df$Upper[i_up]
} else {
  peak_median_age <- peak_median_value <- NA_real_
  peak_lower_age  <- peak_lower_value  <- NA_real_
  peak_upper_age  <- peak_upper_value  <- NA_real_
}

if (nrow(peak_input_df) > 0) {
  peak_shannon_min    <- min(peak_input_df$Shannon,    na.rm = TRUE)
  peak_shannon_mean   <- mean(peak_input_df$Shannon,   na.rm = TRUE)
  peak_shannon_median <- median(peak_input_df$Shannon, na.rm = TRUE)
  peak_shannon_max    <- max(peak_input_df$Shannon,    na.rm = TRUE)
} else {
  peak_shannon_min <- peak_shannon_mean <-
    peak_shannon_median <- peak_shannon_max <- NA_real_
}

peak_summary_df <- data.frame(
  diagnostic_age_center             = peak_age_center,
  diagnostic_age_half_width         = peak_age_half_width,
  diagnostic_age_min                = peak_age_min,
  diagnostic_age_max                = peak_age_max,
  n_input_rows_in_window            = nrow(peak_input_df),
  shannon_min_in_window             = peak_shannon_min,
  shannon_mean_in_window            = peak_shannon_mean,
  shannon_median_in_window          = peak_shannon_median,
  shannon_max_in_window             = peak_shannon_max,
  median_curve_peak_age_in_window   = peak_median_age,
  median_curve_peak_value_in_window = peak_median_value,
  lower_curve_peak_age_in_window    = peak_lower_age,
  lower_curve_peak_value_in_window  = peak_lower_value,
  upper_curve_peak_age_in_window    = peak_upper_age,
  upper_curve_peak_value_in_window  = peak_upper_value
)
write_csv(peak_summary_df, output_peak_summary_csv)
cat("Peak-18 summary exported to:", output_peak_summary_csv, "\n")

# Per-sample predictions inside the 18yo window (using each sample's own
# site/region/depth/project, but excluding the BioProject random offset).
if (nrow(peak_input_df) > 0) {
  pred_input <- peak_input_df[, c("age_cleaned", "site_cleaned",
                                  "sequencing_region_cleaned",
                                  "total_reads", "bioproject_id")]
  peak_input_df$peak18_predicted_median <-
    as.numeric(predict(fit_50, pred_input, exclude = "s(bioproject_id)"))
  peak_input_df$peak18_predicted_lower  <-
    as.numeric(predict(fit_05, pred_input, exclude = "s(bioproject_id)"))
  peak_input_df$peak18_predicted_upper  <-
    as.numeric(predict(fit_95, pred_input, exclude = "s(bioproject_id)"))
  peak_input_df$peak18_residual_to_median <-
    peak_input_df$Shannon - peak_input_df$peak18_predicted_median
  peak_input_df$peak18_rank_by_shannon_in_window <-
    rank(-peak_input_df$Shannon, ties.method = "min")
  peak_input_df$peak18_shannon_percentile_in_window <-
    100 * rank(peak_input_df$Shannon, ties.method = "average") / nrow(peak_input_df)
  global_shannon_ecdf <- ecdf(df$Shannon)
  peak_input_df$peak18_shannon_percentile_global <-
    100 * global_shannon_ecdf(peak_input_df$Shannon)

  peak_input_df <- peak_input_df[order(-peak_input_df$Shannon,
                                       peak_input_df$age_cleaned), ]
  peak_diag_cols <- c(
    diagnosis_line_col,
    "age_cleaned", "Shannon",
    "site_cleaned", "sequencing_region_cleaned",
    "total_reads", "bioproject_id",
    "peak18_rank_by_shannon_in_window",
    "peak18_shannon_percentile_in_window",
    "peak18_shannon_percentile_global",
    "peak18_predicted_median",
    "peak18_predicted_lower",
    "peak18_predicted_upper",
    "peak18_residual_to_median"
  )
  peak_diag_cols <- intersect(peak_diag_cols, names(peak_input_df))
  peak_input_df <- peak_input_df[, c(peak_diag_cols,
                                     setdiff(names(peak_input_df), peak_diag_cols))]
}
write_csv(peak_input_df, output_peak_driver_csv)
cat("Peak-18 driver rows exported to:", output_peak_driver_csv, "\n\n")

# ============================================================
# 5. Plot (three curves, same look as original script)
# ============================================================
cat("Generating plot...\n")
n_samples <- nrow(df)

p <- ggplot() +
  # A. Background scatter
  geom_point(data = df, aes(x = age_cleaned, y = Shannon),
             color = "grey70", alpha = 0.1, size = 0.8) +
  # B. Median trajectory (solid red)
  geom_line(data = pred_df, aes(x = age_cleaned, y = Median),
            color = "#D55E00", linewidth = 1.2) +
  # C. 5th and 95th adjusted-quantile boundaries (dashed blue)
  geom_line(data = pred_df, aes(x = age_cleaned, y = Lower),
            color = "#0072B2", linetype = "dashed", linewidth = 1) +
  geom_line(data = pred_df, aes(x = age_cleaned, y = Upper),
            color = "#0072B2", linetype = "dashed", linewidth = 1) +
  labs(
    title    = paste0("Adjusted Oral Microbiome Diversity Growth Chart (n=", n_samples, ")"),
    subtitle = "Solid Red: Adjusted Median (50%); Dashed Blue: 5%-95% Adjusted Quantile Range  |  Population-averaged across oral sites",
    x = "Host Age (Years)",
    y = "Shannon Diversity Index"
  ) +
  scale_x_continuous(breaks = seq(0, 110, 20)) +
  theme_bw() +
  theme(
    plot.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    axis.title    = element_text(size = 14),
    axis.text     = element_text(size = 12)
  )

ggsave(output_img, plot = p, width = 10, height = 6)
cat("\nFigure saved to:", output_img, "\n")

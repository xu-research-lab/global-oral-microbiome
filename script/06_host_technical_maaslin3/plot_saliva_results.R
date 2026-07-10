#!/usr/bin/env Rscript
# ============================================================================
# MaAsLin3 v2 — Saliva Results Plots
# References original plot_saliva_heatmaps.R for heatmap style
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

setwd("/beegfs/huangxiaochang/global_oral/new_analysis/maaslin3_v2")

Q_THRESHOLD <- 0.1
TAXONOMY <- "/beegfs/huangxiaochang/global_oral/taxonomy.tsv"

OUT_DIR <- "saliva_plots"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

geography_order <- c("North America", "Asia", "Europe", "South America", "Oceania")
age_order <- c("0-6 months", "6-36 months", "3-12 years", "12-18 years", "18-60 years", ">60 years")

# ============================================================================
# Helper functions (from original plot_saliva_heatmaps.R)
# ============================================================================
clean_genus <- function(x) {
  x <- trimws(as.character(x))
  has_genus <- grepl("g__", x)
  extracted <- regmatches(x, regexpr("g__[^;]*", x))
  x[has_genus] <- sub("^g__", "", extracted)
  x <- sub("_.+$", "", x)
  ifelse(is.na(x) | x == "" | x == "g__", NA_character_, x)
}

extract_rank <- function(taxon, prefix) {
  parts <- strsplit(as.character(taxon), ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  hit <- parts[startsWith(parts, prefix)]
  if (length(hit) == 0 || hit[1] == prefix) return(NA_character_)
  hit[1]
}

build_family_map <- function() {
  tax <- read.delim(TAXONOMY, stringsAsFactors = FALSE, check.names = FALSE)
  tax_col <- if ("Taxon" %in% colnames(tax)) "Taxon" else colnames(tax)[2]
  fam <- vapply(tax[[tax_col]], extract_rank, character(1), prefix = "f__")
  gen <- vapply(tax[[tax_col]], extract_rank, character(1), prefix = "g__")
  gen_clean <- clean_genus(gen)
  df <- data.frame(Genus = gen_clean, Family = fam, stringsAsFactors = FALSE) %>%
    filter(!is.na(Genus), !is.na(Family), Family != "f__") %>%
    count(Genus, Family, name = "n") %>%
    group_by(Genus) %>%
    arrange(desc(n), Family, .by_group = TRUE) %>%
    slice(1) %>% ungroup()
  setNames(df$Family, df$Genus)
}

family_map <- build_family_map()

make_axis_labels <- function(genera, family_map) {
  family <- family_map[genera]
  family[is.na(family)] <- "f__"
  genus <- paste0("g__", genera)
  family_width <- max(nchar(family), na.rm = TRUE) + 2
  labels <- paste0(str_pad(family, width = family_width, side = "right"), genus)
  names(labels) <- genera
  labels
}

order_features <- function(data) {
  data %>%
    mutate(sig_signed_logq = ifelse(qval_individual < Q_THRESHOLD, signed_logq, 0)) %>%
    group_by(Genus) %>%
    summarise(
      mean_signed_logq = mean(sig_signed_logq, na.rm = TRUE),
      min_q = min(qval_individual, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_signed_logq, min_q, Genus) %>%
    pull(Genus)
}

# ============================================================================
# Plot 1: Age heatmap (single model, no facets)
# ============================================================================
age_res_file <- "age_saliva_results/all_results.tsv"
if (file.exists(age_res_file)) {
  cat("=== Plot: Saliva Age-associated microbiota ===\n")

  age_all <- read.delim(age_res_file, stringsAsFactors = FALSE)
  age_all$Genus <- clean_genus(age_all$feature)
  age_all <- age_all %>%
    filter(metadata == "age_cleaned_group", !is.na(qval_individual), !is.na(Genus))

  # Prepare: ALL results (sig + non-sig), with sig_label derived from q
  age_data <- age_all %>%
    group_by(Genus, value) %>%
    arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
    slice(1) %>% ungroup() %>%
    mutate(
      q_plot = ifelse(qval_individual <= 0, 1e-300, qval_individual),
      signed_logq = pmin(-log10(q_plot), 10) * sign(coef),
      sig_label = ifelse(qval_individual < Q_THRESHOLD, ifelse(coef > 0, "+", "-"), ""),
      Direction = ifelse(coef > 0, "+", "-")
    )

  sig_rows <- age_data %>% filter(qval_individual < Q_THRESHOLD)
  if (nrow(sig_rows) > 0) {
    # Keep only age values where at least one genus is significant
    kept_values <- unique(sig_rows$value)
    kept_genera <- unique(sig_rows$Genus)

    # Filter data to only significant genera, but include ALL age values (even non-sig)
    plot_data <- age_data %>%
      filter(Genus %in% kept_genera, value %in% kept_values)

    genus_order <- order_features(plot_data)
    value_levels <- age_order[age_order %in% kept_values]

    # Build complete grid (sig + non-sig cells)
    grid <- expand.grid(
      Genus = genus_order, value = value_levels, stringsAsFactors = FALSE
    ) %>%
      left_join(plot_data %>% select(Genus, value, coef, qval_individual, signed_logq, sig_label),
                by = c("Genus", "value")) %>%
      mutate(
        signed_logq = ifelse(is.na(signed_logq), 0, signed_logq),
        sig_label   = ifelse(is.na(sig_label), "", sig_label),
        Genus = factor(Genus, levels = rev(genus_order)),
        value = factor(value, levels = value_levels)
      )

    axis_labels <- make_axis_labels(genus_order, family_map)
    write.table(grid, file.path(OUT_DIR, "saliva_age_heatmap_data.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    width  <- max(6.2, length(value_levels) * 0.30)
    height <- max(5.2, length(genus_order) * 0.18)

    p <- ggplot(grid, aes(x = value, y = Genus, fill = signed_logq)) +
      geom_tile(color = "grey85", linewidth = 0.14) +
      geom_text(aes(label = sig_label), size = 2.2, fontface = "bold") +
      scale_y_discrete(labels = axis_labels) +
      scale_x_discrete(drop = TRUE) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "#F7F7F7", high = "#D73027",
        midpoint = 0, limits = c(-10, 10), name = "signed\n-log10(q)"
      ) +
      labs(
        title = "Saliva Age-Associated Microbiota",
        subtitle = paste0("Reference: 18-60 years; q < ", Q_THRESHOLD,
                          "; +/− indicates coefficient direction"),
        x = NULL, y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5, face = "bold"),
        axis.text.y = element_text(size = 8.2, family = "mono", hjust = 0),
        strip.text = element_text(size = 8.2, face = "bold"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        legend.position = "bottom",
        legend.key.width = unit(0.75, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        plot.title = element_text(hjust = 0.5, size = 11.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 8.0),
        plot.margin = margin(6, 6, 6, 6)
      )

    ggsave(file.path(OUT_DIR, "saliva_age_heatmap.pdf"), p,
           width = width, height = height, limitsize = FALSE)
    ggsave(file.path(OUT_DIR, "saliva_age_heatmap.png"), p,
           width = width, height = height, limitsize = FALSE, dpi = 300)
    cat(sprintf("  Age heatmap: %d genera, %d age levels (%.1f x %.1f in)\n",
                length(genus_order), length(value_levels), width, height))
  }
}

# ============================================================================
# Plot 2 & 3: Sex-by-Geography and Region-by-Geography
# ============================================================================
geo_results_file <- "by_geography_saliva/all_results_by_geography.tsv"
if (!file.exists(geo_results_file)) {
  cat("by_geography_saliva/all_results_by_geography.tsv not found.\n")
  quit(save = "no", status = 0)
}

cat("\n=== Loading by-geography results ===\n")
geo_res <- read.delim(geo_results_file, stringsAsFactors = FALSE)
geo_res$Genus <- clean_genus(geo_res$feature)
geo_res <- geo_res %>% filter(!is.na(Genus))

variable_labels <- c(
  "sex_cleaned" = "Sex",
  "sequencing_region_cleaned" = "16S Region"
)
variable_refs <- c(
  "sex_cleaned" = "Reference: F",
  "sequencing_region_cleaned" = "Reference: V4"
)
variable_titles <- c(
  "sex_cleaned" = "Sex-Associated",
  "sequencing_region_cleaned" = "Primer-Region-Associated"
)

plot_by_geo_heatmap <- function(results, variable, consistent = FALSE,
                                  output_suffix = "heatmap",
                                  title_suffix = "") {
  var_label <- variable_labels[variable]
  if (is.na(var_label)) var_label <- variable

  # ALL results (sig + non-sig)
  var_data <- results %>%
    filter(metadata == variable, !is.na(qval_individual)) %>%
    group_by(Genus, continent_group, value) %>%
    arrange(qval_individual, desc(abs(coef)), .by_group = TRUE) %>%
    slice(1) %>% ungroup() %>%
    mutate(
      q_plot = ifelse(qval_individual <= 0, 1e-300, qval_individual),
      signed_logq = pmin(-log10(q_plot), 10) * sign(coef),
      sig_label = ifelse(qval_individual < Q_THRESHOLD, ifelse(coef > 0, "+", "-"), ""),
      Direction = ifelse(coef > 0, "+", "-")
    )

  sig_rows <- var_data %>% filter(qval_individual < Q_THRESHOLD)
  if (nrow(sig_rows) == 0) {
    cat(sprintf("  No significant results for %s\n", variable)); return(invisible(NULL))
  }

  # Consistent mode: keep only genera significant in >=2 geographies with same direction
  if (consistent) {
    kept <- sig_rows %>%
      distinct(value, Genus, continent_group, Direction) %>%
      count(value, Genus, Direction, name = "SignificantGeographies") %>%
      filter(SignificantGeographies >= 2)
    if (nrow(kept) == 0) {
      cat(sprintf("  No consistent genera for %s\n", variable)); return(invisible(NULL))
    }
    var_data <- var_data %>%
      inner_join(kept %>% select(value, Genus, Direction),
                 by = c("value", "Genus", "Direction"))
    # re-derive sig_rows from filtered var_data
    sig_rows <- var_data %>% filter(qval_individual < Q_THRESHOLD)
  }

  # Per facet, keep only geos with at least one sig genus
  kept_pairs <- sig_rows %>% distinct(value, continent_group)
  kept_genera <- sig_rows %>% distinct(value, Genus)
  plot_data <- var_data %>%
    inner_join(kept_pairs, by = c("value", "continent_group")) %>%
    inner_join(kept_genera, by = c("value", "Genus"))

  genus_order <- order_features(plot_data)
  value_levels <- unique(plot_data$value[order(plot_data$value)])

  # Build per-facet grid
  grid_parts <- list()
  for (v in value_levels) {
    fd <- plot_data %>% filter(value == v)
    fg <- geography_order[geography_order %in% unique(fd$continent_group)]
    fgen <- genus_order[genus_order %in% unique(fd$Genus)]
    if (length(fg) == 0 || length(fgen) == 0) next
    grid_parts[[as.character(v)]] <- expand.grid(
      Genus = fgen, continent_group = fg, value = v, stringsAsFactors = FALSE)
  }
  grid <- bind_rows(grid_parts) %>%
    left_join(plot_data %>% select(Genus, continent_group, value, signed_logq, sig_label),
              by = c("Genus", "continent_group", "value")) %>%
    mutate(
      signed_logq = ifelse(is.na(signed_logq), 0, signed_logq),
      sig_label   = ifelse(is.na(sig_label), "", sig_label),
      Genus = factor(Genus, levels = rev(genus_order)),
      continent_group = factor(continent_group, levels = geography_order),
      value = factor(value, levels = value_levels)
    )

  axis_labels <- make_axis_labels(genus_order, family_map)
  safe_name <- tolower(gsub("_", "-", variable))
  write.table(grid, file.path(OUT_DIR, paste0("saliva_", safe_name, "_", output_suffix, "_data.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  max_geos <- grid %>% distinct(value, continent_group) %>% count(value) %>% pull(n) %>% max()
  width  <- max(6.2, length(value_levels) * max_geos * 0.30)
  height <- max(5.2, length(genus_order) * 0.18)

  p <- ggplot(grid, aes(x = continent_group, y = Genus, fill = signed_logq)) +
    geom_tile(color = "grey85", linewidth = 0.14) +
    geom_text(aes(label = sig_label), size = 2.2, fontface = "bold") +
    facet_wrap(~ value, nrow = 1, scales = "free_x") +
    scale_y_discrete(labels = axis_labels) +
    scale_x_discrete(drop = TRUE) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#D73027",
      midpoint = 0, limits = c(-10, 10), name = "signed\n-log10(q)"
    ) +
    labs(
      title = paste0("Saliva ", variable_titles[variable], " Microbiota (NA downsampled)", title_suffix),
      subtitle = paste0(variable_refs[variable], "; q < ", Q_THRESHOLD,
                        "; +/− indicates coefficient direction"),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.5, face = "bold"),
      axis.text.y = element_text(size = 8.2, family = "mono", hjust = 0),
      strip.text = element_text(size = 8.2, face = "bold"),
      panel.grid = element_blank(),
      panel.spacing = unit(0.2, "lines"),
      legend.position = "bottom",
      legend.key.width = unit(0.75, "cm"),
      legend.key.height = unit(0.25, "cm"),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      plot.title = element_text(hjust = 0.5, size = 11.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 8.0),
      plot.margin = margin(6, 6, 6, 6)
    )

  ggsave(file.path(OUT_DIR, paste0("saliva_", safe_name, "_", output_suffix, ".pdf")), p,
         width = width, height = height, limitsize = FALSE)
  ggsave(file.path(OUT_DIR, paste0("saliva_", safe_name, "_", output_suffix, ".png")), p,
         width = width, height = height, limitsize = FALSE, dpi = 300)
  cat(sprintf("  %s %s: %d genera, %d facets (%.1f x %.1f in)\n",
              var_label, output_suffix, length(genus_order), length(value_levels), width, height))
}

cat("\n=== Plot: Sex by Geography ===\n")
plot_by_geo_heatmap(geo_res, "sex_cleaned")

cat("\n=== Plot: Region by Geography (all) ===\n")
plot_by_geo_heatmap(geo_res, "sequencing_region_cleaned")

cat("\n=== Plot: Region by Geography (consistent: >=2 geos, same direction) ===\n")
plot_by_geo_heatmap(geo_res, "sequencing_region_cleaned",
                    consistent = TRUE,
                    output_suffix = "consistent_heatmap",
                    title_suffix = " consistent")

# Summary
sig_by_var_geo <- geo_res %>%
  filter(!is.na(qval_individual), qval_individual < Q_THRESHOLD) %>%
  group_by(metadata, continent_group) %>%
  summarise(n_assoc = n(), n_genus = n_distinct(Genus), .groups = "drop")
write.csv(sig_by_var_geo, file.path(OUT_DIR, "saliva_significant_summary_by_geo.csv"), row.names = FALSE)

cat("\n=== Saliva plots done ===\n")

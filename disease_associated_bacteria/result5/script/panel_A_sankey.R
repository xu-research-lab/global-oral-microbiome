#!/usr/bin/env Rscript
# ============================================================================
# Panel A: cohort provenance Sankey (BioProject-weighted).
# Ribbon width = BioProject count; shade = sample-count rank within group.
# Tuned for Cell Press full width (174 mm).
#
# Three outputs from one builder:
#   panel_A_sankey                 Disease -> Disease category -> Continent -> 16S region
#   panel_A_sankey_no_category     Disease -> Continent -> 16S region
#                                  (ordered and coloured by anatomical site)
#   panel_A_sankey_AnatomicalGroup Disease -> Anatomical site -> Continent -> 16S region
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggforce)
})

source(file.path(Sys.getenv("DA_WORK_DIR", "."), "script", "pub_theme.R"))

DATA_PATH <- file.path(COMP_DIR, "cohort_overview", "cohort_overview_data.csv")
if (!file.exists(DATA_PATH)) stop("Missing: ", DATA_PATH)

CATEGORY_MEMBERS <- list(
  "Oral and dental diseases" = c(
    "Periodontal Disease", "Dental Caries", "Aphthous", "Oral Lichen Planus"
  ),
  "Malignancies" = c(
    "OSCC", "Head and Neck Cancer", "Colorectal Cancer", "Lung Cancer",
    "Breast Cancer", "Esophageal Cancer", "Thyroid Cancer"
  ),
  "Infectious and immune-mediated diseases" = c("COVID-19", "HIV", "RA"),
  "Cardiometabolic and sleep-related disorders" = c(
    "Obesity", "Diabetes", "Hypertension", "Obstructive Sleep Apnea"
  ),
  "Gastrointestinal diseases" = c("IBD", "Gastroesophageal Reflux Disease")
)

CATEGORY_SHORT <- c(
  "Oral and dental diseases" = "Oral / dental",
  "Malignancies" = "Malignancies",
  "Infectious and immune-mediated diseases" = "Infectious / immune",
  "Cardiometabolic and sleep-related disorders" = "Cardiometabolic / sleep",
  "Gastrointestinal diseases" = "GI"
)

CATEGORY_COLOURS <- setNames(
  c("#2C7FB8", "#B2182B", "#7B3294", "#E08214", "#1B7837"),
  names(CATEGORY_MEMBERS)
)

# Anatomical site of disease. List order is the fixed top-to-bottom order.
ANATOMICAL_MEMBERS <- list(
  "Oral cavity" = c(
    "Aphthous", "Dental Caries", "Oral Lichen Planus",
    "Periodontal Disease", "OSCC"
  ),
  "Head/neck & upper airway" = c(
    "Head and Neck Cancer", "Obstructive Sleep Apnea", "Thyroid Cancer"
  ),
  "Gastrointestinal tract" = c(
    "Esophageal Cancer", "Gastroesophageal Reflux Disease",
    "Colorectal Cancer", "IBD"
  ),
  "Respiratory / thoracic" = c("Lung Cancer", "COVID-19", "Breast Cancer"),
  "Systemic" = c("Diabetes", "Hypertension", "Obesity", "HIV", "RA")
)

ANATOMICAL_SHORT <- c(
  "Oral cavity" = "Oral cavity",
  "Head/neck & upper airway" = "Head / neck &\nupper airway",
  "Gastrointestinal tract" = "GI tract",
  "Respiratory / thoracic" = "Respiratory /\nthoracic",
  "Systemic" = "Systemic"
)

# Same five colours as the category version; oral stays blue.
ANATOMICAL_COLOURS <- c(
  "Oral cavity" = "#2C7FB8",
  "Head/neck & upper airway" = "#7B3294",
  "Gastrointestinal tract" = "#1B7837",
  "Respiratory / thoracic" = "#B2182B",
  "Systemic" = "#E08214"
)

DISPLAY_NAME <- c(
  "Aphthous" = "RAS",
  "Colorectal Cancer" = "CRC",
  "Head and Neck Cancer" = "H&N Cancer",
  "Obstructive Sleep Apnea" = "OSA",
  "Gastroesophageal Reflux Disease" = "GERD",
  "Periodontal Disease" = "Periodontal",
  "Esophageal Cancer" = "Esophageal Ca",
  "Breast Cancer" = "Breast Ca",
  "Lung Cancer" = "Lung Ca",
  "Thyroid Cancer" = "Thyroid Ca",
  "Dental Caries" = "Dental Caries",
  "Oral Lichen Planus" = "OLP"
)

display_disease <- function(x) {
  out <- as.character(x)
  hit <- out %in% names(DISPLAY_NAME)
  out[hit] <- unname(DISPLAY_NAME[out[hit]])
  out
}

# Compact geometry for 174 mm width.
SANKEY_SEP <- 0.035
SANKEY_AXIS_W <- 0.10
LABEL_PT_OTHER <- 5.5
LABEL_PT_DISEASE <- 5.5

cohort <- read.csv(DATA_PATH, stringsAsFactors = FALSE, check.names = FALSE)

continent_levels <- cohort %>% count(Continent, sort = TRUE) %>% pull(Continent)
region_levels <- cohort %>% count(Region, sort = TRUE) %>% pull(Region)

# Used only by the compact by-disease figure.
CONTINENT_ABBREV <- c(
  "Asia" = "AS",
  "Europe" = "EU",
  "Africa" = "AF",
  "Oceania" = "OC",
  "North America" = "NA",
  "South America" = "SA"
)

# Case + healthy-control sample counts per (disease, BioProject, sampling site),
# aggregated from the sample metadata so the raw metadata is not needed here.
SITE_COUNTS_PATH <- file.path(COMP_DIR, "cohort_overview", "site_sample_counts.csv")
if (!file.exists(SITE_COUNTS_PATH)) stop("Missing: ", SITE_COUNTS_PATH)

# Sampling-site share of each (disease, BioProject) pair: cases and healthy
# controls of that pair only. Fractions sum to 1 per pair, so a BioProject
# that mixes sites still contributes one BioProject of ribbon width.
load_site_fractions <- function() {
  site_counts <- read.csv(SITE_COUNTS_PATH, stringsAsFactors = FALSE,
                          check.names = FALSE)
  frac <- site_counts %>%
    group_by(Disease, bioproject_id) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()
  check <- frac %>%
    group_by(Disease, bioproject_id) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    inner_join(cohort, by = c("Disease", "bioproject_id"))
  if (nrow(check) != nrow(cohort) || any(check$n != check$n_samples)) {
    stop("Site sample counts do not match cohort_overview_data.csv")
  }
  site_levels <- site_counts %>% count(Site, wt = n, sort = TRUE) %>% pull(Site)
  list(frac = frac, levels = site_levels)
}

shade_ramp <- function(base, n) {
  if (n <= 1) return(base)
  dark <- grDevices::colorRampPalette(c(base, "#1A1A1A"))(10)[3]
  light <- grDevices::colorRampPalette(c(base, "white"))(10)[6]
  grDevices::colorRampPalette(c(dark, light))(n)
}

# Twenty medium-dark colours, ordered so a stride through the list jumps hue.
# Light yellows and greys are avoided: at alpha 0.72 they disappear on white.
DISEASE_PALETTE <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
  "#1199B8", "#F032E6", "#7CB518", "#0E8A7D", "#9A6324",
  "#800000", "#000075", "#E56AAA", "#6C4AB6", "#C4A000",
  "#008B8B", "#E09A4A", "#2E7D4F", "#5C4B8A", "#B85C38"
)

# Every cohort disease must sit in exactly one group.
group_lookup <- function(group_members) {
  long <- stack(group_members)
  dup <- unique(as.character(long$values[duplicated(long$values)]))
  if (length(dup)) stop("Disease in more than one group: ", paste(dup, collapse = ", "))
  map <- setNames(as.character(long$ind), as.character(long$values))
  cohort_dz <- unique(cohort$Disease)
  missing <- setdiff(cohort_dz, names(map))
  extra <- setdiff(names(map), cohort_dz)
  if (length(missing)) stop("Diseases without a group: ", paste(missing, collapse = ", "))
  if (length(extra)) stop("Grouped diseases absent from cohort: ", paste(extra, collapse = ", "))
  map
}

build_sankey <- function(group_members, group_short, group_axis_lab,
                         show_group_axis, group_colours, stem,
                         fixed_order = FALSE, show_group_legend = FALSE,
                         abbreviate_cancers = FALSE,
                         colour_by = c("group_shade", "disease"),
                         compact = FALSE,
                         extra_labels = NULL,
                         tighter = FALSE,
                         show_site = FALSE) {
  colour_by <- match.arg(colour_by)
  label_pt_disease <- if (compact) 7 else LABEL_PT_DISEASE
  # 6 pt sits in the gap above each box; 6.5 pt is taller than that gap.
  label_pt_other <- if (compact) 6 else LABEL_PT_OTHER
  axis_w <- if (compact) 0.08 else SANKEY_AXIS_W
  # Wider gaps on the compact figure so a label fits above each box.
  flow_sep <- if (compact) 0.036 else SANKEY_SEP
  label_fn <- if (abbreviate_cancers) disease_abbrev else display_disease
  if (!is.null(extra_labels)) {
    base_fn <- label_fn
    label_fn <- function(x) {
      out <- base_fn(x)
      hit <- out %in% names(extra_labels)
      out[hit] <- unname(extra_labels[out[hit]])
      out
    }
  }
  group_levels <- names(group_members)
  group_map <- group_lookup(group_members)

  dat <- cohort
  cont_levels_use <- as.character(continent_levels)
  region_levels_use <- as.character(region_levels)
  if (compact) {
    # Top to bottom. V3 is not in the requested list; it sits before V3-V4.
    region_order <- c(
      "V1-V2", "V2-V3", "V1-V3", "V3", "V3-V4",
      "V4", "V4-V5", "V6-V8", "Unknown"
    )
    missing <- setdiff(region_levels_use, region_order)
    if (length(missing)) {
      stop("Region missing from display order: ", paste(missing, collapse = ", "))
    }
    region_levels_use <- intersect(region_order, region_levels_use)
    cont <- as.character(dat$Continent)
    mapped <- CONTINENT_ABBREV[cont]
    dat$Continent <- ifelse(is.na(mapped), cont, unname(mapped))
    mapped_lv <- CONTINENT_ABBREV[cont_levels_use]
    cont_levels_use <- ifelse(is.na(mapped_lv), cont_levels_use, unname(mapped_lv))
  }
  dat$Group <- factor(unname(group_map[dat$Disease]), levels = group_levels)

  disease_stats <- dat %>%
    group_by(Group, Disease) %>%
    summarise(n_bp = n(), n_samples = sum(n_samples), .groups = "drop")

  sankey_group_levels <- if (fixed_order) {
    group_levels
  } else {
    dat %>%
      count(Group, name = "n_bp_pairs") %>%
      mutate(
        Group = as.character(Group),
        n_diseases = vapply(
          Group, function(g) length(group_members[[g]]), integer(1)
        )
      ) %>%
      arrange(desc(n_diseases), desc(n_bp_pairs), Group) %>%
      pull(Group)
  }

  disease_levels <- unlist(lapply(sankey_group_levels, function(g) {
    disease_stats %>%
      filter(as.character(Group) == g) %>%
      mutate(disp = display_disease(Disease)) %>%
      arrange(desc(n_bp), desc(n_samples), disp, Disease) %>%
      pull(Disease)
  }), use.names = FALSE)

  if (colour_by == "group_shade") {
    disease_colour_map <- unlist(lapply(group_levels, function(g) {
      df <- disease_stats %>%
        filter(as.character(Group) == g) %>%
        mutate(disp = display_disease(Disease)) %>%
        arrange(desc(n_samples), disp, Disease)
      stats::setNames(
        shade_ramp(group_colours[[g]], nrow(df)),
        as.character(df$Disease)
      )
    }), use.names = TRUE)
  } else {
    # Stride coprime to the palette length, so diseases that sit next to each
    # other on the axis do not receive neighbouring hues.
    n_dz_col <- length(disease_levels)
    if (n_dz_col > length(DISEASE_PALETTE)) {
      stop("Need ", n_dz_col, " disease colours, palette has ",
           length(DISEASE_PALETTE))
    }
    stride <- 7L
    idx <- ((seq_len(n_dz_col) - 1L) * stride) %% length(DISEASE_PALETTE) + 1L
    disease_colour_map <- stats::setNames(DISEASE_PALETTE[idx], disease_levels)
    if (anyDuplicated(unname(disease_colour_map))) {
      stop("Disease colours are not unique")
    }
  }

  # Compact labels use a narrower space before "(" so the closing
  # parenthesis does not sit against the disease bar.
  label_fmt <- if (compact) {
    "%s<span style=\"font-size:2pt\"> </span>(%s)"
  } else {
    "%s (%s)"
  }
  disease_label_map <- stats::setNames(
    sprintf(
      label_fmt,
      label_fn(as.character(disease_stats$Disease)),
      format(disease_stats$n_bp, big.mark = ",", trim = TRUE)
    ),
    as.character(disease_stats$Disease)
  )

  if (show_site) {
    site <- load_site_fractions()
    agg <- dat %>%
      inner_join(site$frac, by = c("Disease", "bioproject_id")) %>%
      group_by(Disease, Group, Continent, Site, Region) %>%
      summarise(value = sum(frac), .groups = "drop") %>%
      mutate(Site = factor(Site, levels = site$levels))
  } else {
    agg <- dat %>%
      group_by(Disease, Group, Continent, Region) %>%
      summarise(value = n(), .groups = "drop")
  }
  agg <- agg %>%
    mutate(
      DiseaseLabel = factor(
        unname(disease_label_map[as.character(Disease)]),
        levels = unname(disease_label_map[as.character(disease_levels)])
      ),
      GroupLabel = factor(
        unname(group_short[as.character(Group)]),
        levels = unname(group_short[sankey_group_levels])
      ),
      Continent = factor(Continent, levels = cont_levels_use),
      Region = factor(Region, levels = region_levels_use),
      Disease = factor(as.character(Disease), levels = disease_levels)
    )

  if (show_group_axis) {
    axis_cols <- c("DiseaseLabel", "GroupLabel", "Continent", "Region")
    axis_labs <- c("Disease", group_axis_lab, "Continent", "16S region")
  } else if (show_site) {
    axis_cols <- c("DiseaseLabel", "Continent", "Site", "Region")
    axis_labs <- c("Disease", "Continent", "Site", "16S region")
  } else {
    axis_cols <- c("DiseaseLabel", "Continent", "Region")
    axis_labs <- c("Disease", "Continent", "16S region")
  }
  n_axis <- length(axis_cols)
  all_levels <- unlist(lapply(axis_cols, function(col) levels(agg[[col]])),
                       use.names = FALSE)
  if (anyDuplicated(all_levels)) {
    stop("Category shared by two axes: ",
         paste(unique(all_levels[duplicated(all_levels)]), collapse = ", "))
  }
  agg <- agg %>%
    select(all_of(axis_cols), Disease, value) %>%
    as.data.frame()

  gathered <- gather_set_data(agg, seq_len(n_axis))
  if (is.numeric(gathered$x) ||
      all(as.character(gathered$x) %in% as.character(seq_len(n_axis)))) {
    gathered$x <- factor(as.integer(as.character(gathered$x)),
                         levels = seq_len(n_axis), labels = axis_labs)
  } else {
    gathered$x <- factor(as.character(gathered$x),
                         levels = axis_cols, labels = axis_labs)
  }
  gathered$y <- factor(
    as.character(gathered$y),
    levels = unlist(lapply(axis_cols, function(col) levels(agg[[col]])),
                    use.names = FALSE)
  )
  gathered$Disease <- factor(as.character(gathered$Disease), levels = disease_levels)

  gathered_disease <- gathered[gathered$x == "Disease", , drop = FALSE]
  gathered_other <- gathered[gathered$x != "Disease", , drop = FALSE]

  # Padding tuned for 4 axes (1.55 left for disease labels, 0.2 right).
  # Keep the same fraction of the width with fewer axes.
  # Half the bar is 0.04; -0.035 put the closing parenthesis on the bar.
  label_nudge <- if (compact) -0.048 else -0.12
  if (compact) {
    # Left pad is only as wide as the longest disease label. The 18 mm
    # margin used by the full-width figures is not added on top.
    fig_w <- if (show_site) 120 else if (tighter) 94 else 120
    margin_l <- 1.0
    margin_t <- if (tighter) 0.3 else 2.0
    margin_b <- if (tighter) 1.0 else 1.5
    # Half of the "16S region" axis title hangs past the last column.
    margin_r <- if (tighter) 9 else 8
    grDevices::pdf(NULL)
    lab_mm <- max(vapply(unname(disease_label_map), function(s) {
      grob <- gridtext::richtext_grob(
        s,
        gp = grid::gpar(
          fontsize = label_pt_disease, fontface = "bold", fontfamily = FONT
        ),
        padding = grid::unit(0, "pt")
      )
      as.numeric(grid::convertWidth(grid::grobWidth(grob), "mm"))
    }, numeric(1))) + 0.3
    grDevices::dev.off()
    panel_w <- fig_w - margin_l - margin_r
    pad_r <- 0.12
    pad_l <- lab_mm * ((n_axis - 1) + pad_r) / (panel_w - lab_mm) - label_nudge
    x_pad <- c(pad_l, pad_r)
  } else if (n_axis == 4L) {
    x_pad <- c(1.55, 0.2)
  } else {
    span_4 <- 3 + 1.55 + 0.2
    x_pad <- c(1.55, 0.2) * (n_axis - 1) / (span_4 - 1.55 - 0.2)
  }

  p <- ggplot(gathered, aes(x = x, id = id, split = y, value = value)) +
    geom_parallel_sets(
      aes(fill = Disease),
      alpha = 0.72,
      axis.width = axis_w,
      sep = flow_sep,
      n = 50,
      strength = 0.55
    ) +
    geom_parallel_sets_axes(
      axis.width = axis_w,
      sep = flow_sep,
      fill = "grey88",
      colour = "grey55",
      linewidth = 0.2
    ) +
    {
      if (compact) {
        # Names sit in the white gap above each box. Inside the box they
        # land on the ribbons and are hard to read.
        geom_text(
          data = gathered_other,
          stat = "parallel_sets_axes",
          aes(label = after_stat(label), y = after_stat(ymax)),
          axis.width = axis_w,
          sep = flow_sep,
          hjust = 0.5,
          vjust = 0,
          nudge_y = 0.15,
          colour = "grey10",
          size = label_pt_other / ggplot2::.pt,
          lineheight = 0.9,
          fontface = "bold",
          family = FONT
        )
      } else {
        geom_parallel_sets_labels(
          data = gathered_other,
          axis.width = axis_w,
          sep = flow_sep,
          colour = "grey10",
          size = label_pt_other / ggplot2::.pt,
          angle = 0,
          hjust = 0.5,
          lineheight = 0.9,
          fontface = "bold",
          family = FONT
        )
      }
    } +
    {
      if (compact) {
        ggtext::geom_richtext(
          data = gathered_disease,
          stat = "parallel_sets_axes",
          aes(label = after_stat(label)),
          axis.width = axis_w,
          sep = flow_sep,
          hjust = 1,
          nudge_x = label_nudge,
          colour = "grey10",
          size = label_pt_disease / ggplot2::.pt,
          fontface = "bold",
          family = FONT,
          fill = NA,
          label.colour = NA,
          label.size = 0,
          label.padding = grid::unit(0, "pt"),
          label.r = grid::unit(0, "pt")
        )
      } else {
        geom_parallel_sets_labels(
          data = gathered_disease,
          axis.width = axis_w,
          sep = flow_sep,
          colour = "grey10",
          size = label_pt_disease / ggplot2::.pt,
          angle = 0,
          hjust = 1,
          nudge_x = label_nudge,
          fontface = "bold",
          family = FONT
        )
      }
    } +
    scale_x_discrete(expand = expansion(add = x_pad)) +
    {
      if (tighter) {
        # Drop the default headroom above the top bar and below the bottom
        # bar so the axis titles sit close to the plot.
        scale_y_continuous(expand = expansion(mult = c(0.018, 0.002)))
      }
    } +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = disease_colour_map, guide = "none", drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_cell() +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = if (tighter) {
        element_text(
          face = "bold", size = 7, colour = "black",
          margin = margin(t = 0.4, r = 0, b = 0, l = 0)
        )
      } else {
        element_text(
          face = "bold",
          size = if (compact) 7 else BASE_PT,
          colour = "black"
        )
      },
      plot.margin = if (compact) {
        margin(margin_t, margin_r, margin_b, margin_l)
      } else {
        margin(2, 4, 2, 18)
      }
    )

  # Ribbons are filled per disease; this legend shows the group base colour.
  if (show_group_legend) {
    legend_df <- data.frame(
      x = factor("Disease", levels = axis_labs),
      y = 0,
      Group = factor(sankey_group_levels, levels = sankey_group_levels)
    )
    p <- p +
      geom_point(
        data = legend_df, aes(x = x, y = y, colour = Group),
        inherit.aes = FALSE, alpha = 0, size = 0
      ) +
      scale_colour_manual(
        values = group_colours[sankey_group_levels],
        name = group_axis_lab,
        guide = guide_legend(
          nrow = 1, title.position = "left",
          override.aes = list(alpha = 1, size = 3.2, shape = 15)
        )
      ) +
      theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_text(size = 6.5, face = "bold"),
        legend.text = element_text(size = 6.2),
        legend.key.size = unit(3.2, "mm"),
        legend.margin = margin(4, 0, 0, 0)
      )
  }

  n_dz <- n_distinct(dat$Disease)
  if (compact) {
    # Same aspect as the 174 x 190 mm figure, with the empty left gutter removed.
    # The abbrev figure drops the top and bottom headroom and narrows the
    # gaps between Disease, Continent and 16S region.
    width_mm <- if (show_site) 120 else if (tighter) 94 else 120
    height_mm <- if (tighter) 116 else 131
  } else {
    width_mm <- W_FULL_MM
    height_mm <- min(MAX_H_MM, max(160, 20 + n_dz * 8.5))
  }
  save_panel(p, stem, width_mm = width_mm, height_mm = height_mm)
  cat(sprintf(
    "%s: %d BioProjects | %d diseases | %d axes | height %.0f mm\n",
    stem, n_distinct(dat$bioproject_id), n_dz, n_axis, height_mm
  ))
  invisible(p)
}

build_sankey(
  CATEGORY_MEMBERS, CATEGORY_SHORT, "Category",
  show_group_axis = TRUE, group_colours = CATEGORY_COLOURS,
  stem = "panel_A_sankey"
)
build_sankey(
  ANATOMICAL_MEMBERS, ANATOMICAL_SHORT, "Anatomical site",
  show_group_axis = FALSE, group_colours = ANATOMICAL_COLOURS,
  stem = "panel_A_sankey_no_category", fixed_order = TRUE,
  show_group_legend = TRUE, abbreviate_cancers = TRUE
)
build_sankey(
  ANATOMICAL_MEMBERS, ANATOMICAL_SHORT, "Anatomical site",
  show_group_axis = FALSE, group_colours = ANATOMICAL_COLOURS,
  stem = "panel_A_sankey_by_disease", fixed_order = TRUE,
  abbreviate_cancers = TRUE, colour_by = "disease", compact = TRUE
)
build_sankey(
  ANATOMICAL_MEMBERS, ANATOMICAL_SHORT, "Anatomical site",
  show_group_axis = FALSE, group_colours = ANATOMICAL_COLOURS,
  stem = "panel_A_sankey_by_disease_abbrev", fixed_order = TRUE,
  abbreviate_cancers = TRUE, colour_by = "disease", compact = TRUE,
  tighter = TRUE, show_site = TRUE,
  extra_labels = c(
    "Periodontal" = "PD",
    "Dental Caries" = "DC",
    "Hypertension" = "HTN"
  )
)
build_sankey(
  ANATOMICAL_MEMBERS, ANATOMICAL_SHORT, "Anatomical site",
  show_group_axis = TRUE, group_colours = ANATOMICAL_COLOURS,
  stem = "panel_A_sankey_AnatomicalGroup", fixed_order = TRUE
)

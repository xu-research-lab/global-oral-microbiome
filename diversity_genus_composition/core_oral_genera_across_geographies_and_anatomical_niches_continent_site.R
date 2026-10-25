suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

find_script_dir <- function() {
  if (interactive()) return(normalizePath(getwd()))
  args <- commandArgs(trailingOnly = FALSE)
  f <- args[grep("^--file=", args)]
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f))))
  normalizePath(getwd())
}
script_dir <- find_script_dir()

PROJECT_ROOT <- normalizePath(file.path(script_dir, ".."))
DATA_DIR     <- file.path(PROJECT_ROOT, "data", "metadata-and-biom-v2-qc")

BASE_DIR  <- PROJECT_ROOT
PANEL_DIR <- Sys.getenv(
  "FIG3_PANEL_DIR",
  unset = file.path(script_dir, "results", "Fig3")
)
INPUT_CACHE_DIR <- file.path(PANEL_DIR, ".input_cache_0922")

INPUT_BIOM <- file.path(DATA_DIR, "biom_oral_qc_dedup.biom")
INPUT_TAXONOMY <- file.path(DATA_DIR, "taxonomy.tsv")

metadata_path <- Sys.getenv("FIG1_METADATA",
                            unset = file.path(DATA_DIR,
                                              "metadata_oral_qc_dedup.tsv"))
abund_file <- Sys.getenv(
  "FIG1_PHYLUM_TABLE",
  unset = file.path(INPUT_CACHE_DIR, "genus-relative-table_0922.tsv")
)

species_file <- Sys.getenv(
  "FIG1_SPECIES_TABLE",
  unset = file.path(INPUT_CACHE_DIR, "species-relative-table_0922.tsv")
)

PREV_THRESH <- 0.0001
TOP_N       <- 7L
FIG_DPI     <- 600L
BASE_SZ     <- 18
OUT_PREFIX  <- "Fig3"

SITE_DISPLAY <- c("Saliva", "Tooth", "Oral Mucosa", "Oral Swab", "Tongue",
                  "Gingiva")
SITE_LOWER   <- tolower(SITE_DISPLAY)

CONTINENTS <- c("North America", "South America", "Europe", "Asia", "Africa",
                "Oceania")

NAME_TO_ABBR <- c("North America" = "NA", "South America" = "SA",
                  "Europe" = "EU", "Asia" = "AS", "Africa" = "AF",
                  "Oceania" = "OC")

CONT_ABBR <- c("AS", "EU", "AF", "OC", "NA", "SA")
CONTINENT_FACET_ABBR <- c("NA", "EU", "AS", "OC", "AF", "SA")
AGE_ORDER <- c("0-6 months", "6-36 months", "3-12 years", "12-18 years",
               "18-60 years", ">60 years")
AGE_SHORT <- c("0-6 mo", "6-36 mo", "3-12 yr", "12-18 yr", "18-60 yr", ">60 yr")

GENUS_TOP7 <- c("Streptococcus", "Veillonella_A", "Prevotella", "Rothia",
                "Granulicatella", "Haemophilus_D", "Neisseria")
DETECTION_BAR_ORDER <- c("Streptococcus", "Veillonella_A", "Prevotella",
                         "Rothia", "Granulicatella", "Haemophilus_D",
                         "Neisseria")

LOW_FILL  <- "#FFFFFF"
HIGH_FILL <- "#084594"
NA_FILL   <- "#E8E8E8"

PALETTE_AB <- c("#C8102E", "#1F77B4", "#2EA043", "#8E44AD", "#E8732C",
                "#00A2A2", "#C99E10", "#C2185B", "#5A7D9A", "#9C4A1E",
                "#3CB371", "#8B6914", "#0072B2", "#CC79A7", "#009E73",
                "#D55E00")

pick_meta_col <- function(md, preferred, fallback = NULL) {
  if (preferred %in% names(md)) return(preferred)
  if (!is.null(fallback) && fallback %in% names(md)) return(fallback)
  preferred
}

short_taxon_label <- function(x) sub("_[0-9]+$", "", x)

to_title <- function(x) {
  vapply(strsplit(x, " "), function(parts) {
    paste0(toupper(substr(parts, 1, 1)), substr(parts, 2, nchar(parts)),
           collapse = " ")
  }, character(1))
}

fmt_k0 <- function(x) {
  out <- ifelse(
    is.na(x), NA_character_,
    ifelse(x == 0, "0",
           ifelse(x %% 1000 == 0,
                  paste0(formatC(x / 1000, format = "d"), "K"),
                  paste0(formatC(x / 1000, format = "f", digits = 1), "K")))
  )
  unname(out)
}

ensure_prepared_inputs <- function() {
  input_files <- c(INPUT_BIOM, metadata_path, INPUT_TAXONOMY)
  missing_inputs <- input_files[!file.exists(input_files)]
  if (length(missing_inputs) > 0) {
    stop("Figure 3 input not found: ", paste(missing_inputs, collapse = ", "))
  }

  newest_input <- max(file.info(input_files)$mtime, na.rm = TRUE)
  stale <- function(path) {
    !file.exists(path) || file.info(path)$mtime < newest_input
  }
  if (!stale(abund_file) && !stale(species_file)) return(invisible(NULL))

  prep_script <- file.path(script_dir, "figure3_prepare_inputs_0922.py")
  if (!file.exists(prep_script)) {
    stop("input conversion script not found: ", prep_script)
  }

  python_bin <- Sys.getenv("FIG3_PYTHON", unset = "")
  if (!nzchar(python_bin)) {
    found <- Sys.which(c("python3", "python"))
    found <- found[nzchar(found)]
    if (length(found) == 0) stop("Python is required to convert the BIOM input")
    python_bin <- unname(found[1])
  }

  cat("-- rebuilding Figure 3 input cache --\n")
  status <- system2(
    python_bin,
    args = c(
      shQuote(prep_script),
      "--biom", shQuote(INPUT_BIOM),
      "--metadata", shQuote(metadata_path),
      "--taxonomy", shQuote(INPUT_TAXONOMY),
      "--out-dir", shQuote(INPUT_CACHE_DIR)
    )
  )
  if (!identical(status, 0L)) stop("Figure 3 input conversion failed")
  invisible(NULL)
}

prepare_data <- function() {
  ensure_prepared_inputs()
  if (!file.exists(metadata_path)) stop("metadata not found: ", metadata_path)
  if (!file.exists(abund_file))    stop("genus table not found: ", abund_file)

  md <- data.table::fread(metadata_path, header = TRUE, sep = "\t")
  col_site <- pick_meta_col(md, "site_cleaned", "site")
  col_ageg <- pick_meta_col(md, "age_cleaned_group", "age_group")
  col_cont <- pick_meta_col(md, "continent_group", "continent")

  rel <- data.table::fread(abund_file, header = TRUE, sep = "\t", skip = 1)
  data.table::setnames(rel, 1, "Taxon")
  sample_cols <- setdiff(names(rel), "Taxon")
  rel_long <- data.table::melt(rel, id.vars = "Taxon",
                               measure.vars = sample_cols,
                               variable.name = "sample_id",
                               value.name = "Abundance")
  rel_long[, sample_id := as.character(sample_id)]
  rel_long[, Genus := Taxon]
  rel_long <- rel_long[!is.na(Genus) & Genus != ""]
  rel_long[, Short := short_taxon_label(Genus)]
  rel_long[, Present := !is.na(Abundance) & Abundance >= PREV_THRESH]
  rm(rel)

  total_samples <- data.table::uniqueN(rel_long$sample_id)
  rank_dt <- rel_long[Present == TRUE,
                      .(Detected = data.table::uniqueN(sample_id)),
                      by = Short]
  rank_dt[, Prevalence := Detected / total_samples]
  rank_dt <- rank_dt[order(-Prevalence, Short)]
  top9 <- head(rank_dt$Short, TOP_N)
  meta_sub <- md[, .(sample_id,
                     site = get(col_site),
                     cont = get(col_cont),
                     age  = get(col_ageg))]
  meta_sub[, sample_id := as.character(sample_id)]
  meta_sub <- unique(meta_sub, by = "sample_id")

  clean_site_age_cont <- function(d) {
    d[, site := tolower(trimws(as.character(site)))]
    d[, cont := trimws(as.character(cont))]
    d[, age  := gsub("＞", ">", as.character(age))]
    bad_site <- c("other", "others", "other (specify)", "others (specify)")
    d <- d[!is.na(site) & site != "" & !(site %in% bad_site)]
    d
  }
  meta_clean <- clean_site_age_cont(copy(meta_sub))
  present_long <- merge(rel_long[Present == TRUE, .(sample_id, Short)],
                        meta_sub, by = "sample_id", all.x = FALSE)
  present_long <- clean_site_age_cont(present_long)
  present_long[, site_display := factor(to_title(site), levels = SITE_DISPLAY)]
  N_by_site_cont <- meta_clean[!is.na(cont) & cont != "",
                               .(N = .N), by = .(site, cont)]
  N_by_site_age  <- meta_clean[!is.na(age) & age != "",
                               .(N = .N), by = .(site, age)]

  count_present <- function(pd, strata_col) {
    cnt <- pd[, .(N_present = data.table::uniqueN(sample_id)),
              by = c("site", strata_col, "Short")]
    data.table::setnames(cnt, strata_col, "Stratum")
    cnt
  }

  pd_a_obs <- present_long[cont %in% CONTINENTS]
  pd_a_obs[, Contin := NAME_TO_ABBR[cont]]
  cnt_a <- count_present(pd_a_obs, "Contin")
  Nc <- N_by_site_cont[cont %in% CONTINENTS]
  Nc[, Contin := NAME_TO_ABBR[cont]]

  grid_a <- data.table::CJ(site = SITE_LOWER, Stratum = CONT_ABBR, Short = top9)
  pd_a <- merge(grid_a, cnt_a, by = c("site", "Stratum", "Short"), all.x = TRUE)
  pd_a <- merge(pd_a, Nc[, .(site, Stratum = Contin, N)],
                by = c("site", "Stratum"), all.x = TRUE)
  pd_a[, Prevalence := fifelse(N > 0, fifelse(is.na(N_present), 0, N_present) / N,
                              NA_real_)]
  pd_a[, site_display := factor(to_title(site), levels = SITE_DISPLAY)]

  pd_b_obs <- present_long[age %in% AGE_ORDER]
  pd_b_obs[, Ageg := factor(age, levels = AGE_ORDER, labels = AGE_SHORT)]
  cnt_b <- count_present(pd_b_obs, "Ageg")
  Nb <- N_by_site_age[age %in% AGE_ORDER]
  Nb[, Ageg := factor(age, levels = AGE_ORDER, labels = AGE_SHORT)]
  grid_b <- data.table::CJ(site = SITE_LOWER, Stratum = AGE_SHORT, Short = top9)
  pd_b <- merge(grid_b, cnt_b, by = c("site", "Stratum", "Short"), all.x = TRUE)
  pd_b <- merge(pd_b, Nb[, .(site, Stratum = Ageg, N)],
                by = c("site", "Stratum"), all.x = TRUE)
  pd_b[, Prevalence := fifelse(N > 0, fifelse(is.na(N_present), 0, N_present) / N,
                              NA_real_)]
  pd_b[, site_display := factor(to_title(site), levels = SITE_DISPLAY)]

  N_by_continent <- meta_clean[cont %in% CONTINENTS, .(Sample_N = .N),
                               by = .(Continent = cont)]
  continent_count_obs <- present_long[cont %in% CONTINENTS,
                                      .(Count = data.table::uniqueN(sample_id)),
                                      by = .(Continent = cont, Short)]
  continent_count_obs[, Rank := data.table::frank(-Count, ties.method = "first"),
                      by = Continent]
  pd_continent_count <- merge(
    continent_count_obs[, .(Continent, Short, Count, Rank)],
    N_by_continent, by = "Continent", all.x = TRUE
  )
  pd_continent_count[, Prevalence := Count / Sample_N]
  pd_continent_count[, ContLabel := factor(NAME_TO_ABBR[Continent],
                                           levels = CONTINENT_FACET_ABBR)]

  N_by_site <- meta_clean[, .(Sample_N = .N), by = .(Site = site)]
  site_count_obs <- present_long[, .(Count = data.table::uniqueN(sample_id)),
                                 by = .(Site = site, Short)]
  site_count_obs[, Rank := data.table::frank(-Count, ties.method = "first"),
                 by = Site]
  pd_site_count <- merge(
    site_count_obs[, .(Site, Short, Count, Rank)],
    N_by_site, by = "Site", all.x = TRUE
  )
  pd_site_count[, Prevalence := Count / Sample_N]
  pd_site_count[, SiteLabel := factor(to_title(Site), levels = SITE_DISPLAY)]

  prev_a_obs <- pd_a[!is.na(Prevalence)]
  prev_b_obs <- pd_b[!is.na(Prevalence)]

  list(rank_dt = rank_dt, top9 = top9, total_samples = total_samples,
       pd_a = pd_a, pd_b = pd_b,
       pd_continent_count = pd_continent_count,
       pd_site_count = pd_site_count,
       prev_a_obs = prev_a_obs, prev_b_obs = prev_b_obs,
       N_by_site_cont = N_by_site_cont, N_by_site_age = N_by_site_age,
       meta_clean = meta_clean)
}

write_out <- function(plot, stem, width, height, dpi = FIG_DPI) {
  png_path <- file.path(PANEL_DIR, paste0(OUT_PREFIX, "_", stem, "_0922.png"))
  grDevices::png(png_path, width = width, height = height,
                 units = "in", res = dpi, bg = "white")
  print(plot)
  grDevices::dev.off()
  cat("  wrote", basename(png_path), "\n")

  jpeg_path <- file.path(PANEL_DIR, paste0(OUT_PREFIX, "_", stem, "_0922.jpg"))
  grDevices::jpeg(jpeg_path, width = width, height = height,
                  units = "in", res = dpi, quality = 95)
  print(plot)
  grDevices::dev.off()
  cat("  wrote", basename(jpeg_path), "\n")
}

build_detection_count_bars <- function(pd, group_col, group_levels, colour_map,
                                         base = BASE_SZ, ncol = 6L,
                                         bar_width = 0.86,
                                         facet_labeller = ggplot2::label_value,
                                         ref_order = NULL) {  pd <- data.table::copy(pd)
  pd[, Group := factor(as.character(get(group_col)), levels = group_levels)]

  if (!is.null(ref_order)) {
    ref_top9 <- ref_order
  } else {
    ref_level <- as.character(group_levels[1L])
    ref_dt <- pd[as.character(Group) == ref_level]
    data.table::setorder(ref_dt, -Count)
    ref_top9 <- ref_dt$Short[seq_len(min(9L, nrow(ref_dt)))]
    if (length(ref_top9) == 0L) stop("Reference facet has no data")
  }

  pd <- pd[Short %in% ref_top9]
  grid <- data.table::CJ(
    Group = factor(group_levels, levels = group_levels),
    Short = factor(ref_top9, levels = ref_top9)
  )
  pd <- merge(grid, pd, by = c("Group", "Short"), all.x = TRUE)
  pd[is.na(Count), Count := 0L]

  pd[, Ykey := factor(Short, levels = rev(ref_top9))]
  y_lab <- stats::setNames(rev(ref_top9), rev(ref_top9))

  pd[, AxisMax := {
    s <- ceiling(max(Count) / 3000) * 1000
    3L * s
  }, by = Group]
  if (identical(as.character(group_levels[1L]), "NA")) {
    pd[as.character(Group) == "SA",  AxisMax := 1000L]
    pd[as.character(Group) == "AF",  AxisMax := 2000L]
    pd[as.character(Group) == "OC",  AxisMax := 2000L]
  } else if (identical(as.character(group_levels[1L]), "Saliva")) {
    pd[as.character(Group) == "Oral Swab", AxisMax := 5000L]
    pd[as.character(Group) == "Gingiva",   AxisMax := 1000L]
  }

  anchor <- data.table::rbindlist(list(
    pd[, .(Ykey = Ykey[1L], AxisMax = 0L), by = Group],
    pd[, .(Ykey = Ykey[1L], AxisMax = AxisMax[1L]), by = Group]
  ))

  p <- ggplot(pd, aes(y = Ykey)) +
    geom_blank(aes(x = AxisMax, y = Ykey), data = anchor, na.rm = TRUE) +
    geom_col(aes(x = Count, fill = Short), width = bar_width,
             color = "#4A4A4A", linewidth = 0.2) +
    facet_grid(cols = vars(Group), scales = "free_x", space = "fixed",
               labeller = facet_labeller) +
    scale_fill_manual(values = colour_map, guide = "none") +
    scale_x_continuous(labels = fmt_k0,
                       breaks = function(lim) {
                         mx <- lim[2]
                         if (abs(mx - 1000) <  1)  return(c(0, 1000))
                         if (abs(mx - 2000) <  1)  return(c(0, 1000, 2000))
                         if (abs(mx - 5000) <  1)  return(c(0, 2500, 5000))
                         seq(lim[1], lim[2], length.out = 4)
                       },
                       expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(labels = y_lab, limits = rev(ref_top9)) +
    labs(x = "Number of samples", y = NULL) +
    theme_classic(base_size = base) +
    theme(
      panel.background   = element_rect(fill = "#F8F5F0", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "#E0DCD5", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.border       = element_rect(color = "#4A4A4A", fill = NA,
                                        linewidth = 0.5),
      axis.line          = element_blank(),
      axis.text.y        = element_text(size = rel(1.15), face = "bold",
                                        color = "#000000"),
      axis.ticks.y       = element_blank(),
      axis.text.x        = element_text(size = rel(0.95), face = "bold",
                                        color = "#303030"),
      axis.title.x       = element_text(size = rel(1.15), face = "bold",
                                        color = "#303030",
                                        margin = margin(t = 7)),
      strip.text         = element_text(size = rel(1.05), face = "bold",
                                        color = "#303030",
                                        lineheight = 0.9,
                                        margin = margin(b = 5, t = 3)),
      strip.background   = element_blank(),
      panel.spacing      = unit(1.2, "lines"),
      plot.margin        = margin(8, 14, 6, 6),
      plot.background    = element_rect(fill = "white", color = NA)
    )
  p
}

main <- function() {
  if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)
  cat("Figure 3 panel B renderer (genus x site)\n")
  cat("  script dir       :", script_dir, "\n")
  cat("  dedup BIOM       :", INPUT_BIOM, "\n")
  cat("  dedup metadata   :", metadata_path, "\n")
  cat("  taxonomy         :", INPUT_TAXONOMY, "\n")
  cat("  cached genus tab:", abund_file, "\n")
  cat("  PANEL_DIR        :", PANEL_DIR, "\n")
  cat("  detection cutoff : RA >= ", PREV_THRESH, " (",
      PREV_THRESH * 100, "%)\n", sep = "")

  dat <- prepare_data()

  rank_dt <- data.table::copy(dat$rank_dt)
  rank_dt[, Rank := seq_len(.N)]
  cat("\n-- Global detection-count top", TOP_N, "genera (n =",
      dat$total_samples, "samples) --\n")
  print(rank_dt[Rank <= TOP_N, .(Rank, Genus = Short,
                                 GlobalDetectionCount = Detected,
                                 GlobalPrevalence = Prevalence)],
        digits = 4)

  mc <- dat$meta_clean
  cat("\n-- Samples by site --\n")
  print(mc[, .(Samples = .N), by = .(Site = to_title(site))][order(-Samples)])

  data.table::fwrite(
    dat$pd_site_count[Short %in% DETECTION_BAR_ORDER][
      order(match(Site, SITE_LOWER), match(Short, DETECTION_BAR_ORDER)),
      .(Site = to_title(Site), Facet = as.character(SiteLabel),
        Rank, Genus = Short, DetectedSamples = Count,
        Sample_N, DetectionPrevalence = Prevalence)
    ],
    file.path(PANEL_DIR,
              "Fig3_panelB_genus_site_sample_count_stats_0922.csv")
  )

  cat("\n-- Panel B: genus detection counts by site --\n")
  ab_genera <- DETECTION_BAR_ORDER
  col_order <- c(dat$rank_dt$Short[dat$rank_dt$Short %in% ab_genera],
                 setdiff(ab_genera, dat$rank_dt$Short))
  ab_cols <- stats::setNames(PALETTE_AB[seq_along(col_order)], col_order)
  p_b <- build_detection_count_bars(
    dat$pd_site_count, "SiteLabel", SITE_DISPLAY, ab_cols, ncol = 6L,
    ref_order = DETECTION_BAR_ORDER
  )
  write_out(p_b, "panelB_genus_site_sample_count", 14, 3.5, dpi = 900)

  cat("\nPanel B written to:", PANEL_DIR, "\n")
}

if (!identical(Sys.getenv("FIG3_NO_MAIN"), "1")) main()

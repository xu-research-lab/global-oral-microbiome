# Analysis-specific settings and helpers: continent

MIN_SAMPLES_PER_GROUP <- 1L

CONTINENT_ORDER <- c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")

SITE_ORDER <- c("saliva")

valid_groups <- function(meta, min_n = MIN_SAMPLES_PER_GROUP) {
  counts <- table(factor(meta[[COL_CONTINENT]], levels = CONTINENT_ORDER))
  names(counts)[counts >= min_n]
}

valid_age_groups <- function(meta, site = NA, min_n = MIN_SAMPLES_PER_GROUP) {
  valid_groups(meta, min_n)
}

make_one_ova <- function(valid_conts) {
  lapply(valid_conts, function(g) {
    list(site = NA_character_, g1 = "Other", g2 = g, other_groups = setdiff(valid_conts, g), comparison_name = paste0(
      sanitize(g),
      "_enrich"
    ))
  })
}

make_site_age_ova <- function(site, valid_groups) make_one_ova(valid_groups)

make_age_ova_comparisons <- function(meta, sites = NULL, min_n = MIN_SAMPLES_PER_GROUP) {
  vg <- valid_groups(meta, min_n)
  make_one_ova(vg)
}

parse_ova_comparison <- function(comp) {
  for (g in CONTINENT_ORDER) {
    if (identical(paste0(sanitize(g), "_enrich"), comp)) {
      return(list(site = NA_character_, target = g))
    }
  }
  list(site = NA_character_, target = NA_character_)
}

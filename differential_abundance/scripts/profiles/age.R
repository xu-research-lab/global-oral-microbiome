# Analysis-specific settings and helpers: age

MIN_SAMPLES_PER_GROUP <- 10L

SITE_ORDER <- c("saliva")

valid_age_groups <- function(meta, site, min_n = MIN_SAMPLES_PER_GROUP) {
  keep <- meta[[COL_SITE]] == site & meta[[COL_AGE]] %in% AGE_ORDER
  counts <- table(factor(meta[[COL_AGE]][keep], levels = AGE_ORDER))
  names(counts)[counts >= min_n]
}

make_site_age_ova <- function(site, valid_groups) {
  lapply(valid_groups, function(g) {
    list(site = site, g1 = "Other", g2 = g, other_groups = setdiff(valid_groups, g), comparison_name = paste0(
      sanitize(site),
      "_enrich_", sanitize_age_level(g)
    ))
  })
}

make_age_ova_comparisons <- function(meta, sites, min_n = MIN_SAMPLES_PER_GROUP) {
  out <- list()
  for (site in sites) {
    vg <- valid_age_groups(meta, site, min_n)
    out <- c(out, make_site_age_ova(site, vg))
  }
  out
}

parse_ova_comparison <- function(comp) {
  for (site in SITE_ORDER) {
    pre <- paste0(sanitize(site), "_enrich_")
    if (startsWith(comp, pre)) {
      rest <- substr(comp, nchar(pre) + 1L, nchar(comp))
      for (g in AGE_ORDER) {
        if (identical(sanitize_age_level(g), rest)) {
          return(list(site = site, target = g))
        }
      }
    }
  }
  list(site = NA_character_, target = NA_character_)
}

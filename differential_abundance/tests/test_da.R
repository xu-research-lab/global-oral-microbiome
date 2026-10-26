# Run from any directory: Rscript /path/to/differential_abundance/tests/test_da.R
file_arg <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1])
root <- normalizePath(file.path(dirname(file_arg), "../.."))
source(file.path(root, "differential_abundance/scripts/config.R"))
checks <- 0L
check <- function(value) {
  stopifnot(value)
  checks <<- checks + 1L
}

env <- da_profile("site", root)
sys.source(file.path(env$script_dir, "selection.R"), envir = env)
features <- c("abundant", "strong_effect", "tied", "wrong_direction", "no_maaslin", "at_cutoff")
enh <- data.frame(
  feature = features, target = "A", n_methods = 3L,
  linda_lfc = c(1, 1, 1, -1, 1, 1), maaslin2_lfc = c(2, 4, 3, 3, 3, 3),
  ancombc2_lfc = c(1, 1, 1, -1, 1, 1), linda_padj = 0.01,
  maaslin2_padj = c(0.01, 0.01, 0.01, 0.01, 0.2, 0.05), ancombc2_padj = 0.01
)
means <- expand.grid(feature = features, target = c("A", "B"), stringsAsFactors = FALSE)
means$target_mean_abundance <- c(4, 3, 1, 2, 2, 2, rep(1, 6))
selected <- env$select_target_topn(enh, 2L, means, 2L, c("A", "B"))
check(identical(as.character(selected), c("abundant", "strong_effect")))
selected_one <- env$select_target_topn(enh, 2L, means, 1L, c("A", "B"))
check(identical(as.character(selected_one), "strong_effect"))
check(!"tied" %in% attr(selected, "selection")$feature)
single <- env$select_target_topn(enh[1, ], 2L, means[means$feature == features[1], ], 1L, c("A", "B"))
check(identical(as.character(single), "abundant"))
bad_means <- means[-1, ]
check(inherits(try(env$select_target_topn(enh, 2L, bad_means, 2L, c("A", "B")), silent = TRUE), "try-error"))

# Model output keeps the complete contrast and records the configuration once.
tmp <- tempfile("da-tests-")
dir.create(tmp)
all_res <- data.frame(
  feature = c("f1", "f2", "covar"), variable = c("siteX", "siteX", "covarY"),
  padj = c(0.02, 0.01, 0.001), log2FoldChange = c(-1, 2, 8), enriched_in = c("Other", "X", "Y")
)
res <- env$write_comparison_outputs(
  all_res, file.path(tmp, "fit"), "X_enrich", "Other", "X",
  "Site1", "Site2", "^site", data.frame(N = 2L), data.frame(seed = 1L), "site_results.tsv"
)
check(identical(res, list(n_sig = 2L, n_up_g1 = 1L, n_up_g2 = 1L)))
check(setequal(list.files(file.path(tmp, "fit")), c("site_results.tsv", "summary.tsv", "run_config.tsv", "analysis.done")))
check(identical(data.table::fread(file.path(tmp, "fit/site_results.tsv"))$feature, c("f2", "f1")))

# A missing method must not count as a nonsignificant vote or invent a direction.
sparse <- da_profile("sex", root, file.path(tmp, "sparse"))
sparse$.maalind2_rev_map_cache <- c(f1 = "f1", f2 = "f2", f3 = "f3")
for (m in c("linda", "maaslin2")) {
  d <- file.path(sparse$method_dirs[[m]], "saliva_pooled_F_vs_M")
  dir.create(d, recursive = TRUE)
  df <- if (m == "linda") {
    data.frame(feature = c("f1", "f2"), log2FoldChange = c(1, -1), padj = c(.01, .2), enriched_in = c("M", "none"))
  } else {
    data.frame(feature = c("f1", "f3"), log2FoldChange = c(2, -2), padj = c(.02, .01), enriched_in = c("M", "F"))
  }
  data.table::fwrite(df, file.path(d, "sex_results.tsv"), sep = "\t")
  writeLines("done", file.path(d, "analysis.done"))
}
sys.source(file.path(sparse$script_dir, "consensus.R"), envir = sparse)
rob <- data.table::fread(file.path(sparse$out_dir, "consensus_robust.tsv"))
check(identical(rob$feature, "f1") && identical(rob$direction_consensus, "agreed:M"))
full <- data.table::fread(file.path(sparse$out_dir, "full_table.tsv"))
check(all(is.na(full$ancombc2_padj)) && identical(full$n_methods_sig, c(2L, 0L, 1L)))

# Load profiles together to detect accidental sharing of main effects or scopes.
profiles <- lapply(c("age", "site", "continent", "sex"), da_profile, project_root = root)
check(identical(vapply(profiles, function(e) e$scope, ""), c("hc", "all", "all", "hc")))
check(identical(profiles[[1]]$SITE_ORDER, "saliva") && length(profiles[[2]]$SITE_ORDER) == 6L)
check(identical(profiles[[1]]$sanitize_age_level("＞60 years"), "gt60_years"))
unlink(tmp, recursive = TRUE)
cat(checks, "checks passed\n")

#!/usr/bin/env Rscript
# netmoss_get_network.R — call the official NetMoss2::getNetwork on a manifest
# of per-cohort SparCC correlation matrices (one group per invocation).
#
# Usage:
#   Rscript netmoss_get_network.R --manifest <csv> --out <dir>
#
# manifest CSV columns: study,matrix_path,n        (n = cohort sample size)
# Outputs in <dir>:
#   pooled_union.csv   — zero-padded pooled correlation matrix over the union taxa
#   netmoss_run.json   — run metadata (engine, n_studies, n_taxa, n_nonfinite, versions)
# Last stdout line: "NETMOSS_OK <out_dir>" (success marker for the Python caller).
#
# getNetwork is called self-paired (case == control) so the union space is
# exactly the group's own taxa; the two returned matrices must be identical.
# Patch provenance: install_netmoss2.sh (same directory) patches the official
# NetMoss2 source (4 scalar loop headers etc.) and emits LOCAL_PATCH.diff.
# NetMoss2 v0.2.0 directory mode additionally reads the data files (melt/dcast
# via a first column named "X") — hence the "X"-headed zero stubs below — and
# returns 4 elements of which [1]/[2] are the pooled case/control matrices
# (same elements/semantics as the pooled result of NetMoss v1's getNetwork).

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(flag, " needs a value")
  args[i + 1]
}
manifest_path <- get_arg("--manifest")
out_dir       <- get_arg("--out")
if (is.null(manifest_path) || is.null(out_dir))
  stop("usage: Rscript netmoss_get_network.R --manifest <csv> --out <dir>")

suppressPackageStartupMessages(library(NetMoss2))

man <- read.csv(manifest_path, stringsAsFactors = FALSE)
stopifnot(all(c("study", "matrix_path", "n") %in% names(man)),
          !anyDuplicated(man$study),
          all(man$n >= 2),
          all(file.exists(man$matrix_path)))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg))
  sub("^--file=", "", script_arg) else "netmoss_get_network.R"
patch_desc <- file.path(dirname(script_path), "install_netmoss2.sh")

oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)

staging <- file.path(tempdir(), "netmoss_run")
unlink(staging, recursive = TRUE)
dir.create(file.path(staging, "net"), recursive = TRUE)
dir.create(file.path(staging, "data"), recursive = TRUE)

# 17 significant digits so the matrix round-trips exactly
mat17 <- function(m) matrix(format(c(m), digits = 17, trim = TRUE),
                            nrow = nrow(m), ncol = ncol(m),
                            dimnames = dimnames(m))

n_check_failures <- 0L
for (k in seq_len(nrow(man))) {
  m <- as.matrix(read.csv(man$matrix_path[k], row.names = 1,
                          check.names = FALSE))
  storage.mode(m) <- "numeric"
  if (!is.matrix(m) || nrow(m) != ncol(m) ||
      !identical(rownames(m), colnames(m)) ||
      max(abs(diag(m) - 1), na.rm = TRUE) > 1e-12) {
    stop("matrix check failed for ", man$study[k])
  }

  tf <- file.path(staging, "net", paste0(man$study[k], ".txt"))
  write.table(mat17(m), tf, sep = "\t", quote = FALSE,
              row.names = TRUE, col.names = NA)

  # data stub: NetMoss2 derives n from length(colnames) - 1 and also
  # melt/dcast's the data with a first column named "X"
  sf <- file.path(staging, "data", paste0(man$study[k], ".txt"))
  cat(paste(c("X", paste0("s", seq_len(man$n[k]))), collapse = "\t"),
      "\n", file = sf, sep = "")
  cat(paste(c("x", rep("0", man$n[k])), collapse = "\t"), "\n",
      file = sf, sep = "", append = TRUE)

  # read back exactly as getNetwork does; validate
  chk_dat <- tryCatch(
    length(colnames(read.table(sf, header = TRUE, sep = "\t"))) - 1,
    error = function(e) -1)
  chk_ok <- tryCatch({
    z <- read.table(tf, header = TRUE, sep = "\t", row.names = 1)
    isTRUE(nrow(z) == ncol(z) &&
             identical(rownames(z), as.character(rownames(m))) &&
             all.equal(as.numeric(as.matrix(z)), c(m), tolerance = 1e-12))
  }, error = function(e) FALSE)
  ok <- isTRUE(chk_dat == man$n[k]) && isTRUE(chk_ok)
  n_check_failures <- n_check_failures + if (ok) 0L else 1L
  cat(sprintf("  [%d/%d] %s n=%d taxa=%d -> TSV%s\n",
              k, nrow(man), man$study[k], man$n[k], nrow(m),
              if (ok) "" else "  CHECK-FAIL"))
}
if (n_check_failures > 0L)
  stop(n_check_failures, " staged file check failures; aborting")

t0 <- Sys.time()
res <- NetMoss2::getNetwork(
  case_dir        = file.path(staging, "data"),
  control_dir     = file.path(staging, "data"),
  net_case_dir    = file.path(staging, "net"),
  net_control_dir = file.path(staging, "net"))
setwd(oldwd)
P <- res[[1]]
if (!isTRUE(all.equal(P, res[[2]], tolerance = 1e-12)))
  stop("case/control sides of the same input diverged (alignment bug?)")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(mat17(P), file.path(out_dir, "pooled_union.csv"), sep = ",",
            quote = TRUE, row.names = TRUE, col.names = NA)

meta <- list(
  engine = "netmoss2",
  package_version = as.character(packageVersion("NetMoss2")),
  getnetwork_patched = TRUE,
  patched_source = patch_desc,
  n_studies = nrow(man),
  n_taxa = nrow(P),
  n_nonfinite = as.integer(sum(!is.finite(P) & upper.tri(P))),
  n_data_check_failures = as.integer(n_check_failures),
  wall_sec = as.integer(difftime(Sys.time(), t0, units = "secs")))

lines <- "{"          # hand-written JSON, no JSON package dependency
items <- vapply(seq_along(meta), function(i) {
  v <- meta[[i]]
  val <- if (is.logical(v)) tolower(as.character(v))
         else if (is.character(v)) paste0('"', v, '"') else as.character(v)
  paste0('  "', names(meta)[i], '": ', val)
}, character(1))
if (length(items))
  lines <- c(lines, paste0(items,
                           ifelse(seq_along(items) == length(items), "", ",")),
             "}")
writeLines(lines, file.path(out_dir, "netmoss_run.json"))
unlink(staging, recursive = TRUE)

cat("NETMOSS_OK", out_dir, "\n")
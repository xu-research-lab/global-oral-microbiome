# Verify that directory-mode NetMoss pooling includes both cohorts.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
test_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
wrapper <- file.path(test_dir, "..", "scripts", "netmoss_get_network.R")
root <- tempfile("netmoss_test_")
dir.create(root)
on.exit(unlink(root, recursive = TRUE), add = TRUE)
taxa <- c("a", "b", "c", "d")
for (study in c("A", "B")) {
  matrix <- diag(4)
  dimnames(matrix) <- list(taxa, taxa)
  matrix[1, 2] <- matrix[2, 1] <- if (study == "A") 0.5 else -0.2
  write.csv(matrix, file.path(root, paste0(study, ".csv")))
}
manifest <- data.frame(study = c("A", "B"), n = c(30, 50),
                       matrix_path = file.path(root, c("A.csv", "B.csv")))
manifest_path <- file.path(root, "manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
output <- file.path(root, "pooled")
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c(shQuote(wrapper), "--manifest", shQuote(manifest_path),
                    "--out", shQuote(output)))
stopifnot(status == 0L)
actual <- read.csv(file.path(output, "pooled_union.csv"), row.names = 1,
                    check.names = FALSE)
stopifnot(isTRUE(all.equal(actual[1, 2], 0.1017185322805387, tolerance = 1e-12)))
cat("PASS: both cohorts contribute to NetMoss pooling\n")

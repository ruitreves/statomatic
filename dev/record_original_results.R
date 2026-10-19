# Record the results of the original (pre-rebuild) statomatic functions that
# the old-vs-new tests compare against. Run from the package root on a commit
# that still contains the original functions:
#     Rscript dev/record_original_results.R
# Writes tests/testthat/fixtures/original_results.rds (see original() in
# tests/testthat/helper-new.R).

fixture <- "tests/testthat/fixtures/original_results.rds"
if (file.exists(fixture)) file.remove(fixture)

Sys.setenv(STATOMATIC_RECORD_ORIGINAL = "true")
res <- as.data.frame(devtools::test(filter = "new_", reporter = "summary", stop_on_failure = FALSE))
Sys.unsetenv("STATOMATIC_RECORD_ORIGINAL")

saved <- readRDS(fixture)
commit <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) NA)
attr(saved, "recorded_from") <- list(commit = commit, date = Sys.Date(),
                                     r_version = R.version.string)
saveRDS(saved, fixture)

cat("Saved", length(saved), "original results to", fixture, "from commit", commit, "\n")
cat("Tests failed:", sum(res$failed), " errors:", sum(res$error), "\n")

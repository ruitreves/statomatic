# Helpers for the statomatic tests.

# Evaluate expr with a fixed seed, then restore the previous RNG state so the
# existing (unseeded) tests behave exactly as before.
with_seed <- function(seed, expr) {
    had_seed <- exists(".Random.seed", envir = globalenv())
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit({
        if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
        else rm(".Random.seed", envir = globalenv())
    })
    set.seed(seed)
    expr
}

# dummy_data() with a fixed seed.
seeded_dummy_data <- function(seed, ...) {
    d <- with_seed(seed, dummy_data(...))
    list(x = d[[1]], groups = d[[2]]$groups)
}

# Run an old run_* function without printing its progress messages.
quietly <- function(expr) {
    utils::capture.output(res <- expr)
    res
}

# Results of the original (pre-rebuild) statomatic functions, saved once so the
# old-vs-new tests keep working after the original code is removed.
# original("key", expr) returns the saved result for `key`. To (re)record,
# run dev/record_original_results.R on a commit that still has the original
# functions; it sets STATOMATIC_RECORD_ORIGINAL=true, so `expr` (a call to an
# original function) is evaluated and saved. If a test's input data changes,
# its saved results must be recorded again.
original_file <- testthat::test_path("fixtures", "original_results.rds")
original_results <- if (file.exists(original_file)) readRDS(original_file) else list()

original <- function(key, expr) {
    if (key %in% names(original_results)) return(original_results[[key]])
    if (!identical(Sys.getenv("STATOMATIC_RECORD_ORIGINAL"), "true")) {
        stop("No saved result of the original statomatic for '", key, "'. ",
             "Record it with dev/record_original_results.R.", call. = FALSE)
    }
    value <- expr
    original_results[[key]] <<- value
    dir.create(dirname(original_file), showWarnings = FALSE, recursive = TRUE)
    saveRDS(original_results, original_file)
    value
}

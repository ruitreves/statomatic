# Helpers for the new_* rebuild tests.

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

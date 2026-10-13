# new_* multi-group tests must give the same p-values as the original run_* functions.
# Column names differ on purpose: old = Pval_<comparison>, new = pvalue_<comparison>.

d <- seeded_dummy_data(202, group_min = 3, group_max = 6, dim_min = 5, dim_max = 10)

old_pvals <- function(old) {
    old <- old[, colnames(old) != "test", drop = FALSE]
    colnames(old) <- sub("^Pval_", "", colnames(old))
    old
}

new_pvals <- function(new) {
    new <- new[, grepl("^pvalue_", colnames(new)), drop = FALSE]
    colnames(new) <- sub("^pvalue_", "", colnames(new))
    new
}

test_that("new_run_tukey matches run_tukey", {
    old <- old_pvals(quietly(run_tukey(d$x, d$groups)))
    new <- new_pvals(new_run_tukey(d$x, d$groups, padj = FALSE))
    expect_equal(new, old)
})

test_that("new_run_dunnett matches run_dunnett", {
    # Dunnett T3 p-values come from Monte Carlo integration (mvtnorm::pmvt),
    # so repeated runs differ slightly. Use the same seed for both.
    old <- old_pvals(with_seed(1, quietly(run_dunnett(d$x, d$groups))))
    new <- new_pvals(with_seed(1, new_run_dunnett(d$x, d$groups, padj = FALSE)))
    expect_equal(new, old)
})

test_that("padj columns are BH-adjusted pvalue columns", {
    res <- new_run_tukey(d$x, d$groups)
    p <- res[, grepl("^pvalue_", colnames(res))]
    q <- res[, grepl("^padj_", colnames(res))]
    expect_equal(unname(as.matrix(q)), unname(sapply(p, stats::p.adjust, method = "BH")))
})

test_that("sort_columns pairs pvalue and padj per comparison", {
    res <- new_run_tukey(d$x, d$groups, sort_columns = TRUE)
    cols <- setdiff(colnames(res), "test")
    expect_equal(sub("^(pvalue|padj)_", "", cols[c(TRUE, FALSE)]),
                 sub("^(pvalue|padj)_", "", cols[c(FALSE, TRUE)]))
    expect_true(all(grepl("^pvalue_", cols[c(TRUE, FALSE)])))
})

test_that("sort_columns is ignored when padj = FALSE", {
    a <- new_run_tukey(d$x, d$groups, padj = FALSE, sort_columns = TRUE)
    b <- new_run_tukey(d$x, d$groups, padj = FALSE, sort_columns = FALSE)
    expect_identical(a, b)
})

test_that("multi-group tests work on a single gene with padj", {
    one <- d$x[1, , drop = FALSE]
    expect_equal(nrow(new_run_tukey(one, d$groups)), 1)
    expect_equal(nrow(new_run_dunnett(one, d$groups)), 1)
})

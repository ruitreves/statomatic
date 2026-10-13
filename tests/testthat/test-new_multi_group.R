# new_* multi-group tests must give the same p-values as the original run_* functions.
# Column names differ on purpose: old = Pval_<comparison>, new = pvalue_<comparison>.

d <- seeded_dummy_data(202, group_min = 3, group_max = 6, dim_min = 5, dim_max = 10)

# Same data with group names whose sort order depends on the locale
# (e.g. "A_vs_AB" vs "AB_vs_B"), to check that comparisons match by name.
tricky_groups <- c("A", "AB", "B", "C", "CD", "D")[match(d$groups, sort(unique(d$groups)))]

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

# Old and new must have the same comparisons with the same p-values.
# Columns are matched by name, so column order may differ.
expect_same_pvals <- function(new, old) {
    new <- new_pvals(new)
    old <- old_pvals(old)
    expect_setequal(colnames(new), colnames(old))
    expect_equal(new[, colnames(old), drop = FALSE], old)
}

test_that("new_run_tukey matches run_tukey", {
    expect_same_pvals(new_run_tukey(d$x, d$groups, padj = FALSE),
                      quietly(run_tukey(d$x, d$groups)))
    expect_same_pvals(new_run_tukey(d$x, tricky_groups, padj = FALSE),
                      quietly(run_tukey(d$x, tricky_groups)))
})

test_that("new_run_dunnett matches run_dunnett", {
    # Dunnett T3 p-values come from Monte Carlo integration (mvtnorm::pmvt),
    # so repeated runs differ slightly. Use the same seed for both.
    expect_same_pvals(with_seed(1, new_run_dunnett(d$x, d$groups, padj = FALSE)),
                      with_seed(1, quietly(run_dunnett(d$x, d$groups))))
    expect_same_pvals(with_seed(1, new_run_dunnett(d$x, tricky_groups, padj = FALSE)),
                      with_seed(1, quietly(run_dunnett(d$x, tricky_groups))))
})

test_that("new_run_dunn matches run_dunn", {
    expect_same_pvals(new_run_dunn(d$x, d$groups, padj = FALSE),
                      quietly(run_dunn(d$x, d$groups)))
    expect_same_pvals(new_run_dunn(d$x, tricky_groups, padj = FALSE),
                      quietly(run_dunn(d$x, tricky_groups)))
})

test_that("new_run_dunn column order does not depend on locale", {
    old_locale <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", old_locale))
    Sys.setlocale("LC_COLLATE", "C")
    a <- colnames(new_run_dunn(d$x, tricky_groups))
    skip_if(suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8")) == "",
            "en_US.UTF-8 locale not available")
    b <- colnames(new_run_dunn(d$x, tricky_groups))
    expect_identical(a, b)
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
    expect_equal(nrow(new_run_dunn(one, d$groups)), 1)
})

# The two-group tests must give the same p-values as the original run_* functions.
# Column names differ on purpose: old = Pval/Padj, new = pvalue/padj.

d <- seeded_dummy_data(101, group_max = 2, dim_min = 5, dim_max = 10)

test_that("run_ttest matches the original", {
    old <- original("run_ttest", quietly(run_ttest(d$x, d$groups, padj = TRUE)))
    new <- run_ttest(d$x, d$groups)
    expect_equal(new$pvalue, old$Pval)
    expect_equal(new$padj, old$Padj)
    expect_equal(rownames(new), rownames(old))
})

test_that("run_welch matches the original", {
    old <- original("run_welch", quietly(run_welch(d$x, d$groups, padj = TRUE)))
    new <- run_welch(d$x, d$groups)
    expect_equal(new$pvalue, old$Pval)
    expect_equal(new$padj, old$Padj)
    expect_equal(rownames(new), rownames(old))
})

test_that("run_wilcox matches the original", {
    old <- original("run_wilcox", quietly(run_wilcox(d$x, d$groups, padj = TRUE)))
    new <- suppressWarnings(run_wilcox(d$x, d$groups))
    expect_equal(new$pvalue, old$Pval)
    expect_equal(new$padj, old$Padj)
    expect_equal(rownames(new), rownames(old))
})

test_that("two-group output has pvalue, padj, test columns", {
    expect_named(run_ttest(d$x, d$groups), c("pvalue", "padj", "test"))
    expect_named(run_ttest(d$x, d$groups, padj = FALSE), c("pvalue", "test"))
    expect_named(run_ttest(d$x, d$groups, include_test_name = FALSE), c("pvalue", "padj"))
})

test_that("two-group tests work on a single gene", {
    one <- d$x[1, , drop = FALSE]
    res <- run_ttest(one, d$groups)
    expect_equal(nrow(res), 1)
    expect_equal(rownames(res), rownames(one))
    expect_equal(res$padj, res$pvalue)
})

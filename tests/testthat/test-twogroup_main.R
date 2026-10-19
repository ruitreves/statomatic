# Tests for twogroup_main.

# 2 groups of 12. Genes are a mix of normal, skewed and unequal-variance rows
# so that ne, nu and nn are all non-empty.
condition <- rep(c("Ctrl", "Exp"), each = 12)
samples <- paste0("S", 1:24)
x <- with_seed(808, rbind(
    matrix(stats::rnorm(40 * 24, 10, 2), nrow = 40),
    matrix(stats::rexp(20 * 24) + 1, nrow = 20),
    matrix(stats::rnorm(30 * 24, 10, rep(c(1, 3), each = 12)), nrow = 30, byrow = TRUE)
))
dimnames(x) <- list(paste0("gene", 1:90), samples)
# sample IDs in a column (default row names), as read from a CSV
si <- data.frame(sample = samples, condition = condition, sex = rep(c("F", "M"), 12),
                 three = rep(c("A", "B", "C"), each = 8))

run <- function(...) twogroup_main(..., verbose = FALSE)
result_tables <- c("ttest_results", "welch_results", "wilcox_results")
layout <- c("pvalue", "padj", "log2FoldChange_Ctrl_vs_Exp", "Ctrl_mean", "Exp_mean", "test")

test_that("results match the original twogroup_main", {
    old <- original("twogroup_main", quietly(twogroup_main(x, condition, padj = TRUE)))
    new <- run(x, si, "condition")
    for (s in c("ne", "nu", "nn")) {
        expect_gt(nrow(new[[s]]), 0)
        expect_identical(rownames(new[[s]]), rownames(old[[s]]))
    }
    for (k in result_tables) {
        o <- old[[k]]
        n <- new[[k]][rownames(o), ]
        expect_equal(n$pvalue, o$Pval)
        expect_equal(n$padj, o$Padj)
        expect_equal(n$log2FoldChange_Ctrl_vs_Exp, o$log2FoldChange_Ctrl_vs_Exp)
        expect_equal(n[, c("Ctrl_mean", "Exp_mean")], o[, c("Ctrl_mean", "Exp_mean")])
    }
})

test_that("each table has the same layout and names its test", {
    res <- run(x, si, "condition")
    for (k in result_tables) expect_identical(colnames(res[[k]]), layout)
    expect_identical(unique(res$ttest_results$test), "t-test")
    expect_identical(unique(res$welch_results$test), "welch-test")
    expect_identical(unique(res$wilcox_results$test), "wilcox-test")
})

test_that("include_test_name = FALSE and padj = FALSE drop those columns", {
    res <- run(x, si, "condition", include_test_name = FALSE, padj = FALSE)
    for (k in result_tables) {
        expect_identical(colnames(res[[k]]), setdiff(layout, c("padj", "test")))
    }
})

test_that("every gene is in exactly one results table, matching its set", {
    res <- run(x, si, "condition")
    expect_identical(rownames(res$ttest_results), rownames(res$ne))
    expect_identical(rownames(res$welch_results), rownames(res$nu))
    expect_identical(rownames(res$wilcox_results), rownames(res$nn))
    all_genes <- unlist(lapply(res[result_tables], rownames), use.names = FALSE)
    expect_setequal(all_genes, rownames(x))
    expect_false(anyDuplicated(all_genes) > 0)
    expect_identical(rownames(res$foldchanges), rownames(x))
})

test_that("row order of sample_info does not matter", {
    shuffled <- si[with_seed(2, sample(nrow(si))), ]
    expect_identical(run(x, shuffled, "condition"), run(x, si, "condition"))
})

test_that("an empty set gives a 0-row table with the usual columns", {
    no_nu <- run(x, si, "condition", levene_threshold = 0)
    no_nn <- run(x, si, "condition", shapiro_threshold = 0)
    expect_equal(nrow(no_nu$welch_results), 0)
    expect_identical(colnames(no_nu$welch_results), layout)
    expect_equal(nrow(no_nn$wilcox_results), 0)
    expect_identical(colnames(no_nn$wilcox_results), layout)
})

test_that("tied values don't produce warnings", {
    # Before R 4.6.0, wilcox.test(exact = TRUE) can't compute exact p-values with
    # ties and warns for every gene; run_wilcox suppresses that. From R 4.6.0
    # it computes exact conditional p-values and doesn't warn, so this test only
    # exercises the suppression on older R.
    ties <- round(x)
    expect_no_warning(run(ties, si, "condition", shapiro_threshold = 1))   # every gene goes to nn
})

test_that("verbose prints progress, verbose = FALSE is silent", {
    msgs <- capture_messages(twogroup_main(x, si, "condition"))
    expect_true(any(grepl("Running Wilcoxon test", msgs)))
    expect_silent(twogroup_main(x, si, "condition", verbose = FALSE))
})

test_that("more than 2 groups, or more than one factor, stop", {
    expect_error(run(x, si, "three"), "exactly one factor with 2 levels")
    expect_error(run(x, si, c("condition", "sex")), "exactly one factor with 2 levels")
    expect_error(run(x, si, "conditon"), "not found in `sample_info`: conditon")
})

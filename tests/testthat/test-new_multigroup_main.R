# Tests for new_multigroup_main.

# 3 groups x 2 sexes, unequal cell sizes. Genes are a mix of normal, skewed
# and unequal-variance rows so that ne, nu and nn are all non-empty.
group <- rep(c("A", "B", "C"), each = 8)
sex <- c(rep(c("F", "M"), c(5, 3)), rep(c("F", "M"), c(4, 4)), rep(c("F", "M"), c(3, 5)))
samples <- paste0("S", 1:24)
x <- with_seed(707, rbind(
    matrix(stats::rnorm(40 * 24, 10, 2), nrow = 40),
    matrix(stats::rexp(20 * 24) + 1, nrow = 20),
    matrix(stats::rnorm(20 * 24, 10, rep(c(1, 1, 5), each = 8)), nrow = 20, byrow = TRUE)
))
dimnames(x) <- list(paste0("gene", 1:80), samples)
# sample IDs in a column (default row names), as read from a CSV
si <- data.frame(sample = samples, group = group, sex = sex)
two_way <- c("group", "sex")

run <- function(...) new_multigroup_main(..., verbose = FALSE)
result_tables <- c("anova_results", "welch_results", "kw_results")

test_that("one-way results match the original multigroup_main", {
    # Dunnett T3 uses Monte Carlo integration, so both runs use the same seed.
    old <- original("multigroup_main_one_way", with_seed(1, quietly(multigroup_main(x, group, padj = FALSE))))
    new <- with_seed(1, run(x, si, "group", padj = FALSE))
    for (s in c("ne", "nu", "nn")) expect_identical(rownames(new[[s]]), rownames(old[[s]]))
    for (k in result_tables) {
        o <- old[[k]][, colnames(old[[k]]) != "test", drop = FALSE]
        colnames(o) <- sub("^Pval", "pvalue", colnames(o))
        expect_setequal(colnames(o), setdiff(colnames(new[[k]]), "test"))
        expect_equal(new[[k]][rownames(o), colnames(o)], o)
    }
})

test_that("every gene is in exactly one results table, matching its set", {
    res <- run(x, si, two_way)
    expect_gt(nrow(res$ne), 0)
    expect_gt(nrow(res$nu), 0)
    expect_gt(nrow(res$nn), 0)
    expect_identical(rownames(res$anova_results), rownames(res$ne))
    expect_identical(rownames(res$welch_results), rownames(res$nu))
    expect_identical(rownames(res$kw_results), rownames(res$nn))
    all_genes <- unlist(lapply(res[result_tables], rownames), use.names = FALSE)
    expect_setequal(all_genes, rownames(x))
    expect_false(anyDuplicated(all_genes) > 0)
    expect_identical(rownames(res$foldchanges), rownames(x))
    expect_identical(rownames(res$tests), rownames(x))
})

test_that("two-way column order: omnibus, then pvalue/padj/fc per comparison, means, test", {
    res <- run(x, si, two_way)
    comparisons <- c("A_F_vs_A_M", "A_F_vs_B_F", "A_F_vs_B_M", "A_F_vs_C_F", "A_F_vs_C_M",
                     "A_M_vs_B_F", "A_M_vs_B_M", "A_M_vs_C_F", "A_M_vs_C_M", "B_F_vs_B_M",
                     "B_F_vs_C_F", "B_F_vs_C_M", "B_M_vs_C_F", "B_M_vs_C_M", "C_F_vs_C_M")
    per_cmp <- as.vector(rbind(paste0("pvalue_", comparisons), paste0("padj_", comparisons),
                               paste0("log2FoldChange_", comparisons)))
    means <- paste0(c("A_F", "A_M", "B_F", "B_M", "C_F", "C_M"), "_mean")
    anova_cols <- c("pvalue_group", "padj_group", "pvalue_sex", "padj_sex",
                    "pvalue_group:sex", "padj_group:sex")
    expect_identical(colnames(res$anova_results), c(anova_cols, per_cmp, means, "test"))
    expect_identical(colnames(res$welch_results), c("pvalue", "padj", per_cmp, means, "test"))
    expect_identical(colnames(res$kw_results), c("pvalue", "padj", per_cmp, means, "test"))
    expect_identical(unique(res$anova_results$test), "tukey-test")
    expect_identical(unique(res$welch_results$test), "dunnett-test")
    expect_identical(unique(res$kw_results$test), "dunn-test")
})

test_that("results attach the right fold changes and means to each gene", {
    res <- run(x, si, two_way)
    fc_cols <- grep("^log2FoldChange_|_mean$", colnames(res$foldchanges), value = TRUE)
    for (k in result_tables) {
        expect_equal(res[[k]][, fc_cols], res$foldchanges[rownames(res[[k]]), fc_cols])
    }
})

test_that("padj = FALSE leaves out every padj column", {
    res <- run(x, si, two_way, padj = FALSE)
    for (k in result_tables) expect_false(any(grepl("^padj", colnames(res[[k]]))))
})

test_that("row order of sample_info does not matter", {
    shuffled <- si[with_seed(2, sample(nrow(si))), ]
    a <- with_seed(3, run(x, si, two_way))
    b <- with_seed(3, run(x, shuffled, two_way))
    expect_identical(a, b)
})

test_that("an empty set gives a 0-row table with the usual columns", {
    full <- run(x, si, two_way)
    no_nu <- run(x, si, two_way, levene_threshold = 0)     # no Levene p-value is below 0
    no_nn <- run(x, si, two_way, shapiro_threshold = 0)    # no Shapiro p-value is below 0
    expect_equal(nrow(no_nu$welch_results), 0)
    expect_identical(colnames(no_nu$welch_results), colnames(full$welch_results))
    expect_equal(nrow(no_nn$kw_results), 0)
    expect_identical(colnames(no_nn$kw_results), colnames(full$kw_results))
})

test_that("verbose prints progress, verbose = FALSE is silent", {
    msgs <- capture_messages(new_multigroup_main(x, si, two_way))
    expect_true(any(grepl("Sorting genes", msgs)))
    expect_true(any(grepl("ne: 50  nu: 4  nn: 26", msgs)))
    expect_silent(new_multigroup_main(x, si, two_way, verbose = FALSE))
})

test_that("samples are found by row names, a sample column or sample_col", {
    by_column <- with_seed(3, run(x, si, two_way))
    by_rownames <- si[, two_way]
    rownames(by_rownames) <- samples
    expect_identical(with_seed(3, run(x, by_rownames, two_way)), by_column)
    two_id_cols <- cbind(si, id = samples)
    expect_error(run(x, two_id_cols, two_way), "Choose one with `sample_col`")
    expect_identical(with_seed(3, run(x, two_id_cols, two_way, sample_col = "id")), by_column)
})

test_that("extra rows in sample_info are ignored", {
    extra <- rbind(si, data.frame(sample = "S99", group = "D", sex = "F"))
    expect_identical(with_seed(3, run(x, extra, two_way)), with_seed(3, run(x, si, two_way)))
})

test_that("sample_info and factors are checked before any analysis", {
    expect_error(run(x, si, c("group", "sexx")), "not found in `sample_info`: sexx")
    expect_error(run(x, si[-1, ], two_way), "Can't find the samples of x")
    expect_error(run(x, si$group, "group"), "must be a data.frame")
    expect_error(run(x, si, si$group), "not the columns themselves")
    expect_error(run(x, si, 2), "must be the names of columns")
})

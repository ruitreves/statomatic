# Tests for the omnibus tests: run_anova, run_kruskal, and run_welch with
# more than 2 groups. Column names differ on purpose: old = Pval/Padj, new = pvalue/padj.

# 3 groups x 2 sexes with unequal, non-proportional cell sizes, so Type I and
# Type II ANOVA give different main-effect p-values.
group <- rep(c("A", "B", "C"), each = 8)
sex <- c(rep(c("F", "M"), c(5, 3)), rep(c("F", "M"), c(4, 4)), rep(c("F", "M"), c(3, 5)))
batch <- rep(c("b1", "b2"), length.out = 24)
x <- with_seed(303, matrix(stats::rnorm(20 * 24, mean = 10, sd = 2), nrow = 20,
                           dimnames = list(paste0("gene", 1:20), paste0("S", 1:24))))

test_that("run_kruskal matches the original", {
    old <- original("run_kruskal", quietly(run_kruskal(x, group, padj = TRUE)))
    new <- run_kruskal(x, group)
    expect_equal(new$pvalue, unname(old[, "Pval"]))
    expect_equal(new$padj, old$Padj)
    expect_equal(rownames(new), rownames(x))
})

test_that("run_welch matches the original with more than 2 groups", {
    old <- original("run_welch_3_groups", quietly(run_welch(x, group, padj = TRUE)))
    new <- run_welch(x, group)
    expect_equal(new$pvalue, old$Pval)
    expect_equal(new$padj, old$Padj)
})

test_that("one-way run_anova matches the original", {
    # Type I and Type II are identical with a single factor.
    old <- original("run_anova_one_way", quietly(run_anova(x, group, padj = TRUE)))
    new <- run_anova(x, group)
    expect_equal(new$pvalue, old[, 1])
    expect_equal(new$padj, old[, 2])
})

test_that("one-way run_anova returns pvalue/padj columns", {
    expect_named(run_anova(x, group, sort_columns = FALSE), c("pvalue", "padj"))
    expect_named(run_anova(x, group, include_test_name = TRUE, sort_columns = FALSE),
                 c("pvalue", "padj", "test"))
})

test_that("two-way run_anova matches car::Anova type 2", {
    f <- data.frame(group = group, sex = sex)
    new <- run_anova(x, f, padj = FALSE)
    expect_named(new, c("pvalue_group", "pvalue_sex", "pvalue_group:sex"))
    for (i in seq_len(nrow(x))) {
        a <- car::Anova(stats::lm(x[i, ] ~ group * sex), type = 2)
        expect_equal(unlist(new[i, ], use.names = FALSE), a[["Pr(>F)"]][1:3])
    }
})

test_that("two-way run_anova does not depend on factor order", {
    gs <- run_anova(x, data.frame(group = group, sex = sex), padj = FALSE)
    sg <- run_anova(x, data.frame(sex = sex, group = group), padj = FALSE)
    expect_equal(gs$pvalue_group, sg$pvalue_group)
    expect_equal(gs$pvalue_sex, sg$pvalue_sex)
    expect_equal(gs$`pvalue_group:sex`, sg$`pvalue_sex:group`)

    # Sanity check: this design is unbalanced enough that Type I would differ.
    type1_gs <- stats::anova(stats::lm(x[1, ] ~ group * sex))[["Pr(>F)"]][1]
    type1_sg <- stats::anova(stats::lm(x[1, ] ~ sex * group))[["Pr(>F)"]][2]
    expect_false(isTRUE(all.equal(type1_gs, type1_sg)))
})

test_that("three-way run_anova returns all 7 terms", {
    new <- run_anova(x, data.frame(group = group, sex = sex, batch = batch), padj = FALSE)
    expect_equal(ncol(new), 7)
    expect_true(all(c("pvalue_group:sex:batch", "pvalue_sex:batch") %in% colnames(new)))
})

test_that("omnibus tests work on a single gene", {
    one <- x[1, , drop = FALSE]
    expect_equal(nrow(run_kruskal(one, group)), 1)
    expect_equal(nrow(run_anova(one, group)), 1)
    expect_equal(nrow(run_anova(one, data.frame(group = group, sex = sex))), 1)
    expect_equal(nrow(run_welch(one, group)), 1)
})

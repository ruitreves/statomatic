# fold_change must give the same results as the original fold_change on
# data on its original scale. The old function returns a matrix; the new one a data.frame.

group <- rep(c("A", "B", "C"), each = 6)
x <- with_seed(505, matrix(stats::rpois(10 * 18, lambda = rep(c(20, 40, 80), each = 6 * 10)),
                           nrow = 10, byrow = FALSE,
                           dimnames = list(paste0("gene", 1:10), paste0("S", 1:18))))

test_that("fold_change matches the original", {
    old <- original("fold_change", quietly(fold_change(x, group)))
    new <- fold_change(x, group)
    expect_s3_class(new, "data.frame")
    expect_identical(colnames(new), colnames(old))
    expect_equal(as.matrix(new), old)
})

test_that("zero means are handled like fold_change", {
    x0 <- x
    x0["gene1", group == "B"] <- 0                    # mean_B = 0
    x0["gene2", group %in% c("A", "B")] <- 0          # both A and B = 0, C not
    old <- original("fold_change_zero_means", quietly(fold_change(x0, group)))
    new <- fold_change(x0, group)
    expect_equal(as.matrix(new), old)
    expect_equal(new["gene1", "log2FoldChange_A_vs_B"], 1e99)
    expect_equal(new["gene1", "log2FoldChange_B_vs_C"], -1e99)
    expect_equal(new["gene2", "log2FoldChange_A_vs_B"], 0)
})

test_that("log_input = TRUE gives the difference of group means", {
    lx <- log2(x + 1)
    new <- fold_change(lx, group, log_input = TRUE)
    expect_equal(new$log2FoldChange_A_vs_B,
                 unname(rowMeans(lx[, group == "A"]) - rowMeans(lx[, group == "B"])))
    # exactly 2x expression on the original scale is a log2 fold change of 1
    two_x <- log2(matrix(c(4, 4, 2, 2), nrow = 1))
    expect_equal(fold_change(two_x, c("a", "a", "b", "b"), log_input = TRUE)[1, 1], 1)
})

test_that("log_input = TRUE ignores the zero rules", {
    lx <- matrix(c(0, 0, 0, 1, 1, 1), nrow = 1)       # log2(1) = 0 is a real value
    new <- fold_change(lx, c("a", "a", "a", "b", "b", "b"), log_input = TRUE)
    expect_equal(new[1, 1], -1)
})

test_that("works on a single gene and on two groups", {
    one <- fold_change(x[1, , drop = FALSE], group)
    expect_equal(dim(one), c(1, 6))
    two <- fold_change(x[, group != "C"], group[group != "C"])
    expect_named(two, c("log2FoldChange_A_vs_B", "A_mean", "B_mean"))
})

test_that("missing values are ignored in group means", {
    xn <- x
    xn["gene1", 1] <- NA
    new <- fold_change(xn, group)
    expect_false(anyNA(new["gene1", ]))
    expect_equal(new["gene1", "A_mean"], mean(x["gene1", 2:6]))
})

test_that("comparison names use radix order, independent of locale", {
    tricky <- c("A", "AB", "B")[match(group, c("A", "B", "C"))]
    expect_identical(colnames(fold_change(x, tricky))[1:3],
                     c("log2FoldChange_A_vs_AB", "log2FoldChange_A_vs_B", "log2FoldChange_AB_vs_B"))
})

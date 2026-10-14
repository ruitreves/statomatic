# Tests for new_test_norm. Old test_norm returns an unnamed list of 3 subsets;
# the new one returns named ne/nu/nn plus a `tests` table.

group <- rep(c("A", "B", "C"), each = 8)
sex <- c(rep(c("F", "M"), c(5, 3)), rep(c("F", "M"), c(4, 4)), rep(c("F", "M"), c(3, 5)))
x <- with_seed(404, rbind(
    matrix(stats::rnorm(40 * 24, 10, 2), nrow = 40),                                       # mostly ne
    matrix(stats::rexp(20 * 24), nrow = 20),                                               # mostly nn
    matrix(stats::rnorm(20 * 24, 10, rep(c(1, 1, 5), each = 8)), nrow = 20, byrow = TRUE)  # mostly nu
))
dimnames(x) <- list(paste0("gene", 1:80), paste0("S", 1:24))

test_that("one-way new_test_norm matches test_norm", {
    old <- quietly(test_norm(x, group))
    new <- new_test_norm(x, group)
    for (k in 1:3) expect_equal(rownames(new[[k]]), rownames(old[[k]]))
})

test_that("two-way new_test_norm matches test_norm", {
    old <- suppressWarnings(quietly(test_norm(x, group, sex)))
    new <- new_test_norm(x, data.frame(group = group, sex = sex))
    for (k in 1:3) expect_equal(rownames(new[[k]]), rownames(old[[k]]))
})

test_that("new_test_norm returns named subsets and a tests table", {
    new <- new_test_norm(x, group)
    expect_named(new, c("ne", "nu", "nn", "tests"))
    expect_named(new$tests, c("shapiro", "levene", "set"))
    expect_equal(sum(sapply(new[c("ne", "nu", "nn")], nrow)), nrow(x))
    for (s in c("ne", "nu", "nn")) {
        expect_equal(rownames(new[[s]]), rownames(new$tests)[new$tests$set == s])
    }
})

test_that("each threshold controls its own test", {
    base <- new_test_norm(x, group)
    strict_shapiro <- new_test_norm(x, group, shapiro_threshold = 0.2)
    strict_levene <- new_test_norm(x, group, levene_threshold = 0.2)
    expect_gte(nrow(strict_shapiro$nn), nrow(base$nn))
    expect_equal(strict_levene$nn, base$nn)          # Levene does not affect nn
    expect_gte(nrow(strict_levene$nu), nrow(base$nu))
})

test_that("a missing value is dropped for that gene only", {
    x_na <- x
    x_na["gene2", 3] <- NA
    with_na <- new_test_norm(x_na, group)$tests["gene2", ]
    dropped <- new_test_norm(x["gene2", -3, drop = FALSE], group[-3])$tests
    expect_equal(with_na, dropped)
})

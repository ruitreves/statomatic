# Tests for check_factors: matching factors to samples by name, design checks,
# and building the combined group.

samples <- paste0("S", 1:24)
x <- with_seed(606, matrix(stats::rnorm(5 * 24, 10), nrow = 5,
                           dimnames = list(paste0("gene", 1:5), samples)))
f2 <- data.frame(condition = rep(c("Ctrl", "Exp"), each = 12),
                 sex = rep(rep(c("F", "M"), each = 6), 2),
                 row.names = samples)

test_that("two factors: combined group is built and aligned to x", {
    res <- check_factors(x, f2)
    expect_named(res, c("factors", "group"))
    expect_identical(rownames(res$factors), colnames(x))
    expect_identical(as.character(res$group), paste(f2$condition, f2$sex, sep = "_"))
    expect_identical(levels(res$group), c("Ctrl_F", "Ctrl_M", "Exp_F", "Exp_M"))
})

test_that("rows of factors are reordered to match the columns of x", {
    shuffled <- f2[with_seed(1, sample(nrow(f2))), ]
    expect_identical(check_factors(x, shuffled), check_factors(x, f2))
})

test_that("one factor works as a named vector and as a one-column data.frame", {
    three <- rep(c("A", "B", "C"), each = 8)
    vec <- check_factors(x, stats::setNames(three, samples))
    df <- check_factors(x, data.frame(g = three, row.names = samples))
    expect_identical(as.character(vec$group), three)
    expect_identical(vec$group, df$group)
})

test_that("three factors give a three-way combined group", {
    f3 <- cbind(f2, batch = rep(c("b1", "b2"), 12))
    res <- check_factors(x, f3)
    expect_equal(ncol(res$factors), 3)
    expect_equal(nlevels(res$group), 8)
    expect_true("Ctrl_F_b1" %in% levels(res$group))
})

test_that("levels are radix-sorted, so the order does not depend on locale", {
    f <- data.frame(g = rep(c("a", "B", "c"), each = 8), row.names = samples)
    expect_identical(levels(check_factors(x, f)$group), c("B", "a", "c"))
})

test_that("factors that can't be matched to samples by name stop with a fix", {
    unnamed <- f2
    rownames(unnamed) <- NULL
    expect_error(check_factors(x, unnamed), "rownames\\(sample_info\\) <- sample_info\\$sample")
    expect_error(check_factors(x, rep(c("A", "B", "C"), each = 8)), "Can't match")
    no_names <- x
    colnames(no_names) <- NULL
    expect_error(check_factors(no_names, f2), "`x` needs column names")
})

test_that("missing or extra samples stop and are named", {
    expect_error(check_factors(x, f2[-1, ]), "In x but not in factors: S1")
    expect_error(check_factors(x[, -1], f2), "In factors but not in x: S1")
})

test_that("invalid designs stop", {
    expect_error(check_factors(x, f2[, "condition", drop = FALSE]), "Use twogroup_main")
    expect_error(check_factors(x, cbind(f2, a = "p", b = "q")), "at most 3 factors")
    na <- f2
    na$sex[1] <- NA
    expect_error(check_factors(x, na), "missing values in: sex")
    one_level <- cbind(f2, batch = "b1")
    expect_error(check_factors(x, one_level), "Only one level in: batch")
})

test_that("factors that are not fully crossed stop and name the empty cells", {
    no_exp_m <- f2
    no_exp_m$sex[no_exp_m$condition == "Exp"] <- "F"
    expect_error(check_factors(x, no_exp_m), "not fully crossed.*Exp_M")
})

test_that("levels that collide when joined with _ stop", {
    collide <- data.frame(a = rep(c("p_q", "p"), 12), b = rep(c("r", "q_r"), each = 12),
                          row.names = samples)
    expect_error(check_factors(x, collide), "same combined label")
})

test_that("group sizes: fewer than 2 stops, exactly 2 warns", {
    one <- f2
    one$sex[one$condition == "Exp" & one$sex == "M"][-1] <- "F"     # leaves 1 Exp_M sample
    expect_error(check_factors(x, one), "Fewer than 2 in: Exp_M")
    two <- f2
    two$sex[two$condition == "Exp" & two$sex == "M"][-(1:2)] <- "F" # leaves 2 Exp_M samples
    expect_warning(res <- check_factors(x, two), "only 2 samples: Exp_M")
    expect_equal(sum(res$group == "Exp_M"), 2)
})

test_that("numeric columns warn", {
    num <- f2
    num$sex <- ifelse(num$sex == "F", 1, 2)
    expect_warning(check_factors(x, num), "Numeric column\\(s\\) treated as categories: sex")
})

test_that("two_group = TRUE requires exactly one factor with 2 levels", {
    two <- check_factors(x, f2[, "condition", drop = FALSE], two_group = TRUE)
    expect_identical(levels(two$group), c("Ctrl", "Exp"))
    expect_error(check_factors(x, stats::setNames(rep(c("A", "B", "C"), each = 8), samples), two_group = TRUE),
                 "exactly one factor with 2 levels")
    expect_error(check_factors(x, f2, two_group = TRUE), "exactly one factor with 2 levels")
})

test_that("long lists in messages are shortened", {
    f <- data.frame(g = rep(c("A", "B", "C"), each = 8), row.names = paste0("T", 1:24))
    expect_error(check_factors(x, f), "S1, S2, S3, S4, S5 and 19 more")
})

test_that("n_rows uses the singular for one row", {
    expect_identical(n_rows(1), "1 row")
    expect_identical(n_rows(0), "0 rows")
    expect_identical(n_rows(80), "80 rows")
})

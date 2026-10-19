# The multi-group tests must give the same p-values as the original run_* functions.
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

test_that("run_tukey matches the original", {
    expect_same_pvals(run_tukey(d$x, d$groups, padj = FALSE),
                      original("run_tukey", quietly(run_tukey(d$x, d$groups))))
    expect_same_pvals(run_tukey(d$x, tricky_groups, padj = FALSE),
                      original("run_tukey_tricky", quietly(run_tukey(d$x, tricky_groups))))
})

test_that("run_dunnett matches the original", {
    # Dunnett T3 p-values come from Monte Carlo integration (mvtnorm::pmvt),
    # so repeated runs differ slightly. Use the same seed for both.
    expect_same_pvals(with_seed(1, run_dunnett(d$x, d$groups, padj = FALSE)),
                      original("run_dunnett", with_seed(1, quietly(run_dunnett(d$x, d$groups)))))
    expect_same_pvals(with_seed(1, run_dunnett(d$x, tricky_groups, padj = FALSE)),
                      original("run_dunnett_tricky", with_seed(1, quietly(run_dunnett(d$x, tricky_groups)))))
})

test_that("run_dunn matches the original", {
    expect_same_pvals(run_dunn(d$x, d$groups, padj = FALSE),
                      original("run_dunn", quietly(run_dunn(d$x, d$groups))))
    expect_same_pvals(run_dunn(d$x, tricky_groups, padj = FALSE),
                      original("run_dunn_tricky", quietly(run_dunn(d$x, tricky_groups))))
})

test_that("run_dunn column order does not depend on locale", {
    old_locale <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", old_locale))
    Sys.setlocale("LC_COLLATE", "C")
    a <- colnames(run_dunn(d$x, tricky_groups))
    skip_if(suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8")) == "",
            "en_US.UTF-8 locale not available")
    b <- colnames(run_dunn(d$x, tricky_groups))
    expect_identical(a, b)
})

test_that("post-hoc tests and fold changes name comparisons identically in any locale", {
    # Mixed-case names sort differently in C and en_US locales. With levels set
    # once, every function must use the same pair names, direction and order.
    raw <- c("a", "B", "c", "D", "e", "F")[match(d$groups, sort(unique(d$groups)))]
    g <- factor(raw, levels = sort(unique(raw), method = "radix"))
    cmp <- function(res, prefix) sub(paste0("^", prefix), "", grep(paste0("^", prefix), colnames(res), value = TRUE))

    old_locale <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", old_locale))
    for (loc in c("C", "en_US.UTF-8")) {
        if (suppressWarnings(Sys.setlocale("LC_COLLATE", loc)) == "") next
        fc <- cmp(fold_change(d$x, g), "log2FoldChange_")
        expect_identical(cmp(run_tukey(d$x, g, padj = FALSE), "pvalue_"), fc)
        expect_identical(cmp(run_dunnett(d$x, g, padj = FALSE), "pvalue_"), fc)
        expect_identical(cmp(run_dunn(d$x, g, padj = FALSE), "pvalue_"), fc)
    }
})

test_that("padj columns are BH-adjusted pvalue columns", {
    res <- run_tukey(d$x, d$groups)
    p <- res[, grepl("^pvalue_", colnames(res))]
    q <- res[, grepl("^padj_", colnames(res))]
    expect_equal(unname(as.matrix(q)), unname(sapply(p, stats::p.adjust, method = "BH")))
})

test_that("sort_columns pairs pvalue and padj per comparison", {
    res <- run_tukey(d$x, d$groups, sort_columns = TRUE)
    cols <- setdiff(colnames(res), "test")
    expect_equal(sub("^(pvalue|padj)_", "", cols[c(TRUE, FALSE)]),
                 sub("^(pvalue|padj)_", "", cols[c(FALSE, TRUE)]))
    expect_true(all(grepl("^pvalue_", cols[c(TRUE, FALSE)])))
})

test_that("sort_columns is ignored when padj = FALSE", {
    a <- run_tukey(d$x, d$groups, padj = FALSE, sort_columns = TRUE)
    b <- run_tukey(d$x, d$groups, padj = FALSE, sort_columns = FALSE)
    expect_identical(a, b)
})

test_that("multi-group tests work on a single gene with padj", {
    one <- d$x[1, , drop = FALSE]
    expect_equal(nrow(run_tukey(one, d$groups)), 1)
    expect_equal(nrow(run_dunnett(one, d$groups)), 1)
    expect_equal(nrow(run_dunn(one, d$groups)), 1)
})

test_that("dunn_bonferroni matches Adunn.test on random data", {
    # Adunn.test is the modified copy of dunn.test 1.3.6 used by the original
    # statomatic. Results must be identical, including with ties.
    with_seed(42, {
        for (i in 1:200) {
            k <- sample(3:6, 1)
            sizes <- sample(3:12, k, replace = TRUE)
            g <- rep(paste0("g", seq_len(k)), sizes)
            N <- sum(sizes)
            y <- switch(i %% 4 + 1,
                stats::rnorm(N, rep(seq_len(k) / 2, sizes)),        # continuous
                round(stats::rnorm(N, 10, 2)),                      # many ties
                stats::rpois(N, 3),                                 # counts, heavy ties
                stats::rexp(N) * rep(seq_len(k), sizes))            # skewed, unequal spread
            old <- original(paste0("Adunn.test_", i),
                            Adunn.test(y, g, method = "bonferroni", table = FALSE, kw = FALSE, altp = TRUE))
            old_p <- stats::setNames(old$altP.adjusted, sub(" - ", "_vs_", old$comparisons, fixed = TRUE))
            new_p <- dunn_bonferroni(y, g)
            expect_setequal(names(new_p), names(old_p))
            expect_equal(new_p, old_p[names(new_p)])
        }
    })
})

test_that("dunn_bonferroni drops missing values and keeps every pair", {
    y <- with_seed(9, stats::rnorm(18))
    g <- rep(c("A", "B", "C"), each = 6)
    y_na <- y
    y_na[c(2, 9)] <- NA
    expect_equal(dunn_bonferroni(y_na, g), dunn_bonferroni(y[-c(2, 9)], g[-c(2, 9)]))
    expect_named(dunn_bonferroni(y_na, g), c("A_vs_B", "A_vs_C", "B_vs_C"))
})

test_that("dunnett_t3 matches AdunnettT3Test on random data", {
    # AdunnettT3Test is the modified copy of PMCMRplus::dunnettT3Test used by
    # the original statomatic. Its p-values use Monte Carlo integration, so
    # each comparison runs both versions with the same seed.
    seeds <- with_seed(7, sample.int(1e6, 120))
    for (s in seq_along(seeds)) {
        data <- with_seed(seeds[s], {
            k <- sample(3:5, 1)
            sizes <- sample(3:10, k, replace = TRUE)
            N <- sum(sizes)
            g <- rep(paste0("g", seq_len(k)), sizes)
            y <- switch(s %% 4 + 1,
                stats::rnorm(N, rep(seq_len(k) / 2, sizes)),
                stats::rnorm(N, 10, rep(seq_len(k), sizes)),          # unequal variances
                round(stats::rnorm(N, 10, 2)),                         # ties
                { v <- stats::rnorm(N, 5); v[g %in% c("g1", "g2")] <- 3; v })  # two identical constant groups
            list(y = y, g = g)
        })
        old <- original(paste0("AdunnettT3Test_", s),
                        with_seed(s, AdunnettT3Test(data$y, factor(data$g))$p.value))
        cells <- which(!is.na(old), arr.ind = TRUE)
        old_p <- stats::setNames(old[cells], paste0(colnames(old)[cells[, "col"]], "_vs_",
                                                    rownames(old)[cells[, "row"]]))
        new_p <- with_seed(s, dunnett_t3(data$y, data$g))
        expect_setequal(names(new_p), names(old_p))
        expect_equal(new_p, old_p[names(new_p)])
    }
})

test_that("dunnett_t3 gives p = 1 for two identical constant groups", {
    y <- c(rep(3, 5), rep(3, 5), with_seed(3, stats::rnorm(5, 5)))
    p <- with_seed(1, dunnett_t3(y, rep(c("a", "b", "c"), each = 5)))
    expect_equal(p[["a_vs_b"]], 1)
    expect_true(all(p[c("a_vs_c", "b_vs_c")] < 0.05))
})

test_that("dunnett_t3 drops missing values and keeps every pair", {
    y <- with_seed(9, stats::rnorm(18))
    g <- rep(c("A", "B", "C"), each = 6)
    y_na <- y
    y_na[c(2, 9)] <- NA
    expect_equal(with_seed(1, dunnett_t3(y_na, g)), with_seed(1, dunnett_t3(y[-c(2, 9)], g[-c(2, 9)])))
    expect_named(dunnett_t3(y_na, g), c("A_vs_B", "A_vs_C", "B_vs_C"))
})

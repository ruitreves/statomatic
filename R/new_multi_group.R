#' Tukey's HSD test on each row of a data.frame (rebuild of run_tukey)
#' @param x a numeric data.frame or matrix, genes as rows
#' @param var1 group label for each column of x
#' @param padj logical, add Benjamini-Hochberg adjusted p-values
#' @param include_test_name logical, add a `test` column
#' @param sort_columns logical, place each comparison's pvalue and padj columns side by side
#' @return a data.frame of p-values
#' @noRd
new_run_tukey <- function(x, var1, padj = TRUE, include_test_name = TRUE, sort_columns = TRUE) {
    per_gene <- lapply(seq_len(nrow(x)), function(i) {
        model <- stats::aov(unlist(x[i, ]) ~ var1)
        pvalues <- stats::TukeyHSD(model)[[1]][, "p adj"]
        return(pvalues)
    })
    res <- do.call(rbind, per_gene)
    rownames(res) <- rownames(x)
    colnames(res) <- vapply(strsplit(colnames(res), "-", fixed = TRUE), function(pair) {
        paste0(pair[2], "_vs_", pair[1])
    }, character(1))
    res <- make_res(res, padj, sort_columns, include_test_name, test_name = "tukey-test")
    return(res)
}

#' Dunnett's T3 test on each row of a data.frame (rebuild of run_dunnett)
#' @inheritParams new_run_tukey
#' @noRd
new_run_dunnett <- function(x, var1, padj = TRUE, include_test_name = TRUE, sort_columns = TRUE) {
    per_gene <- lapply(seq_len(nrow(x)), function(i) dunnett_t3(unlist(x[i, ]), var1))
    res <- do.call(rbind, per_gene)
    rownames(res) <- rownames(x)
    res <- make_res(res, padj, sort_columns, include_test_name, test_name = "dunnett-test")
    return(res)
}

#' Dunn's test on each row of a data.frame (rebuild of run_dunn)
#' @inheritParams new_run_tukey
#' @noRd
new_run_dunn <- function(x, var1, padj = TRUE, include_test_name = TRUE, sort_columns = TRUE) {
    per_gene <- lapply(seq_len(nrow(x)), function(i) dunn_bonferroni(unlist(x[i, ]), var1))
    res <- do.call(rbind, per_gene)
    rownames(res) <- rownames(x)
    res <- make_res(res, padj, sort_columns, include_test_name, test_name = "dunn-test")
    return(res)
}

#' Dunn's test for one gene: two-sided p-values, Bonferroni-adjusted
#'
#' Gives the same results as Adunn.test(y, g, method = "bonferroni",
#' table = FALSE, kw = FALSE, altp = TRUE), the modified copy of dunn.test
#' 1.3.6 used by the original statomatic:
#'   1. rank all values together (ties get their average rank)
#'   2. z = (mean rank of A - mean rank of B) / sqrt(sigma2 * (1/n_A + 1/n_B)),
#'      sigma2 = N(N+1)/12 - sum(t^3 - t) / (12(N-1)), t = sizes of tied groups
#'   3. p = 2 * (1 - pnorm(|z|)), multiplied by the number of pairs, capped at 1
#' Pairs follow the level order of g, named `<A>_vs_<B>` like the other
#' post-hoc tests. Missing values are dropped.
#'
#' @param y numeric values for one gene
#' @param g group label for each value
#' @return a named numeric vector, one adjusted p-value per pair of groups
#' @noRd
dunn_bonferroni <- function(y, g) {
    g <- factor(g)
    keep <- !is.na(y)
    y <- y[keep]
    g <- g[keep]                                  # keeps all levels, so every gene has the same pairs
    N <- length(y)
    r <- rank(y)
    n <- tabulate(g, nbins = nlevels(g))
    mean_rank <- tapply(r, g, mean)
    ties <- table(y)
    sigma2 <- N * (N + 1) / 12 - sum(ties^3 - ties) / (12 * (N - 1))
    pairs <- utils::combn(nlevels(g), 2)
    z <- (mean_rank[pairs[1, ]] - mean_rank[pairs[2, ]]) /
         sqrt(sigma2 * (1 / n[pairs[1, ]] + 1 / n[pairs[2, ]]))
    p <- pmin(1, 2 * stats::pnorm(-abs(z)) * ncol(pairs))
    stats::setNames(as.numeric(p), paste0(levels(g)[pairs[1, ]], "_vs_", levels(g)[pairs[2, ]]))
}

#' Dunnett's T3 test for one gene
#'
#' Gives the same results as AdunnettT3Test(y, factor(g)), the modified copy
#' of PMCMRplus::dunnettT3Test used by the original statomatic (with the same
#' random seed, see below):
#'   1. Welch t-statistic for each pair: (mean_A - mean_B) / sqrt(s2_A/n_A + s2_B/n_B)
#'   2. Welch-Satterthwaite degrees of freedom, rounded to an integer
#'   3. p = 1 - P(max |T| < |t|) for m independent t variables (the studentized
#'      maximum modulus distribution, m = number of pairs), from mvtnorm::pmvt
#' As in the original's modification, two groups with identical constant
#' values (t and df are NaN) get t = 0 and df = 0, so p = 1 instead of an error.
#' pmvt uses Monte Carlo integration: p-values vary slightly between runs
#' unless the random seed is set. Pairs are computed in the same order as the
#' original, so the same seed gives identical results.
#' Pairs follow the level order of g, named `<A>_vs_<B>`. Missing values are dropped.
#'
#' @param y numeric values for one gene
#' @param g group label for each value
#' @return a named numeric vector, one p-value per pair of groups
#' @noRd
dunnett_t3 <- function(y, g) {
    g <- factor(g)
    keep <- !is.na(y)
    y <- y[keep]
    g <- g[keep]                                  # keeps all levels, so every gene has the same pairs
    n <- tabulate(g, nbins = nlevels(g))
    means <- tapply(y, g, mean)
    vars <- tapply(y, g, stats::var)
    pairs <- utils::combn(nlevels(g), 2)
    i <- pairs[1, ]
    j <- pairs[2, ]
    se2 <- vars[i] / n[i] + vars[j] / n[j]
    t <- (means[i] - means[j]) / sqrt(se2)
    df <- se2^2 / (vars[i]^2 / (n[i]^2 * (n[i] - 1)) + vars[j]^2 / (n[j]^2 * (n[j] - 1)))
    t[is.nan(t)] <- 0                             # identical constant groups: p = 1
    df[is.nan(df)] <- 0                           # df = 0: pmvt uses the normal distribution
    df <- round(df)
    m <- ncol(pairs)
    p <- vapply(seq_len(m), function(k) {
        1 - mvtnorm::pmvt(lower = rep(-abs(t[k]), m), upper = rep(abs(t[k]), m),
                          df = df[k], corr = diag(m))
    }, numeric(1))
    stats::setNames(p, paste0(levels(g)[i], "_vs_", levels(g)[j]))
}

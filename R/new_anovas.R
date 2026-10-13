#' Type II ANOVA on each row of a data.frame (rebuild of run_anova)
#'
#' Fits `y ~ factor1 * factor2 * ...` for each gene and tests every main effect
#' and interaction with a Type II ANOVA (`car::Anova(type = 2)`), so results do
#' not depend on the order of the factors. Any number of factors is supported.
#'
#' @param x a numeric data.frame or matrix, genes as rows
#' @param factors a vector with one group label per column of x (one-way), or a
#'   data.frame with one column per factor (two-way and up). Column names are
#'   used to name the result columns.
#' @param padj logical, add Benjamini-Hochberg adjusted p-values
#' @param include_test_name logical, add a `test` column
#' @param sort_columns logical, place each term's pvalue and padj columns side by side
#' @return a data.frame of p-values. One-way: `pvalue`/`padj`. Two-way and up:
#'   `pvalue_<term>`/`padj_<term>` for each main effect and interaction.
#' @noRd
new_run_anova <- function(x, factors, padj = TRUE, include_test_name = FALSE, sort_columns = TRUE) {
    if (is.null(dim(factors))) factors <- data.frame(var1 = factors)   # a single factor
    factors <- as.data.frame(factors)
    factors[] <- lapply(factors, factor)
    formu <- stats::as.formula(paste("y ~", paste(colnames(factors), collapse = " * ")))

    per_gene <- lapply(seq_len(nrow(x)), function(i) {
        fit <- stats::lm(formu, data = cbind(y = unlist(x[i, ]), factors))
        a <- car::Anova(fit, type = 2)
        stats::setNames(a[["Pr(>F)"]], rownames(a))[rownames(a) != "Residuals"]
    })
    res <- do.call(rbind, per_gene)
    rownames(res) <- rownames(x)
    if (ncol(res) == 1) res <- res[, 1]   # one-way: plain pvalue/padj, like Kruskal-Wallis
    make_res(res, padj, sort_columns, include_test_name, test_name = "anova-test")
}

#' Kruskal-Wallis test on each row of a data.frame (rebuild of run_kruskal)
#' @param x a numeric data.frame or matrix, genes as rows
#' @param var1 group label for each column of x
#' @param padj logical, add Benjamini-Hochberg adjusted p-values
#' @param include_test_name logical, add a `test` column
#' @return a data.frame of p-values with columns `pvalue` and `padj`
#' @noRd
new_run_kw <- function(x, var1, padj = TRUE, include_test_name = FALSE) {
    res <- apply(x, 1, function(x) {
        kw <- stats::kruskal.test(x ~ var1, as.data.frame(x))
        return(kw$p.value)
    })
    res <- make_res(res, padj, sort_columns = FALSE, include_test_name, test_name = "kruskal-wallis-test")
    return(res)
}
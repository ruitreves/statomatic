#' Student's t-test on each row of a data.frame
#' @param x a numeric data.frame or matrix, genes as rows
#' @param var1 group label for each column of x
#' @param padj logical, add Benjamini-Hochberg adjusted p-values
#' @param include_test_name logical, add a `test` column
#' @return a data.frame of p-values
#' @export
run_ttest <- function(x, var1, padj = TRUE, include_test_name = TRUE) {
    res <- apply(x, 1, function(x) {
        test <- stats::t.test(x ~ var1, var.equal = TRUE)
        pvalue <- test$p.value
        return(pvalue)
    })
    res <- make_res(res, padj, sort_columns = FALSE, include_test_name, test_name = "t-test")
    return(res)
}

#' Welch test on each row of a data.frame
#' @inheritParams run_ttest
#' @export
run_welch <- function(x, var1, padj = TRUE, include_test_name = FALSE) {
    res <- apply(x, 1, function(x) {
        welch <- onewaytests::welch.test(x ~ var1, data = as.data.frame(x), verbose = FALSE)
        return(welch$p.value)
    })
    res <- make_res(res, padj, sort_columns = FALSE, include_test_name, test_name = "welch-test")
    return(res)
}

#' Wilcoxon/Mann-Whitney test on each row of a data.frame
#' @inheritParams run_ttest
#' @export
run_wilcox <- function(x, var1, padj = TRUE, include_test_name = TRUE) {
    res <- apply(x, 1, function(x) {
        # with ties an exact p-value can't be computed; R falls back to the normal
        # approximation and warns for each gene, so hide the warning (as run_wilcox did)
        test <- suppressWarnings(stats::wilcox.test(x ~ var1, exact = TRUE))
        pvalue <- test$p.value
        return(pvalue)
    })
    res <- make_res(res, padj, sort_columns = FALSE, include_test_name, test_name = "wilcox-test")
    return(res)
}

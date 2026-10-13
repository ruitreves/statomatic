#' Format per-gene p-values into a results table
#'
#' Shared by all new_* test functions. Two-group tests pass a named vector
#' (one p-value per gene) and get `pvalue`/`padj` columns; multi-group tests
#' pass a matrix (one column per comparison) and get
#' `pvalue_<comparison>`/`padj_<comparison>` columns.
#'
#' @param res named numeric vector, or matrix with genes as rows and comparisons as columns
#' @param padj logical, add Benjamini-Hochberg adjusted p-values
#' @param sort_columns logical, place each comparison's pvalue and padj columns side by side (only used when padj = TRUE)
#' @param include_test_name logical, add a `test` column
#' @param test_name name stored in the `test` column
#' @return a data.frame
#' @noRd
make_res <- function(res, padj, sort_columns, include_test_name, test_name = "") {
    if (is.null(dim(res))) {
        # two group tests
        res <- data.frame(pvalue = res, row.names = names(res))
    }
    else {
        # multi group tests
        res <- as.data.frame(res)
        colnames(res) <- paste0("pvalue_", colnames(res))
    }
    if (padj) {
        res_padj <- as.data.frame(lapply(res, stats::p.adjust, method = "BH"), check.names = FALSE)
        colnames(res_padj) <- sub("^pvalue", "padj", colnames(res_padj))
        res <- cbind(res, res_padj)
        if (sort_columns) {
            idx <- as.vector(matrix(seq_len(ncol(res)), nrow = 2, byrow = TRUE))
            res <- res[, idx]
        }
    }
    if (include_test_name && nzchar(test_name)) {
        res$test <- test_name
    }
    return(res)
}

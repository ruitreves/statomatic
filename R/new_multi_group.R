#' Tukey's HSD test on each row of a data.frame (rebuild of run_tukey)
#' @param x a numeric data.frame or matrix, genes as rows
#' @param var1 group label for each column of x
#' @param padj logical, add Benjamini-Hochberg adjusted p-values
#' @param include_test_name logical, add a `test` column
#' @param sort_columns logical, place each comparison's pvalue and padj columns side by side
#' @return a data.frame of p-values
#' @noRd
new_run_tukey <- function(x, var1, padj = TRUE, include_test_name = TRUE, sort_columns = FALSE) {
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
new_run_dunnett <- function(x, var1, padj = TRUE, include_test_name = TRUE, sort_columns = FALSE) {
    per_gene <- lapply(seq_len(nrow(x)), function(i) {
        d <- AdunnettT3Test(unlist(x[i, ]), g = as.factor(var1))$p.value
        cells <- which(!is.na(d), arr.ind=TRUE)
        labels <- paste0(colnames(d)[cells[, "col"]],
            "_vs_",
            rownames(d)[cells[, "row"]])
        stats::setNames(d[cells], labels)
    })
    res <- do.call(rbind, per_gene)
    rownames(res) <- rownames(x)
    res <- make_res(res, padj, sort_columns, include_test_name, test_name = "dunnett-test")
    return(res)
}

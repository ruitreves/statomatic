#' Fold changes between every pair of groups for each row of a data.frame (rebuild of fold_change)
#'
#' For each gene, computes the mean of each group and the log2 fold change for
#' every pair of groups, named `<A>_vs_<B>` (A relative to B).
#'   log_input = FALSE (data on its original scale): log2(|mean_A / mean_B|).
#'     If both means are 0 the result is 0; if only mean_B is 0 it is 1e99; if
#'     only mean_A is 0 it is -1e99.
#'   log_input = TRUE (data already log-transformed): mean_A - mean_B, in the
#'     units of the input's log base.
#' Missing values are ignored when computing group means.
#'
#' @param x a numeric data.frame or matrix, genes as rows
#' @param var1 group label for each column of x
#' @param log_input logical, TRUE if x is already log-transformed
#' @return a data.frame with one `log2FoldChange_<A>_vs_<B>` column per pair of
#'   groups, followed by one `<group>_mean` column per group
#' @noRd
new_fold_change <- function(x, var1, log_input = FALSE) {
    groups <- sort(unique(var1), method = "radix")
    means <- vapply(groups, function(g) rowMeans(x[, var1 == g, drop = FALSE], na.rm = TRUE),
                    numeric(nrow(x)))
    means <- matrix(means, nrow = nrow(x), dimnames = list(rownames(x), groups))   # keeps shape for one gene
    pairs <- utils::combn(groups, 2)

    fc <- apply(pairs, 2, function(p) {
        m1 <- means[, p[1]]
        m2 <- means[, p[2]]
        if (log_input) return(m1 - m2)
        out <- log2(abs(m1 / m2))
        out[m1 == 0 & m2 == 0] <- 0
        out[m1 != 0 & m2 == 0] <- 1e99
        out[m1 == 0 & m2 != 0] <- -1e99
        out
    })
    fc <- matrix(fc, nrow = nrow(x))   # keeps shape for one gene or one pair
    colnames(fc) <- paste0("log2FoldChange_", pairs[1, ], "_vs_", pairs[2, ])
    colnames(means) <- paste0(groups, "_mean")
    res <- as.data.frame(cbind(fc, means), check.names = FALSE)
    rownames(res) <- rownames(x)
    return(res)
}

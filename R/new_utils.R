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

#' Check and align the experimental factors for the multigroup workflow
#'
#' Matches `factors` to the samples (columns) of `x` by name, validates the
#' design, and builds the combined group used for post-hoc tests and fold
#' changes (e.g. condition = Control, sex = F -> "Control_F").
#'
#' @param x a numeric data.frame or matrix, genes as rows, samples as named columns
#' @param factors a named vector (one factor), or a data.frame with one column
#'   per factor (up to 3) and row names that match colnames(x)
#' @return a list: `factors`, a data.frame of factors (radix-sorted levels) in
#'   the column order of x, and `group`, the combined group as a factor
#' @noRd
check_factors <- function(x, factors) {
    show <- function(v, n = 5) {   # list up to n names in messages
        more <- if (length(v) > n) paste0(" and ", length(v) - n, " more") else ""
        paste0(paste(utils::head(v, n), collapse = ", "), more)
    }
    # --- match samples by name ---
    how_to <- paste0(
        "  `factors` needs row names that match colnames(x), e.g.\n",
        "      rownames(sample_info) <- sample_info$sample\n",
        "  For a single factor you can also pass a named vector, e.g.\n",
        "      setNames(sample_info$condition, sample_info$sample)")
    if (is.null(colnames(x))) {
        stop("`x` needs column names (sample IDs) so `factors` can be matched to its samples.", call. = FALSE)
    }
    if (is.null(dim(factors))) {
        ids <- names(factors)
        factors <- data.frame(var1 = unname(factors), stringsAsFactors = FALSE)
    } else {
        factors <- as.data.frame(factors)
        ids <- rownames(factors)
        if (identical(ids, as.character(seq_len(nrow(factors))))) ids <- NULL   # default row names
    }
    if (is.null(ids)) {
        stop("Can't match `factors` to the samples (columns) of `x`.\n", how_to, call. = FALSE)
    }
    if (anyDuplicated(ids)) {
        stop("Sample names in `factors` must be unique. Duplicated: ",
             show(unique(ids[duplicated(ids)])), call. = FALSE)
    }
    not_in_factors <- setdiff(colnames(x), ids)
    not_in_x <- setdiff(ids, colnames(x))
    if (length(not_in_factors) || length(not_in_x)) {
        stop("Samples in `x` and `factors` don't match.",
             if (length(not_in_factors)) paste0("\n  In x but not in factors: ", show(not_in_factors)),
             if (length(not_in_x)) paste0("\n  In factors but not in x: ", show(not_in_x)),
             call. = FALSE)
    }
    rownames(factors) <- ids
    factors <- factors[colnames(x), , drop = FALSE]   # same order as the columns of x

    # --- check the factors themselves ---
    if (ncol(factors) > 3) {
        stop("`factors` has ", ncol(factors), " columns; at most 3 factors are supported.", call. = FALSE)
    }
    has_na <- names(factors)[vapply(factors, anyNA, logical(1))]
    if (length(has_na)) {
        stop("`factors` has missing values in: ", paste(has_na, collapse = ", "), call. = FALSE)
    }
    is_num <- names(factors)[vapply(factors, is.numeric, logical(1))]
    if (length(is_num)) {
        warning("Numeric column(s) treated as categories: ", paste(is_num, collapse = ", "),
                ". Continuous covariates are not supported.", call. = FALSE)
    }
    factors[] <- lapply(factors, function(f) {
        f <- as.character(f)
        factor(f, levels = sort(unique(f), method = "radix"))
    })
    n_levels <- vapply(factors, nlevels, integer(1))
    if (ncol(factors) == 1 && n_levels < 3) {
        stop("A single factor needs at least 3 levels for the multigroup workflow; it has ",
             n_levels, ". Use twogroup_main for 2 groups.", call. = FALSE)
    }
    if (any(n_levels < 2)) {
        stop("Each factor needs at least 2 levels. Only one level in: ",
             paste(names(factors)[n_levels < 2], collapse = ", "), call. = FALSE)
    }

    # --- build the combined group and check the design ---
    labels <- do.call(paste, c(unname(as.list(factors)), sep = "_"))
    group <- factor(labels, levels = sort(unique(labels), method = "radix"))
    if (nlevels(group) != nlevels(interaction(factors, drop = TRUE))) {
        stop("Some factor levels produce the same combined label when joined with \"_\". ",
             "Rename levels so that the combinations are unique.", call. = FALSE)
    }
    if (nlevels(group) < prod(n_levels)) {
        all_cells <- do.call(paste, c(unname(as.list(expand.grid(lapply(factors, levels)))), sep = "_"))
        stop("The factors are not fully crossed, so interactions can't be estimated. ",
             "No samples for: ", show(setdiff(all_cells, levels(group))), call. = FALSE)
    }
    cell_n <- table(group)
    if (any(cell_n < 2)) {
        stop("Every group needs at least 2 samples. Fewer than 2 in: ",
             paste(names(cell_n)[cell_n < 2], collapse = ", "), call. = FALSE)
    }
    if (any(cell_n < 3)) {
        warning("Groups with only 2 samples: ", paste(names(cell_n)[cell_n < 3], collapse = ", "),
                ". Normality and variance checks will be very unreliable.", call. = FALSE)
    }
    list(factors = factors, group = group)
}

#' Attach fold changes to post-hoc p-values, matched by comparison name
#'
#' For each comparison in `pvals`, orders the columns as pvalue, padj (if
#' present), log2FoldChange, then appends the group means and any other
#' columns of `pvals` (e.g. test). Genes are matched by row name.
#'
#' @param pvals a data.frame of post-hoc results with `pvalue_<A>_vs_<B>` columns
#' @param fcs a data.frame from new_fold_change, covering at least the genes in pvals
#' @return a data.frame
#' @noRd
add_fold_changes <- function(pvals, fcs) {
    cmp <- sub("^pvalue_", "", grep("^pvalue_", colnames(pvals), value = TRUE))
    fc_cols <- paste0("log2FoldChange_", cmp)
    missing <- setdiff(fc_cols, colnames(fcs))
    if (length(missing)) stop("No fold changes for: ", paste(missing, collapse = ", "), call. = FALSE)
    fcs <- fcs[rownames(pvals), , drop = FALSE]
    per_cmp <- unlist(lapply(cmp, function(k) {
        c(paste0("pvalue_", k),
          if (paste0("padj_", k) %in% colnames(pvals)) paste0("padj_", k),
          paste0("log2FoldChange_", k))
    }))
    means <- grep("_mean$", colnames(fcs), value = TRUE)
    other <- setdiff(colnames(pvals), grep("^(pvalue|padj)_", colnames(pvals), value = TRUE))
    cbind(pvals, fcs)[, c(per_cmp, means, other), drop = FALSE]
}

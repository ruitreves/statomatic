#' The two-group statomatic workflow (rebuild of twogroup_main)
#'
#' Sorts genes by residual normality and equality of variance, then runs the
#' matching test on each set:
#'   ne (normal, equal variance):   Student's t-test
#'   nu (normal, unequal variance): Welch test
#'   nn (not normal):               Wilcoxon rank-sum test
#' Benjamini-Hochberg adjustment (padj) is done within each set.
#'
#' @param x a numeric data.frame or matrix, genes as rows, samples as named columns
#' @param sample_info a data.frame with one row per sample. Samples are matched
#'   to the columns of x by name: through `sample_col`, the row names, or the
#'   column whose values are the sample IDs (found automatically). Rows for
#'   samples that are not in x are ignored.
#' @param factors the name of the column of sample_info with the two groups,
#'   e.g. "condition"
#' @param sample_col optional name of the column of sample_info with the sample
#'   IDs, needed only if several columns contain them
#' @param padj logical, add Benjamini-Hochberg adjusted p-values
#' @param log_input logical, TRUE if x is already log-transformed (fold changes
#'   are then differences of group means)
#' @param include_test_name logical, add a `test` column. Defaults to TRUE
#'   because, unlike the multigroup workflow, there is no post-hoc test name to
#'   show which test a gene got.
#' @param shapiro_threshold significance level for the Shapiro-Wilk test
#' @param levene_threshold significance level for the Levene test
#' @param verbose logical, print a message for each step
#' @return a list:
#'   `ttest_results`, `welch_results`, `wilcox_results`: pvalue, padj,
#'     log2FoldChange, group means and test. A set with no genes gives a 0-row
#'     table with the same columns.
#'   `foldchanges`: fold changes and group means for all genes.
#'   `ne`, `nu`, `nn`: the subsets of x.
#'   `tests`: per-gene Shapiro and Levene p-values and the assigned set.
#' @noRd
new_twogroup_main <- function(x, sample_info, factors, sample_col = NULL,
                              padj = TRUE, log_input = FALSE, include_test_name = TRUE,
                              shapiro_threshold = 0.05, levene_threshold = 0.05, verbose = TRUE) {
    say <- function(...) if (verbose) message(...)

    checked <- check_factors(x, select_factors(x, sample_info, factors, sample_col), two_group = TRUE)
    group <- checked$group
    x <- as.data.frame(x)
    say("Analyzing ", nrow(x), " genes in 2 groups: ", paste(levels(group), collapse = ", "))

    say("Sorting genes by normality and equality of variance...")
    sets <- new_test_norm(x, checked$factors, shapiro_threshold = shapiro_threshold,
                          levene_threshold = levene_threshold)
    say("  ne: ", nrow(sets$ne), "  nu: ", nrow(sets$nu), "  nn: ", nrow(sets$nn))

    say("Calculating fold changes...")
    fcs <- new_fold_change(x, group, log_input = log_input)

    run_set <- function(dat, label, test) {
        run <- function(d) {
            p <- test(d)
            other <- setdiff(colnames(p), c("pvalue", "padj"))   # e.g. test
            cbind(p[, intersect(c("pvalue", "padj"), colnames(p)), drop = FALSE],
                  fcs[rownames(d), , drop = FALSE], p[, other, drop = FALSE])
        }
        if (nrow(dat) == 0) {
            # no genes: run on one gene to get the columns, then keep 0 rows
            return(run(x[1, , drop = FALSE])[0, , drop = FALSE])
        }
        say("Running ", label, " on ", nrow(dat), " genes...")
        run(dat)
    }
    res <- list(
        ttest_results  = run_set(sets$ne, "t-test",
                                 function(d) new_ttest(d, group, padj = padj, include_test_name = include_test_name)),
        welch_results  = run_set(sets$nu, "Welch test",
                                 function(d) new_run_welch(d, group, padj = padj, include_test_name = include_test_name)),
        wilcox_results = run_set(sets$nn, "Wilcoxon test",
                                 function(d) new_run_wilcox(d, group, padj = padj, include_test_name = include_test_name)),
        foldchanges = fcs,
        ne = sets$ne, nu = sets$nu, nn = sets$nn,
        tests = sets$tests)
    say("Done.")
    return(res)
}

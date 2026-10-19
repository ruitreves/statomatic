#' The multigroup statomatic workflow (rebuild of multigroup_main)
#'
#' Sorts genes by residual normality and equality of variance, then runs the
#' matching omnibus and post-hoc tests on each set:
#'   ne (normal, equal variance):   Type II ANOVA + Tukey HSD
#'   nu (normal, unequal variance): Welch test + Dunnett T3
#'   nn (not normal):               Kruskal-Wallis + Dunn
#' Post-hoc tests and fold changes compare every pair of groups, where a group
#' is one combination of the factors (e.g. condition_sex). The ANOVA tests every
#' main effect and interaction; Welch and Kruskal-Wallis are one-way tests on
#' the combined groups. Benjamini-Hochberg adjustment (padj) is done within
#' each set.
#'
#' @param x a numeric data.frame or matrix, genes as rows, samples as named columns
#' @param sample_info a data.frame with one row per sample. Samples are matched
#'   to the columns of x by name: through `sample_col`, the row names, or the
#'   column whose values are the sample IDs (found automatically). Rows for
#'   samples that are not in x are ignored.
#' @param factors names of the columns of sample_info to use as factors: one
#'   (with 3+ levels, one-way) to three (three-way), e.g. c("condition", "sex")
#' @param sample_col optional name of the column of sample_info with the sample
#'   IDs, needed only if several columns contain them
#' @param padj logical, add Benjamini-Hochberg adjusted p-values
#' @param log_input logical, TRUE if x is already log-transformed (fold changes
#'   are then differences of group means)
#' @param shapiro_threshold significance level for the Shapiro-Wilk test
#' @param levene_threshold significance level for the Levene test
#' @param verbose logical, print a message for each step
#' @return a list:
#'   `anova_results`, `welch_results`, `kw_results`: omnibus p-values, then for
#'     each comparison pvalue, padj and log2FoldChange, then group means and
#'     the post-hoc test name. A set with no genes gives a 0-row table with the
#'     same columns.
#'   `foldchanges`: fold changes and group means for all genes.
#'   `ne`, `nu`, `nn`: the subsets of x.
#'   `tests`: per-gene Shapiro and Levene p-values and the assigned set.
#' @noRd
new_multigroup_main <- function(x, sample_info, factors, sample_col = NULL,
                                padj = TRUE, log_input = FALSE,
                                shapiro_threshold = 0.05, levene_threshold = 0.05,
                                verbose = TRUE) {
    say <- function(...) if (verbose) message(...)

    checked <- check_factors(x, select_factors(x, sample_info, factors, sample_col))
    factors <- checked$factors
    group <- checked$group
    x <- as.data.frame(x)
    say("Analyzing ", nrow(x), " genes in ", nlevels(group), " groups: ",
        paste(levels(group), collapse = ", "))

    say("Sorting genes by normality and equality of variance...")
    sets <- new_test_norm(x, factors, shapiro_threshold = shapiro_threshold,
                          levene_threshold = levene_threshold)
    say("  ne: ", nrow(sets$ne), "  nu: ", nrow(sets$nu), "  nn: ", nrow(sets$nn))

    say("Calculating fold changes...")
    fcs <- new_fold_change(x, group, log_input = log_input)

    run_set <- function(dat, label, omnibus, posthoc) {
        run <- function(d) cbind(omnibus(d), add_fold_changes(posthoc(d), fcs))
        if (nrow(dat) == 0) {
            # no genes: run on one gene to get the columns, then keep 0 rows
            return(run(x[1, , drop = FALSE])[0, , drop = FALSE])
        }
        say("Running ", label, " on ", nrow(dat), " genes...")
        run(dat)
    }
    res <- list(
        anova_results = run_set(sets$ne, "ANOVA + Tukey",
                                function(d) new_run_anova(d, factors, padj = padj),
                                function(d) new_run_tukey(d, group, padj = padj)),
        welch_results = run_set(sets$nu, "Welch + Dunnett T3",
                                function(d) new_run_welch(d, group, padj = padj),
                                function(d) new_run_dunnett(d, group, padj = padj)),
        kw_results    = run_set(sets$nn, "Kruskal-Wallis + Dunn",
                                function(d) new_run_kw(d, group, padj = padj),
                                function(d) new_run_dunn(d, group, padj = padj)),
        foldchanges = fcs,
        ne = sets$ne, nu = sets$nu, nn = sets$nn,
        tests = sets$tests)
    say("Done.")
    return(res)
}

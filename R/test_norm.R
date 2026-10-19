#' Sort genes by residual normality and equality of variance
#'
#' For each gene, fits `y ~ factor1 * factor2 * ...` and tests the residuals for
#' normality (Shapiro-Wilk) and the groups for equal variance (Levene, median-
#' centered, groups = every combination of the factors). Genes are sorted into:
#'   ne: normal, equal variance    (Shapiro p >= shapiro_threshold, Levene p >= levene_threshold)
#'   nu: normal, unequal variance  (Shapiro p >= shapiro_threshold, Levene p <  levene_threshold)
#'   nn: not normal                (Shapiro p <  shapiro_threshold)
#' Missing values are dropped per gene. Shapiro-Wilk has little power with small
#' groups (~8 samples per group), so many non-normal genes will still be sorted
#' as normal. Constant genes should be removed beforehand.
#'
#' @param x a numeric data.frame or matrix, genes as rows
#' @param factors a vector with one group label per column of x (one-way), or a
#'   data.frame with one column per factor (two-way and up)
#' @param shapiro_threshold significance level for the Shapiro-Wilk test
#' @param levene_threshold significance level for the Levene test
#' @return a list with data.frames `ne`, `nu`, `nn` (subsets of x, original row
#'   order) and `tests`, a data.frame with per-gene `shapiro` and `levene`
#'   p-values and the assigned `set`. Levene is reported for every gene, but is
#'   only used to sort genes that pass Shapiro.
#' @export
test_norm <- function(x, factors, shapiro_threshold = 0.05, levene_threshold = 0.05) {
    if (is.null(dim(factors))) factors <- data.frame(var1 = factors)   # a single factor
    factors <- as.data.frame(factors)
    factors[] <- lapply(factors, factor)
    formu <- stats::as.formula(paste("y ~", paste(colnames(factors), collapse = " * ")))
    cells <- interaction(factors, drop = TRUE)

    per_gene <- lapply(seq_len(nrow(x)), function(i) {
        fit <- stats::lm(formu, data = cbind(y = unlist(x[i, ]), factors), na.action = stats::na.exclude)
        resids <- stats::residuals(fit)
        c(shapiro = stats::shapiro.test(resids)$p.value,
          levene = car::leveneTest(resids, cells)[["Pr(>F)"]][1])
    })
    tests <- as.data.frame(do.call(rbind, per_gene), row.names = rownames(x))
    tests$set <- factor(ifelse(tests$shapiro < shapiro_threshold, "nn",
                        ifelse(tests$levene < levene_threshold, "nu", "ne")),
                        levels = c("ne", "nu", "nn"))
    res <- split(as.data.frame(x), tests$set)
    res$tests <- tests
    return(res)
}

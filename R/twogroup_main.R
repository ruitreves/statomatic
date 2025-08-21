#' The two group statomatic procedure, pulling together all applicable functions in the correct order. Use this function to analyze a dataset 
#' with exactly two groups
#' @param x a numeric data.frame
#' @param var1 a list which maps columns of x to experimental factors. Must be exactly 2 factors
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @param write_files boolean, whether or not to write csv files of the test results instead of returning. Defaults to FALSE
#' @param file_names a list of three file names to be used when writing files. Specify as a vector using c(). The first will be the name of the t-test results,
#' the second the welch-test results, and the third the wilcox-test results. Only used if write_files = TRUE
#' @return a list containing all test results, fold changes, and subsets of x as sorted by test_norm, unless write_files = TRUE, in which case nothing is returned.
#' @export
twogroup_main <- function(x, var1, padj = FALSE, include_test_name = FALSE, write_files = FALSE, file_names = c("ttest_results.csv", "welch_results.csv", "wilcox_results.csv")) {
    if (write_files == TRUE && length(file_names) != 3) stop("This function creates three files but has not been provided three names")
    nenunn <- test_norm(x, var1)
    ne <- nenunn[[1]]
    nu <- nenunn[[2]]
    nn <- nenunn[[3]]
    ne_folds <- nu_folds <- nn_folds <- NULL
    t_res <- welch_res <- wilcox_res <- NULL
    if (nrow(ne) > 0) {
        t_res <- run_ttest(ne, var1, padj = padj, include_test_name)
        ne_folds <- fold_change(ne, var1)
        t_res <- as.data.frame(cbind(t_res, ne_folds))
    }
    if (nrow(nu) > 0) {
        welch_res <- run_welch(nu, var1, padj, include_test_name)
        nu_folds <- fold_change(nu, var1)
        welch_res <- as.data.frame(cbind(welch_res, nu_folds))
    }
    if (nrow(nn) > 0) {
        wilcox_res <- run_wilcox(nn, var1, padj, include_test_name)
        nn_folds <- fold_change(nn, var1)
        wilcox_res <- as.data.frame(cbind(wilcox_res, nn_folds))
    }

    if (write_files == TRUE) {
        write.csv(t_res, file_names[1])
        write.csv(welch_res, file_names[2])
        write.csv(wilcox_res, file_names[3])
        return(NULL)
    }
    all_fcs <- rbind(ne_folds, nu_folds, nn_folds)
    res <- list("ttest_results" = t_res, "welch_results" = welch_res, "wilcox_results" = wilcox_res, "foldchanges" = all_fcs,
                "ne" = ne, "nu" = nu, "nn" = nn)
    return(res)
}


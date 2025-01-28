#' Calculates fold changes, defaults to log2, for each row of a data frame according to specified factors
#' @param my_data a numeric data.frame
#' @param var1 a list which maps columns of my_data to experimental factors
#' FUN a function which can be passed one numeric value and return one numeric value. Default is log2, but could be log, identity, function(x) {return(x * 2)}, etc. 
#' @return a data.frame of fold changes
#' @export
fold_change <- function(my_data, var1, FUN = log2, ...) {
    if (!is.function(FUN)) stop(paste0(FUN, " is not a function"))
    cat("Calculating fold changes...")
    big_res <- NULL
    big_means_res <- NULL
    grp_list <- unique(var1)
    grp_list <- grp_list[order(grp_list)]
    for (i in rownames(my_data)) {
        x <- my_data[i, ]
        small_res <- NULL
        small_means_res <- NULL
        counter <- 1

        for (grp in grp_list) {
            counter <- counter + 1
            mean1 <- mean(unlist(x)[var1 == grp])
            small_means_res <- cbind(small_means_res, mean1)

            if (counter <= length(grp_list)) {
                for (j in counter:length(grp_list)) {
                    mean2 <- mean(unlist(x)[var1 == grp_list[j]])

                    if (mean2 == 0) {
                        if (mean1 == 0) {
                            fld_cng <- 0
                        }
                        else {
                            fld_cng <- 1e99
                        }
                    }
                    else if (mean1 == 0) {
                        fld_cng <- -1e99
                    }
                    else {
                        av <- abs(mean1 / mean2)
                        fld_cng <- FUN(av)
                    }

                    small_res <- cbind(small_res, fld_cng)
                    
                }
            }
        }
        big_means_res <- rbind(big_means_res, small_means_res)
        big_res <- rbind(big_res, small_res)
    }

    fc_name_list <- c()
    counter <- 1
    for (grp in grp_list) {
        counter <- counter + 1
        if (counter <= length(grp_list)) {
            for (j in counter:length(grp_list)) {
                fc_name <- paste("log2FoldChange_", grp, "_vs_", grp_list[j], sep="")
                fc_name_list <- append(fc_name_list, fc_name)
            }
        }
    }
    mean_names <- paste0(grp_list, "_mean")
    rownames(big_res) <- rownames(my_data)
    colnames(big_res) <- fc_name_list
    colnames(big_means_res) <- mean_names
    result <- cbind(big_res, big_means_res)
    cat("\n")
    return(result)
}

#' This function is a helper function for the statomatic workflow. It is not really designed to be used by itself. Combines columns of p value and fold change data.frames 
#' @param a a data.frame, these values will be on the left
#' @param b a data.frame, these values will be on the right
#' @export
combine_fc_pvals <- function(a, b) {
    out <- data.frame()
    pvals <- a[, -ncol(a), drop = FALSE]
    test <- a[, ncol(a), drop = FALSE]
    fcs <- b[, 1:ncol(pvals), drop = FALSE]
    means <- b[, -c(1:ncol(pvals)), drop = FALSE]
    for (i in 1:ncol(pvals)) {
        x <- cbind(pvals[, i, drop = FALSE], fcs[, i, drop = FALSE])
        out <- rbind(out, t(x))
    }
    out <- t(out)
    out <- cbind(out, means, test)
    return(out)
}

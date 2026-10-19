#' Set rownames to a specified column and remove that column from the dataframe
#' @param dat a data.frame
#' @param column_index the index of the column you want as your rownames, defaults to 1
#' @return the original data frame with the row names set as specified and the specified column removed
#' @export
cf <- function(dat, column_index = 1) {
    rownames(dat) <- dat[, column_index]
    dat[, column_index] <- NULL
    return(dat)
}

#' Produces a dummy data set from a normal distribution and sample_info object. 
#' @param group_min Minimum random number of groups. Defaults to 3
#' @param group_max Maximum random number of groups. Defaults to 20. Note: group_max * dim max must be <= 1000
#' @param dim_min Minimum random dimension of data frame. Defaults to 1
#' @param dim_max Maximum random dimension of data frame. Defaults to 50. Note: group_max * dim max must be <= 1000
#' @param param_min Minimum random parameters for normal distribution (i.e., mean and variance). Defaults to 0
#' @param param_max Maximum random parameters for normal distribution. Defaults to 20
#' @return A list, element one is the randomly generated data frame, element two is sample_info. Access with double brackets: `[[1]]`, `[[2]]`.
#' @export 

dummy_data <- function(group_min = 2, group_max = 20, dim_min = 1, dim_max = 50, param_min = 0, param_max = 20) {
    group_num <- floor(stats::runif(1, min = group_min, max = group_max))

    group_names <- c()
    for (i in 1:group_num) {
        group_name <- paste0("group", i)
        group_names <- append(group_names, group_name)
    }

    dims <- ceiling(stats::runif(2, min = dim_min, max = dim_max))
    col_num <- dims[1] * group_num
    row_num <- dims[2]

    parms <- ceiling(stats::runif(2, min = param_min, max = param_max))
    m <- parms[2]
    v <- parms[1]

    x <- data.frame()
    cols <- c()
    for (i in 1:col_num) {
        n <- stats::rnorm(row_num, m, v)
        c <- paste0("Y", i)
        cols <- append(cols, c)
        x <- rbind(x, n)
    }

    x <- t(x)

    rows <- c()
    for (i in 1:row_num) {
        r <- paste0("X", i)
        rows <- append(rows, r)
    }

    rownames(x) <- rows
    colnames(x) <- cols

    sample_info <- data.frame(samples = colnames(x), groups = group_names)

    dum <- c()
    dum[[1]] <- x
    dum[[2]] <- sample_info

    return(dum)
}

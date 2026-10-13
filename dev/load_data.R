load_data <- function(counts = "dev/counts.csv", sample_info = "dev/sample_info.csv", min_count = 10, min_samples = 10) {
    x <- read.csv(counts, row.names = 1, check.names = FALSE)
    si <- read.csv(sample_info)
    
    missing <- setdiff(si$sample, colnames(x))
    if (length(missing)) stop("Samples missing from counts: ", paste(missing, collapse = ", "))

    x <- x[, si$sample, drop = FALSE]                                   # select + order columns by name
    x <- x[rowSums(x >= min_count) >= min_samples, ]
    list(x = x, sample_info = si)
}
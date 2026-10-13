source("dev/load_data.R")

d <- load_data()

x <- d$x
x <- x[1:10, ]
sample_info <- d$sample_info

### 
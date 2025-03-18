d <- dummy_data()
x <- d[[1]]
sample_info <- d[[2]]

f <- fold_change(x, sample_info$groups)

group_list <- unique(sample_info$groups)
group_list <- group_list[order(group_list)]

sample_info <- sample_info[order(sample_info$groups), ]

x <- x[, match(sample_info$samples, colnames(x))]

for (i in 1:nrow(x)) {
  test <- NULL
  for (g in group_list) {
    temp <- mean(x[i, ][sample_info$groups == g])
    test <- cbind(test, temp)
  }

  a <- nrow(sample_info) / length(unique(sample_info$groups))

  counter <- 1
  means <- NULL
  for (j in 1:length(group_list)) {
    temp <- mean(x[i, counter:(a * j)])
    means <- cbind(means, temp)
    counter <- counter + a
  }

  test_that("means calculated correctly", {
    expect_equal(test, means)
  })
}

#n = length(group_list)
#ncol = sum((n - 1):1)
#first iteration has n-1 elements, then n-2, n-3, ..., 1 iteration

for (i in 1:nrow(x)) {

  test <- NULL
  for (g in group_list) {
    temp <- mean(x[i, ][sample_info$groups == g])
    test <- cbind(test, temp)
  }

  y <- f[i, ]

  means <- c()
  counter <- 1
  for (j in 1:length(test)) {
    mean1 <- test[j]
    counter <- counter + 1
    if (counter <= length(test)) {
      for (k in counter:length(test)) {
        mean2 <- test[k]
        if (mean2 == 0) {
          r <- 1e99
        }
        else if (mean1 == 0) {
          r <- -1e99
        }
        else {
          r <- log2(abs(mean1 / mean2))
        }
        means <- cbind(means, r)
      }
    }
  }
  for (l in 1:length(means)) {
    test_that("fc works", {
    expect_equal(y[[l]], means[l])
    })
  }
}


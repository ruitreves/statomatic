# Welch test on each row of a data.frame

Welch test on each row of a data.frame

## Usage

``` r
run_welch(x, var1, padj = TRUE, include_test_name = FALSE)
```

## Arguments

- x:

  a numeric data.frame or matrix, genes as rows

- var1:

  group label for each column of x

- padj:

  logical, add Benjamini-Hochberg adjusted p-values

- include_test_name:

  logical, add a `test` column

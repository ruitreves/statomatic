# Student's t-test on each row of a data.frame

Student's t-test on each row of a data.frame

## Usage

``` r
run_ttest(x, var1, padj = TRUE, include_test_name = TRUE)
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

## Value

a data.frame of p-values

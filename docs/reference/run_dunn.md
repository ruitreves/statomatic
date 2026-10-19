# Dunn's test on each row of a data.frame

Dunn's test on each row of a data.frame

## Usage

``` r
run_dunn(x, var1, padj = TRUE, include_test_name = TRUE, sort_columns = TRUE)
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

- sort_columns:

  logical, place each comparison's pvalue and padj columns side by side

# Type II ANOVA on each row of a data.frame

Fits `y ~ factor1 * factor2 * ...` for each gene and tests every main
effect and interaction with a Type II ANOVA (`car::Anova(type = 2)`), so
results do not depend on the order of the factors. Any number of factors
is supported.

## Usage

``` r
run_anova(
  x,
  factors,
  padj = TRUE,
  include_test_name = FALSE,
  sort_columns = TRUE
)
```

## Arguments

- x:

  a numeric data.frame or matrix, genes as rows

- factors:

  a vector with one group label per column of x (one-way), or a
  data.frame with one column per factor (two-way and up). Column names
  are used to name the result columns.

- padj:

  logical, add Benjamini-Hochberg adjusted p-values

- include_test_name:

  logical, add a `test` column

- sort_columns:

  logical, place each term's pvalue and padj columns side by side

## Value

a data.frame of p-values. One-way: `pvalue`/`padj`. Two-way and up:
`pvalue_<term>`/`padj_<term>` for each main effect and interaction.

# Fold changes between every pair of groups for each row of a data.frame

For each gene, computes the mean of each group and the log2 fold change
for every pair of groups, named `<A>_vs_<B>` (A relative to B).
log_input = FALSE (data on its original scale): log2(\|mean_A /
mean_B\|). If both means are 0 the result is 0; if only mean_B is 0 it
is 1e99; if only mean_A is 0 it is -1e99. log_input = TRUE (data already
log-transformed): mean_A - mean_B, in the units of the input's log base.
Missing values are ignored when computing group means.

## Usage

``` r
fold_change(x, var1, log_input = FALSE)
```

## Arguments

- x:

  a numeric data.frame or matrix, genes as rows

- var1:

  group label for each column of x

- log_input:

  logical, TRUE if x is already log-transformed

## Value

a data.frame with one `log2FoldChange_<A>_vs_<B>` column per pair of
groups, followed by one `<group>_mean` column per group

# Attach fold changes to post-hoc p-values, matched by comparison name

Combines the results of a pairwise test (run_tukey, run_dunnett,
run_dunn) with the output of fold_change, as multigroup_main does. For
each comparison in `pvals`, the columns are ordered pvalue, padj (if
present), log2FoldChange; then come the group means and any other
columns of `pvals` (e.g. test). Rows are matched by row name, and
comparisons by name, so the order of rows and columns in the two inputs
doesn't matter.

## Usage

``` r
add_fold_changes(pvals, fcs)
```

## Arguments

- pvals:

  a data.frame of pairwise results with `pvalue_<A>_vs_<B>` columns

- fcs:

  a data.frame from fold_change, with a row for every row of pvals
  (extra rows are ignored)

## Value

a data.frame with the rows of pvals in their original order

## Examples

``` r
x <- matrix(c(5, 6, 7, 10, 11, 12, 20, 21, 22,
              6, 6, 8, 9, 11, 12, 19, 22, 24), nrow = 2, byrow = TRUE,
            dimnames = list(c("gene1", "gene2"), paste0("s", 1:9)))
group <- rep(c("A", "B", "C"), each = 3)
add_fold_changes(run_tukey(x, group), fold_change(x, group))
#>       pvalue_A_vs_B padj_A_vs_B log2FoldChange_A_vs_B pvalue_A_vs_C
#> gene1   0.002101241 0.004202481            -0.8744691  3.781723e-06
#> gene2   0.080882786 0.080882786            -0.6780719  1.377579e-04
#>        padj_A_vs_C log2FoldChange_A_vs_C pvalue_B_vs_C  padj_B_vs_C
#> gene1 7.563447e-06             -1.807355  4.450051e-05 8.900102e-05
#> gene2 1.377579e-04             -1.700440  7.764890e-04 7.764890e-04
#>       log2FoldChange_B_vs_C   A_mean   B_mean   C_mean       test
#> gene1            -0.9328858 6.000000 11.00000 21.00000 tukey-test
#> gene2            -1.0223678 6.666667 10.66667 21.66667 tukey-test
```

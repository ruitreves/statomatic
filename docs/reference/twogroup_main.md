# The two-group statomatic workflow

Sorts genes by residual normality and equality of variance, then runs
the matching test on each set: ne (normal, equal variance): Student's
t-test nu (normal, unequal variance): Welch test nn (not normal):
Wilcoxon rank-sum test Benjamini-Hochberg adjustment (padj) is done
within each set.

## Usage

``` r
twogroup_main(
  x,
  sample_info,
  factors,
  sample_col = NULL,
  padj = TRUE,
  log_input = FALSE,
  include_test_name = TRUE,
  shapiro_threshold = 0.05,
  levene_threshold = 0.05,
  verbose = TRUE
)
```

## Arguments

- x:

  a numeric data.frame or matrix, genes as rows, samples as named
  columns

- sample_info:

  a data.frame with one row per sample. Samples are matched to the
  columns of x by name: through `sample_col`, the row names, or the
  column whose values are the sample IDs (found automatically). Rows for
  samples that are not in x are ignored.

- factors:

  the name of the column of sample_info with the two groups, e.g.
  "condition"

- sample_col:

  optional name of the column of sample_info with the sample IDs, needed
  only if several columns contain them

- padj:

  logical, add Benjamini-Hochberg adjusted p-values

- log_input:

  logical, TRUE if x is already log-transformed (fold changes are then
  differences of group means)

- include_test_name:

  logical, add a `test` column. Defaults to TRUE because, unlike the
  multigroup workflow, there is no post-hoc test name to show which test
  a gene got.

- shapiro_threshold:

  significance level for the Shapiro-Wilk test

- levene_threshold:

  significance level for the Levene test

- verbose:

  logical, print a message for each step

## Value

a list: `ttest_results`, `welch_results`, `wilcox_results`: pvalue,
padj, log2FoldChange, group means and test. A set with no genes gives a
0-row table with the same columns. `foldchanges`: fold changes and group
means for all genes. `ne`, `nu`, `nn`: the subsets of x. `tests`:
per-gene Shapiro and Levene p-values and the assigned set.

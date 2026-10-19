# statomatic 1.0.0

A rebuild of statomatic. The statistical approach is unchanged: each row is
sorted by the normality of its residuals and the equality of its group
variances, and then analyzed with the matching tests. For one-factor designs,
the results are the same as in earlier versions (tested against saved results
of the previous version). This release has breaking changes.

## Breaking changes

* `multigroup_main()` and `twogroup_main()` take the sample table and the names
  of its factor columns: `multigroup_main(x, sample_info, factors = c("condition", "sex"))`.
  This replaces the `var1`, `var2` and `var3` arguments.
* Samples are matched to the columns of `x` by name, through a sample ID column
  of `sample_info` (found automatically), its row names, or `sample_col`.
  The analysis stops with an error if samples can't be matched, instead of
  assuming that the order is the same.
* For two and three factors, the pairwise tests (Tukey, Dunnett T3, Dunn), the
  Welch and Kruskal-Wallis tests and the fold changes use the combinations of
  all factors as groups (e.g. `Female_KO`), built automatically. Previously
  they used one factor (`var3`, by default the first factor).
* ANOVA uses Type II sums of squares (`car::Anova()`), so results don't depend
  on the order of the factors. With two or more factors and unequal group
  sizes, the main-effect p-values differ from earlier versions (Type I).
* Result columns are named `pvalue_*`, `padj_*` and `log2FoldChange_*`
  (previously `Pval_*` and `Padj_*`). In each result table, the columns of each
  comparison are next to each other: pvalue, padj, log2FoldChange.
* `padj = TRUE` is the default for all functions (previously `FALSE`).
* The test functions (`run_ttest()`, `run_tukey()`, ...) still take the data
  and the group labels, but the data argument is now called `x` (previously
  `my_data`), and all of them have `padj` and `include_test_name` arguments.
  `test_norm()` and `run_anova()` take the factors as a vector (one factor) or
  a data.frame (several factors) instead of `var1` and `var2`.
* `test_norm()` returns a named list (`ne`, `nu`, `nn`, `tests`).
* `fold_change()` no longer has a `FUN` argument.

## Removed

* The S4 interface: `make_sds()`, `sds_analyze()`, the `sds_*()` methods and the
  `get_*()` accessors. Use `multigroup_main()` and `twogroup_main()` directly.
  The SummarizedExperiment dependency is removed with it.
* `combine_fc_pvals()`, `chop()` and `takeN()`.
* `Adunn.test()`, `AdunnettT3Test()` and `p.adjustment.methods`, modified copies
  of `dunn.test::dunn.test()` and `PMCMRplus::dunnettT3Test()`. statomatic now
  has its own implementations of Dunn's test and Dunnett's T3 test, which give
  the same results as the copies (current dunn.test versions give different
  p-values than the 1.3.6 copy). The PMCMRplus dependency is removed.
* The `write_files` and `file_names` arguments of the main functions. Use
  `write.csv()` on the results.

## New features

* `test_norm()` and the main functions return a `tests` table with each row's
  Shapiro-Wilk and Levene p-values and the set it was sorted into.
* `shapiro_threshold` and `levene_threshold` set the significance levels used
  to sort rows (both 0.05 by default).
* `log_input = TRUE` in `fold_change()` and the main functions, for data that
  is already log-transformed: fold changes are then differences of group means.
* `run_anova()` and `test_norm()` support any number of factors (the main
  functions support up to three).
* The design is checked before the analysis: missing samples, missing values,
  factors that are not fully crossed, groups with fewer than 2 samples, and
  more than 3 factors give clear errors.
* `verbose = TRUE` prints progress messages in the main functions.
* A set without rows gives a 0-row results table with the usual columns.

## Other changes

* Missing values are dropped per row in the tests and ignored in group means.
* Comparison names (`A_vs_B`) and their order are the same in every locale.
* Groups with identical constant values in Dunnett's T3 test give p = 1
  instead of an error, as in earlier versions.
* From R 4.6.0, `wilcox.test()` computes exact p-values for data with ties, so
  Wilcoxon p-values for tied data differ between R versions.

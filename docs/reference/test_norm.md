# Sort genes by residual normality and equality of variance

For each gene, fits `y ~ factor1 * factor2 * ...` and tests the
residuals for normality (Shapiro-Wilk) and the groups for equal variance
(Levene, median- centered, groups = every combination of the factors).
Genes are sorted into: ne: normal, equal variance (Shapiro p \>=
shapiro_threshold, Levene p \>= levene_threshold) nu: normal, unequal
variance (Shapiro p \>= shapiro_threshold, Levene p \< levene_threshold)
nn: not normal (Shapiro p \< shapiro_threshold) Missing values are
dropped per gene. Shapiro-Wilk has little power with small groups (~8
samples per group), so many non-normal genes will still be sorted as
normal. Constant genes should be removed beforehand.

## Usage

``` r
test_norm(x, factors, shapiro_threshold = 0.05, levene_threshold = 0.05)
```

## Arguments

- x:

  a numeric data.frame or matrix, genes as rows

- factors:

  a vector with one group label per column of x (one-way), or a
  data.frame with one column per factor (two-way and up)

- shapiro_threshold:

  significance level for the Shapiro-Wilk test

- levene_threshold:

  significance level for the Levene test

## Value

a list with data.frames `ne`, `nu`, `nn` (subsets of x, original row
order) and `tests`, a data.frame with per-gene `shapiro` and `levene`
p-values and the assigned `set`. Levene is reported for every gene, but
is only used to sort genes that pass Shapiro.

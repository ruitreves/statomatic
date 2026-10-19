# Approach and limitations

This vignette explains the reasoning behind statomatic’s approach, and
what to keep in mind when using it. It includes simulations that measure
how well statomatic sorts rows, and what happens when a row is sorted
into the wrong set.

## The approach

Tests like the ANOVA rely on assumptions about the data: the residuals
should be normally distributed and the groups should have equal
variances. In high-throughput data, thousands of rows (genes, taxa,
metabolites, …) are tested at once, and some rows will meet these
assumptions while others won’t. Applying one test to every row ignores
this.

Statomatic checks the assumptions for each row separately, and chooses a
test that suits that row:

| Set | Meaning | 2 groups | \> 2 groups | \> 2 groups: pairwise comparisons |
|----|----|----|----|----|
| ne | normal, equal variance | t-test | ANOVA | Tukey HSD |
| nu | normal, unequal variance | Welch test | Welch test | Dunnett T3 |
| nn | not normal | Wilcoxon test | Kruskal-Wallis test | Dunn test |

A row is sorted into nn if the Shapiro-Wilk test on its residuals has p
\< 0.05 (`shapiro_threshold`). Otherwise, it is sorted into nu if
Levene’s test has p \< 0.05 (`levene_threshold`), and into ne if not.

The `tests` element of the results shows the Shapiro-Wilk and Levene
p-values of every row and the set it was sorted into, so you can see how
each row was handled.

## How well does the sorting work?

The Shapiro-Wilk and Levene tests are themselves statistical tests. How
often they detect non-normality or unequal variances depends on the
number of samples: with few samples per group, they often miss them.

To measure this, we simulated rows for which the right set is known, and
sorted them with
[`test_norm()`](https://ruitreves.github.io/statomatic/reference/test_norm.md).
Each scenario has three groups of the same size and 500 simulated rows:

- **normal, equal variance**: normally distributed, all groups with SD 1
  (should be sorted into ne)
- **normal, unequal variance**: normally distributed, one group with an
  SD 3 or 5 times larger (should be sorted into nu)
- **skewed**: log-normally distributed (should be sorted into nn)

The table shows the percentage of rows sorted into the right set:

| Scenario | 4 per group | 5 per group | 6 per group | 8 per group | 12 per group | 20 per group |
|:---|---:|---:|---:|---:|---:|---:|
| normal, equal variance | 90 | 94 | 92 | 90 | 88 | 92 |
| normal, unequal variance (SD 3:1) | 23 | 11 | 31 | 40 | 38 | 22 |
| normal, unequal variance (SD 5:1) | 29 | 20 | 29 | 25 | 7 | 1 |
| skewed | 31 | 42 | 57 | 78 | 96 | 100 |

Percentage of rows sorted into the right set

Where the rows that were sorted into the wrong set went:

| Scenario | Samples per group | % sorted into ne | % sorted into nu | % sorted into nn |
|:---|---:|---:|---:|---:|
| normal, unequal variance (SD 3:1) | 4 | 60 | 23 | 17 |
| normal, unequal variance (SD 3:1) | 5 | 70 | 11 | 20 |
| normal, unequal variance (SD 3:1) | 6 | 42 | 31 | 27 |
| normal, unequal variance (SD 3:1) | 8 | 26 | 40 | 35 |
| normal, unequal variance (SD 3:1) | 12 | 7 | 38 | 55 |
| normal, unequal variance (SD 3:1) | 20 | 1 | 22 | 77 |
| normal, unequal variance (SD 5:1) | 4 | 32 | 29 | 39 |
| normal, unequal variance (SD 5:1) | 5 | 37 | 20 | 43 |
| normal, unequal variance (SD 5:1) | 6 | 15 | 29 | 55 |
| normal, unequal variance (SD 5:1) | 8 | 2 | 25 | 73 |
| normal, unequal variance (SD 5:1) | 12 | 0 | 7 | 92 |
| normal, unequal variance (SD 5:1) | 20 | 0 | 1 | 99 |
| skewed | 4 | 59 | 10 | 31 |
| skewed | 5 | 57 | 1 | 42 |
| skewed | 6 | 36 | 7 | 57 |
| skewed | 8 | 19 | 3 | 78 |
| skewed | 12 | 2 | 1 | 96 |
| skewed | 20 | 0 | 0 | 100 |

Two things stand out:

1.  **With few samples, most rows are sorted into ne.** Normal rows with
    equal variance are sorted correctly about 90% of the time at every
    group size, which is expected: with a threshold of 0.05, each of the
    two checks wrongly flags about 5% of these rows. But with 4 to 6
    samples per group, which is common, many skewed rows and many rows
    with unequal variances are also sorted into ne, because the checks
    don’t have enough samples to detect the problem.

2.  **Rows with unequal variances are often sorted into nn.** The
    Shapiro-Wilk test is applied to the residuals of all groups
    together. When one group is much more spread out than the others,
    these pooled residuals have heavier tails than a normal
    distribution, so the row is sorted into nn before Levene’s test is
    reached, even though each group is normally distributed. This
    happens more often with more samples, because the Shapiro-Wilk test
    then detects the heavier tails more easily.

## What happens when a row is sorted into the wrong set?

A row in the wrong set is analyzed with a test that doesn’t suit it. To
see what this means, we simulated rows with **equal group means**, so
that every p-value below 0.05 is a false positive, and tested them with
[`run_anova()`](https://ruitreves.github.io/statomatic/reference/run_anova.md),
[`run_welch()`](https://ruitreves.github.io/statomatic/reference/run_welch.md)
and
[`run_kruskal()`](https://ruitreves.github.io/statomatic/reference/run_kruskal.md).
The third group has an SD of 1, 3 or 5 times that of the other groups.
Each scenario has 5,000 simulated rows. The false-positive rate should
be 5%:

| Group sizes | Group SDs | ANOVA (%) | Welch (%) | Kruskal-Wallis (%) |
|:------------|:----------|----------:|----------:|-------------------:|
| 4/4/4       | 1/1/1     |       5.5 |       4.5 |                4.3 |
| 4/4/4       | 1/1/3     |       7.9 |       4.4 |                6.1 |
| 4/4/4       | 1/1/5     |      10.8 |       5.5 |                8.8 |
| 6/6/6       | 1/1/1     |       4.5 |       4.2 |                4.0 |
| 6/6/6       | 1/1/3     |       7.7 |       4.4 |                5.0 |
| 6/6/6       | 1/1/5     |       9.7 |       5.6 |                6.1 |
| 6/5/4       | 1/1/1     |       5.0 |       4.5 |                4.1 |
| 6/5/4       | 1/1/3     |      13.2 |       5.5 |                7.9 |
| 6/5/4       | 1/1/5     |      16.8 |       5.9 |               10.3 |
| 12/8/5      | 1/1/1     |       5.2 |       5.4 |                5.0 |
| 12/8/5      | 1/1/3     |      18.2 |       5.2 |                8.2 |
| 12/8/5      | 1/1/5     |      22.8 |       5.4 |                9.3 |

False-positive rate at p \< 0.05, with equal group means

- **The Welch test** stays close to 5% in every scenario, as it is
  designed to.
- **The ANOVA** has too many false positives when the variances are
  unequal: up to 23% when the group with the largest variance is also
  the smallest. This is why rows with unequal variances should not be
  analyzed by an ANOVA, and it is what happens to the unequal-variance
  rows that are sorted into ne.
- **The Kruskal-Wallis test** also has too many false positives when the
  variances are unequal (up to 10%), but fewer than the ANOVA. It
  assumes that the groups have the same shape and spread, so it is not a
  safe choice for rows with unequal variances either. This affects the
  unequal-variance rows that are sorted into nn.

## The cost of a more robust test

When the assumptions of the ANOVA are met, it detects real differences
more often than the Welch and Kruskal-Wallis tests. We simulated rows
with equal variances in which the third group’s mean is 2 SD higher, and
measured how often each test detects the difference (p \< 0.05):

| Samples per group | ANOVA (%) | Welch (%) | Kruskal-Wallis (%) |
|------------------:|----------:|----------:|-------------------:|
|                 4 |      69.0 |      57.8 |               59.9 |
|                 5 |      82.7 |      75.5 |               77.1 |
|                 6 |      90.0 |      85.5 |               85.6 |
|                 8 |      97.5 |      96.3 |               96.5 |

Power with equal variances (third group mean 2 SD higher)

With 4 samples per group, the ANOVA detects the difference about 10
percentage points more often than the other two tests. The difference
becomes small with 8 or more samples per group. This trade-off between
power and robustness is the reason for choosing the test per row.

## Other things to keep in mind

- **Interactions.** Interaction p-values come from the ANOVA, so they
  are only available for rows sorted into ne. The Welch and
  Kruskal-Wallis tests are one-way tests on the combined groups
  (e.g. Female_KO, Female_WT, Male_KO, Male_WT).
- **Two kinds of p-value adjustment.** The pairwise p-values (Tukey,
  Dunnett T3, Dunn) are adjusted for comparing several pairs of groups
  within one row. The `padj` columns additionally adjust for testing
  many rows (Benjamini-Hochberg), separately within each set.
- **Reproducibility of Dunnett T3.** Its p-values are calculated with a
  random (Monte Carlo) method and can change very slightly between runs.
  Use [`set.seed()`](https://rdrr.io/r/base/Random.html) before the
  analysis to make them reproducible.
- **Wilcoxon p-values and R versions.** From R 4.6.0,
  [`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) calculates
  exact p-values for data with ties; earlier versions use an
  approximation. Wilcoxon p-values for data with ties (common in count
  data) therefore depend on the R version. Report the R version you
  used.
- **Rows with constant values.** Rows where every value is the same
  can’t be tested and should be removed before the analysis.
- **Specialized tools.** For RNA-seq count data, tools such as DESeq2,
  edgeR and limma model the properties of count data and share
  information between genes, which matters most with few samples. They
  also give interaction p-values for every gene. Statomatic is not a
  replacement for them.

## Practical suggestions

- Look at `res$tests` and the number of rows in each set. With few
  samples per group, keep in mind that the sorting is uncertain.
- More samples per group make the sorting more reliable.
- Report which tests were applied, the thresholds, the R version and the
  seed in your methods.

## Reproducing the simulations

The simulations were run with statomatic 1.0.0 (R version 4.6.1
(2026-06-24)). The script is installed with the package:

``` r
system.file("simulations", "sorting_simulations.R", package = "statomatic")
```

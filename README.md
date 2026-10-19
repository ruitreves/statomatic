# statom@ic (statomatic)
This package was developed to provide a straight-forward and high-throughput method to perform statistical analyses on numeric data.

## Citation
Treves, R.S.; Gripshover, T.C.; Hardesty, J.E. STATom@ic: R Package for Automated Statistical Analysis of Omic Datasets. Stats 2025, 8, 18. https://doi.org/10.3390/stats8010018

In R, `citation("statomatic")` gives the same reference.

## Installation
Statomatic is available from this GitHub site. Install it with the remotes package:

```r
install.packages("remotes")
remotes::install_github("ruitreves/statomatic")
```

All of statomatic's dependencies are on CRAN and are installed automatically.

## Introduction
Certain experimental designs are often suited to the preplanned statistical analysis that will take place after data is collected. It is vital, however, to apply the appropriate statistical techniques in order to ensure the validity of the results. The purpose of this package is to provide a high-throughput way to generate statistical results based on characteristics of the data itself, rather than on assumptions about how certain types of data "usually" are.

As an example, consider a researcher who hopes to perform a two-way ANOVA on a dataset. The researcher may plug their data in and generate statistical results, but without consideration for the mathematical assumptions of the two-way ANOVA, their results could be erroneous. To be suitable for a two-way ANOVA, the residuals of the data must be normally distributed and the data must be homoscedastic. The statomatic package is equipped to address these assumptions and sort the data based on which tests can be appropriately applied, and apply them.

### Background
The statomatic package is specifically designed to handle many separate hypothesis tests where the aim is to look for statistical differences between the observations (or groups of observations) for a given variable. For instance, for some gene X (the row), one could be interested in how it is expressed by individuals (the columns) subject to different experimental factors. This is just one hypothesis test, and in this scenario one would likely want to test for significance for several thousand genes. This is the purpose of statomatic, to be a high-throughput hypothesis testing package.

What sets statomatic apart from other similar programs is its ability to sort data based on which statistical test is mathematically appropriate for the variable being tested. This means that a different statistical test may be applied to different variables in the dataset. Statomatic does not apply a "one-size-fits-all" approach to an entire data set, but rather handles each individual variable separately.

See the details section for in-depth information on how statomatic sorts and analyzes data.

## Details

The aim of the program is to apply the appropriate statistical test to each variable in a data set. Each row is assumed to be a variable, and columns are assumed to be observations of that variable.

1. For each row in a data set, a linear model is fitted with the specified factor(s) (groups, etc.) and all their interactions.
   The residuals are tested with a Shapiro-Wilk normality test (alpha = 0.05 by default, `shapiro_threshold`).
   Rows with normally distributed residuals are then tested with Levene's test for equality of variances (homoscedasticity) between the groups (alpha = 0.05 by default, `levene_threshold`).
   - Rows that are both normally distributed and homoscedastic form the set "ne" (for normal, equal variance).
   - Rows that are normally distributed but heteroscedastic form the set "nu" (for normal, unequal variance).
   - Rows that are not normally distributed form the set "nn" (for non-normal).

2. In the case of multiple groups (> 2 groups), with `multigroup_main()`:
   - ne is analyzed by a one-way, two-way or three-way ANOVA, depending on the number of experimental factors (Type II sums of squares, so the order of the factors doesn't matter). A p-value is reported for each factor and each interaction. A Tukey multiple-comparison test follows, and a p-value is reported for each pairwise comparison of groups.
   - nu is analyzed by a Welch test, followed by a Dunnett T3 multiple-comparison test.
   - nn is analyzed by a Kruskal-Wallis rank sum test, followed by a Dunn multiple-comparison test.

   With more than one factor, the groups for the pairwise comparisons, the Welch and Kruskal-Wallis tests and the fold changes are the combinations of the factors, which statomatic builds automatically (e.g. sex and genotype give the groups Female_KO, Female_WT, Male_KO and Male_WT). With one factor, its levels are the groups.

   In the case of two groups, with `twogroup_main()`:
   - ne is analyzed by a Student's t-test.
   - nu is analyzed by a Welch's t-test.
   - nn is analyzed by a Wilcoxon/Mann-Whitney test.

3. P-values are adjusted for testing many rows with the Benjamini-Hochberg method (`padj` columns, within each set). The pairwise p-values (Tukey, Dunnett T3, Dunn) are in addition already adjusted for comparing several pairs of groups within a row.

4. Log2 fold changes are calculated between every pair of groups. The groups are sorted alphabetically (by character code, so upper case comes before lower case). The fold changes are reported as `log2FoldChange_GroupA_vs_GroupB`, which means log2(mean(GroupA) / mean(GroupB)). For data that is already log-transformed, use `log_input = TRUE`: the fold change is then mean(GroupA) - mean(GroupB).

5. Your results are ready.

## Vignettes

Basic workflow vignette: https://ruitreves.github.io/statomatic/articles/basic_workflow.html

Additional workflows vignette: https://ruitreves.github.io/statomatic/articles/additional_workflows.html

Approach and limitations: https://ruitreves.github.io/statomatic/articles/approach_and_limitations.html

## Quick Start Guide

Load the package:

```r
library(statomatic)
```

Read in data, with the first column (the gene ids) as row names:

```r
x <- read.csv("normalized_counts.csv", row.names = 1)
head(x)
```

Statomatic's `cf()` does the same for data that is already loaded: `x <- cf(read.csv("normalized_counts.csv"))`.

Create sample metadata to define the experimental factors, with one row per sample:

```r
sample_info <- data.frame(sample = colnames(x),
                          sex = c(rep("m", 10), rep("f", 10)),
                          genotype = c(rep("wt", 5), rep("ko", 5), rep("wt", 5), rep("ko", 5)))
sample_info
```

Samples are matched to the columns of `x` by name, using the column of `sample_info` that holds the sample names (found automatically) or its row names.

Run the analysis, giving the names of the factor columns:

```r
res <- multigroup_main(x, sample_info, factors = c("sex", "genotype"))
names(res)
head(res$anova_results)
```

For two groups, use `twogroup_main()`, e.g. `twogroup_main(x, sample_info, factors = "sex")`.

Additionally, statomatic includes stand-alone functions which can be used on data.frames. These functions treat each row as a separate variable and each column as an observation. Unlike the main functions, they take the group labels as a vector, which must be in the same order as the columns of `x`:

```r
#group labels in the column order of x
idx <- match(colnames(x), sample_info$sample)
sex <- sample_info$sex[idx]
genotype <- sample_info$genotype[idx]
group <- paste(sex, genotype, sep = "_")

#sort data according to statomatic's procedure, one factor or several
res1 <- test_norm(x, group)
res2 <- test_norm(x, data.frame(sex = sex, genotype = genotype))

#anova, one factor or several
anova_res1 <- run_anova(x, group)
anova_res2 <- run_anova(x, data.frame(sex = sex, genotype = genotype))

#welch and kruskal-wallis tests, one factor
welch_res <- run_welch(x, group)
kw_res <- run_kruskal(x, group)

#multiple-comparison tests: Tukey HSD, Dunnett T3 and Dunn's test, one factor
tukey_res <- run_tukey(x, group)
dunnett_res <- run_dunnett(x, group)
dunn_res <- run_dunn(x, group)

#two-group tests, one factor with exactly two groups
t_res <- run_ttest(x, sex)
wilcox_res <- run_wilcox(x, sex)

#fold changes, one factor
fc <- fold_change(x, group)
```

## Limitations

Statomatic aims to make addressing the mathematical assumptions of hypothesis tests as easy as possible for users, but be aware that not all kinds of data can be treated in the same way. In multiple examples shown here, bulk RNA-sequencing data is used to demonstrate how to use statomatic. However, there are highly specialized tools specifically designed to handle RNA-seq data, such as DESeq2, which perform many operations and checks on the data that statomatic does not. Although this package draws much inspiration from DESeq2, statomatic should not be seen as a replacement for the highly specialized procedures performed by other tools designed to handle specific data challenges.

Some limitations of the sorting approach itself:

- With few samples per group (4 to 6 is common), the Shapiro-Wilk and Levene tests have little power, so many rows are sorted as normal with equal variance even when they are not.
- Rows with unequal variances are often sorted as not normal, because unequal variances make the pooled residuals look non-normal.
- Interaction p-values are only available for rows sorted as normal with equal variance (ANOVA); the Welch and Kruskal-Wallis tests are one-way tests.
- Dunnett T3 p-values are calculated with a random (Monte Carlo) method; use `set.seed()` for reproducible results.

The approach and limitations vignette describes these in detail, with simulations of how well statomatic sorts rows.

## Changes in version 1.0.0

Version 1.0.0 is a rebuild of statomatic with a new interface: `multigroup_main()` and `twogroup_main()` take the sample table and the names of its factor columns, and the statomatic data set (`make_sds()`, `sds_analyze()`) was removed. See [NEWS.md](NEWS.md) for all changes.

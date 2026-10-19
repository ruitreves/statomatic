# Basic workflow

This vignette demonstrates the basic usage of the statomatic package to
analyze data from an experiment.

It includes three analyses, showing how statomatic handles different
experimental designs: a multi-group (\> 2) two-factor design, a
two-group one-factor design, and a multi-group one-factor design.

First, load the statomatic package.

``` r
library(statomatic)
```

## How statomatic chooses a test

Statomatic tests each row of the data (each gene, taxon, metabolite, …)
separately. For every row it first checks whether the residuals of a
linear model are normally distributed (Shapiro-Wilk test) and whether
the groups have equal variances (Levene test). This sorts the rows into
three sets:

| Set | Meaning | 2 groups | \> 2 groups | \> 2 groups: pairwise comparisons |
|----|----|----|----|----|
| ne | normal, equal variance | t-test | ANOVA | Tukey HSD |
| nu | normal, unequal variance | Welch test | Welch test | Dunnett T3 |
| nn | not normal | Wilcoxon test | Kruskal-Wallis test | Dunn test |

Each set is then analyzed with the tests suited to it, and fold changes
are calculated between every pair of groups.

## Part 1: Multi-group, two-factor design

The data analyzed here is from an RNA sequencing experiment which aimed
to identify differentially expressed genes across several experimental
conditions. The details of this experiment won’t be discussed here.

Load the data to be analyzed. Here we analyze the normalized counts, as
created by DESeq2.

``` r
norm_counts <- read.csv("normalized_counts.csv")
head(norm_counts)
#>                    X Sample_600 Sample_601 Sample_602 Sample_616 Sample_617
#> 1 ENSMUSG00000051285  169.21935  309.83231 417.709088  214.03944  202.50044
#> 2 ENSMUSG00000057363   47.53353   66.64997  54.805802   41.23696   49.30862
#> 3 ENSMUSG00000100555   23.76676   55.84187   0.000000   48.10978   18.19153
#> 4 ENSMUSG00000078184   50.38554   59.44457   5.924952   12.76382   24.41495
#> 5 ENSMUSG00000061024  273.79311  293.62015 232.554350  144.32935  360.00079
#> 6 ENSMUSG00000025911  154.95929  525.99438 309.578721  208.14845  256.59631
#>   Sample_603 Sample_604 Sample_605 Sample_606 Sample_621 Sample_607 Sample_608
#> 1  110.24045  171.38160   61.79141 190.642528  427.33305  228.73359  410.89116
#> 2   57.16171    0.00000    0.00000 103.558904   58.42021   63.76202   88.17407
#> 3   28.58086  151.21906   68.65713   5.884029   30.29196   21.25401   15.87133
#> 4   32.66384   10.08127    0.00000  44.718618   37.86495   19.22981   18.51655
#> 5  145.62627   55.44699  109.85140 222.416283  315.90190  253.02388  319.19013
#> 6  155.15322  171.38160    0.00000 264.781289  529.02750  301.60446  336.82494
#>   Sample_609 Sample_610 Sample_623 Sample_611 Sample_612 Sample_613 Sample_614
#> 1  313.05044  412.25499  133.76758  93.303195 248.040131   46.21469  417.35080
#> 2   60.75650   14.45406   99.46820 103.366522  18.445254    0.00000   65.25144
#> 3   45.30993   40.01052   30.86944   5.093783   1.456204  390.25738    0.00000
#> 4   39.13131   36.86833   56.59398  57.025521  49.996347    0.00000   82.04884
#> 5  209.04355  211.36447  221.23100 240.525945 248.040131  302.96296  440.60874
#> 6  375.86649  528.09696  265.82019 148.589375 171.832106   10.26993  323.02694
#>   Sample_615
#> 1  513.48498
#> 2   39.36603
#> 3    0.00000
#> 4   14.61838
#> 5  213.52043
#> 6  646.77697
```

Notice that the row names of norm_counts are numbers, and that the gene
ids are actually inside the data.frame. Statomatic can only analyze
numbers, so the gene ids need to be the row names, and the gene id
column needs to be removed.

Statomatic has a built-in function for this,
[`cf()`](https://ruitreves.github.io/statomatic/reference/cf.md). It
sets a specified column of a data.frame as the new row names, and then
removes that column from the data.frame. The column can be chosen with
`column_index`, which defaults to 1 (the first column).

Using [`cf()`](https://ruitreves.github.io/statomatic/reference/cf.md)
with the default `column_index = 1` is the same as reading the file with
`read.csv("normalized_counts.csv", row.names = 1)`. Either way works;
[`cf()`](https://ruitreves.github.io/statomatic/reference/cf.md) is
handy when the data is already loaded.

``` r
norm_counts <- cf(norm_counts)
head(norm_counts)
#>                    Sample_600 Sample_601 Sample_602 Sample_616 Sample_617
#> ENSMUSG00000051285  169.21935  309.83231 417.709088  214.03944  202.50044
#> ENSMUSG00000057363   47.53353   66.64997  54.805802   41.23696   49.30862
#> ENSMUSG00000100555   23.76676   55.84187   0.000000   48.10978   18.19153
#> ENSMUSG00000078184   50.38554   59.44457   5.924952   12.76382   24.41495
#> ENSMUSG00000061024  273.79311  293.62015 232.554350  144.32935  360.00079
#> ENSMUSG00000025911  154.95929  525.99438 309.578721  208.14845  256.59631
#>                    Sample_603 Sample_604 Sample_605 Sample_606 Sample_621
#> ENSMUSG00000051285  110.24045  171.38160   61.79141 190.642528  427.33305
#> ENSMUSG00000057363   57.16171    0.00000    0.00000 103.558904   58.42021
#> ENSMUSG00000100555   28.58086  151.21906   68.65713   5.884029   30.29196
#> ENSMUSG00000078184   32.66384   10.08127    0.00000  44.718618   37.86495
#> ENSMUSG00000061024  145.62627   55.44699  109.85140 222.416283  315.90190
#> ENSMUSG00000025911  155.15322  171.38160    0.00000 264.781289  529.02750
#>                    Sample_607 Sample_608 Sample_609 Sample_610 Sample_623
#> ENSMUSG00000051285  228.73359  410.89116  313.05044  412.25499  133.76758
#> ENSMUSG00000057363   63.76202   88.17407   60.75650   14.45406   99.46820
#> ENSMUSG00000100555   21.25401   15.87133   45.30993   40.01052   30.86944
#> ENSMUSG00000078184   19.22981   18.51655   39.13131   36.86833   56.59398
#> ENSMUSG00000061024  253.02388  319.19013  209.04355  211.36447  221.23100
#> ENSMUSG00000025911  301.60446  336.82494  375.86649  528.09696  265.82019
#>                    Sample_611 Sample_612 Sample_613 Sample_614 Sample_615
#> ENSMUSG00000051285  93.303195 248.040131   46.21469  417.35080  513.48498
#> ENSMUSG00000057363 103.366522  18.445254    0.00000   65.25144   39.36603
#> ENSMUSG00000100555   5.093783   1.456204  390.25738    0.00000    0.00000
#> ENSMUSG00000078184  57.025521  49.996347    0.00000   82.04884   14.61838
#> ENSMUSG00000061024 240.525945 248.040131  302.96296  440.60874  213.52043
#> ENSMUSG00000025911 148.589375 171.832106   10.26993  323.02694  646.77697
```

Now the data is ready to be analyzed.

We also need to tell statomatic which experimental conditions each
sample belongs to. We do this with a data.frame, `sample_info`, with one
row per sample: a column with the sample names (the column names of
norm_counts), and one column per experimental condition. The
experimental conditions here are sex and genotype.

``` r
sample_info <- data.frame(sample = colnames(norm_counts),
                          sex = c(rep("Male", 10), rep("Female", 10)),
                          genotype = c(rep("WT", 5), rep("KO", 5), rep("WT", 5), rep("KO", 5)))
sample_info
#>        sample    sex genotype
#> 1  Sample_600   Male       WT
#> 2  Sample_601   Male       WT
#> 3  Sample_602   Male       WT
#> 4  Sample_616   Male       WT
#> 5  Sample_617   Male       WT
#> 6  Sample_603   Male       KO
#> 7  Sample_604   Male       KO
#> 8  Sample_605   Male       KO
#> 9  Sample_606   Male       KO
#> 10 Sample_621   Male       KO
#> 11 Sample_607 Female       WT
#> 12 Sample_608 Female       WT
#> 13 Sample_609 Female       WT
#> 14 Sample_610 Female       WT
#> 15 Sample_623 Female       WT
#> 16 Sample_611 Female       KO
#> 17 Sample_612 Female       KO
#> 18 Sample_613 Female       KO
#> 19 Sample_614 Female       KO
#> 20 Sample_615 Female       KO
```

Statomatic matches the rows of `sample_info` to the columns of the data
by sample name, so the order of the rows doesn’t matter. It finds the
column with the sample names automatically (here, `sample`). If your
sample names are the row names of `sample_info` instead, that works too.

Now we run the analysis with
[`multigroup_main()`](https://ruitreves.github.io/statomatic/reference/multigroup_main.md).
We give it the data, `sample_info`, and the names of the columns of
`sample_info` to use as the experimental factors.

``` r
set.seed(1)
res <- multigroup_main(norm_counts, sample_info, factors = c("sex", "genotype"))
#> Analyzing 100 rows in 4 groups: Female_KO, Female_WT, Male_KO, Male_WT
#> Sorting rows by normality and equality of variance...
#>   ne: 62  nu: 1  nn: 37
#> Calculating fold changes...
#> Running ANOVA + Tukey on 62 rows...
#> Running Welch + Dunnett T3 on 1 row...
#> Running Kruskal-Wallis + Dunn on 37 rows...
#> Done.
```

With two factors, statomatic uses a two-way ANOVA for the normal,
equal-variance (ne) rows. The pairwise comparisons and fold changes
compare every combination of the factors: statomatic combines sex and
genotype into four groups (Female_KO, Female_WT, Male_KO and Male_WT).

The Dunnett T3 test calculates its p-values with a random (Monte Carlo)
method, so its p-values can change very slightly from one run to the
next. Setting a seed with
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the analysis
makes the results reproducible.

The result is a list:

``` r
names(res)
#> [1] "anova_results" "welch_results" "kw_results"    "foldchanges"  
#> [5] "ne"            "nu"            "nn"            "tests"
```

- `anova_results`, `welch_results` and `kw_results` hold the results for
  the ne, nu and nn rows.
- `foldchanges` holds the fold changes and group means of every row.
- `ne`, `nu` and `nn` are the rows of the data in each set.
- `tests` shows, for every row, the Shapiro-Wilk and Levene p-values and
  the set it was sorted into.

Let’s look at the ANOVA results.

``` r
head(res$anova_results)
#>                    pvalue_sex  padj_sex pvalue_genotype padj_genotype
#> ENSMUSG00000051285  0.4181785 0.9384844       0.4267951     0.7520794
#> ENSMUSG00000057363  0.6363190 0.9866797       0.3757298     0.7059165
#> ENSMUSG00000078184  0.3813036 0.9384844       0.9576038     0.9895239
#> ENSMUSG00000061024  0.1879463 0.9384844       0.5524404     0.7611401
#> ENSMUSG00000025911  0.5146527 0.9384844       0.3078191     0.6998790
#> ENSMUSG00000057335  0.9928111 0.9935637       0.1615905     0.5288119
#>                    pvalue_sex:genotype padj_sex:genotype
#> ENSMUSG00000051285          0.79596064         0.9913574
#> ENSMUSG00000057363          0.70346317         0.9913574
#> ENSMUSG00000078184          0.57458343         0.9913574
#> ENSMUSG00000061024          0.08030823         0.9913574
#> ENSMUSG00000025911          0.83172345         0.9913574
#> ENSMUSG00000057335          0.85432131         0.9913574
#>                    pvalue_Female_KO_vs_Female_WT padj_Female_KO_vs_Female_WT
#> ENSMUSG00000051285                     0.9790520                   0.9998409
#> ENSMUSG00000057363                     0.7955281                   0.9998409
#> ENSMUSG00000078184                     0.9699575                   0.9998409
#> ENSMUSG00000061024                     0.8094858                   0.9998409
#> ENSMUSG00000025911                     0.8062227                   0.9998409
#> ENSMUSG00000057335                     0.8020005                   0.9998409
#>                    log2FoldChange_Female_KO_vs_Female_WT
#> ENSMUSG00000051285                            -0.1849281
#> ENSMUSG00000057363                            -0.5285301
#> ENSMUSG00000078184                             0.2579516
#> ENSMUSG00000061024                             0.2521327
#> ENSMUSG00000025911                            -0.4755035
#> ENSMUSG00000057335                            -0.5990134
#>                    pvalue_Female_KO_vs_Male_KO padj_Female_KO_vs_Male_KO
#> ENSMUSG00000051285                   0.8652082                 0.9999797
#> ENSMUSG00000057363                   0.9998882                 0.9999797
#> ENSMUSG00000078184                   0.7280883                 0.9999797
#> ENSMUSG00000061024                   0.1413002                 0.9999797
#> ENSMUSG00000025911                   0.9883949                 0.9999797
#> ENSMUSG00000057335                   0.9992630                 0.9999797
#>                    log2FoldChange_Female_KO_vs_Male_KO
#> ENSMUSG00000051285                           0.4555891
#> ENSMUSG00000057363                           0.0472020
#> ENSMUSG00000078184                           0.7006521
#> ENSMUSG00000061024                           0.7674774
#> ENSMUSG00000025911                           0.2151199
#> ENSMUSG00000057335                           0.1067049
#>                    pvalue_Female_KO_vs_Male_WT padj_Female_KO_vs_Male_WT
#> ENSMUSG00000051285                   0.9999995                         1
#> ENSMUSG00000057363                   0.9899203                         1
#> ENSMUSG00000078184                   0.9050970                         1
#> ENSMUSG00000061024                   0.9469754                         1
#> ENSMUSG00000025911                   0.9925533                         1
#> ENSMUSG00000057335                   0.7267167                         1
#>                    log2FoldChange_Female_KO_vs_Male_WT
#> ENSMUSG00000051285                         0.005584147
#> ENSMUSG00000057363                        -0.196868100
#> ENSMUSG00000078184                         0.413461170
#> ENSMUSG00000061024                         0.148453255
#> ENSMUSG00000025911                        -0.162232725
#> ENSMUSG00000057335                        -0.672042886
#>                    pvalue_Female_WT_vs_Male_KO padj_Female_WT_vs_Male_KO
#> ENSMUSG00000051285                   0.6568296                 0.9999904
#> ENSMUSG00000057363                   0.7597698                 0.9999904
#> ENSMUSG00000078184                   0.9310571                 0.9999904
#> ENSMUSG00000061024                   0.5161829                 0.9999904
#> ENSMUSG00000025911                   0.6260176                 0.9999904
#> ENSMUSG00000057335                   0.7340525                 0.9999904
#>                    log2FoldChange_Female_WT_vs_Male_KO
#> ENSMUSG00000051285                           0.6405172
#> ENSMUSG00000057363                           0.5757321
#> ENSMUSG00000078184                           0.4427005
#> ENSMUSG00000061024                           0.5153447
#> ENSMUSG00000025911                           0.6906234
#> ENSMUSG00000057335                           0.7057183
#>                    pvalue_Female_WT_vs_Male_WT padj_Female_WT_vs_Male_WT
#> ENSMUSG00000051285                   0.9773126                 0.9999998
#> ENSMUSG00000057363                   0.9258447                 0.9999998
#> ENSMUSG00000078184                   0.9954529                 0.9999998
#> ENSMUSG00000061024                   0.9850262                 0.9999998
#> ENSMUSG00000025911                   0.9228793                 0.9999998
#> ENSMUSG00000057335                   0.9990121                 0.9999998
#>                    log2FoldChange_Female_WT_vs_Male_WT
#> ENSMUSG00000051285                          0.19051226
#> ENSMUSG00000057363                          0.33166197
#> ENSMUSG00000078184                          0.15550958
#> ENSMUSG00000061024                         -0.10367949
#> ENSMUSG00000025911                          0.31327074
#> ENSMUSG00000057335                         -0.07302953
#>                    pvalue_Male_KO_vs_Male_WT padj_Male_KO_vs_Male_WT
#> ENSMUSG00000051285                 0.8699999               0.9999869
#> ENSMUSG00000057363                 0.9820468               0.9999869
#> ENSMUSG00000078184                 0.9824965               0.9999869
#> ENSMUSG00000061024                 0.3321566               0.9999869
#> ENSMUSG00000025911                 0.9330056               0.9999869
#> ENSMUSG00000057335                 0.6535879               0.9999869
#>                    log2FoldChange_Male_KO_vs_Male_WT Female_KO_mean
#> ENSMUSG00000051285                        -0.4500049      263.67876
#> ENSMUSG00000057363                        -0.2440701       45.28585
#> ENSMUSG00000078184                        -0.2871910       40.73782
#> ENSMUSG00000061024                        -0.6190242      289.13164
#> ENSMUSG00000025911                        -0.3773527      260.09906
#> ENSMUSG00000057335                        -0.7787478      189.02091
#>                    Female_WT_mean Male_KO_mean Male_WT_mean       test
#> ENSMUSG00000051285      299.73955    192.27781    262.66013 tukey-test
#> ENSMUSG00000057363       65.32297     43.82817     51.90698 tukey-test
#> ENSMUSG00000078184       34.06800     25.06574     30.58677 tukey-test
#> ENSMUSG00000061024      242.77061    169.84857    260.85955 tukey-test
#> ENSMUSG00000025911      361.64261    224.06872    291.05543 tukey-test
#> ENSMUSG00000057335      286.30626    175.54500    301.17223 tukey-test
```

The table starts with the ANOVA p-values for each main effect and their
interaction. Then, for each pair of groups, it has the pairwise (here
Tukey) p-value, the adjusted p-value and the log2 fold change. It ends
with the mean of each group and the name of the pairwise test.

Two kinds of p-value adjustment are involved:

- The pairwise p-values (Tukey, Dunnett T3, Dunn) are already adjusted
  for comparing several pairs of groups within one row.
- The `padj` columns adjust for testing many rows, using the
  Benjamini-Hochberg method. They are calculated separately within each
  set. They can be turned off with `padj = FALSE`.

The Welch and Kruskal-Wallis results have the same layout, except that
they start with a single `pvalue` and `padj`, because those tests
compare all the combined groups at once.

``` r
head(res$welch_results)
#>                       pvalue      padj pvalue_Female_KO_vs_Female_WT
#> ENSMUSG00000004364 0.3519013 0.3519013                      0.939041
#>                    padj_Female_KO_vs_Female_WT
#> ENSMUSG00000004364                    0.939041
#>                    log2FoldChange_Female_KO_vs_Female_WT
#> ENSMUSG00000004364                            -0.3566521
#>                    pvalue_Female_KO_vs_Male_KO padj_Female_KO_vs_Male_KO
#> ENSMUSG00000004364                   0.9999963                 0.9999963
#>                    log2FoldChange_Female_KO_vs_Male_KO
#> ENSMUSG00000004364                          0.07124816
#>                    pvalue_Female_KO_vs_Male_WT padj_Female_KO_vs_Male_WT
#> ENSMUSG00000004364                   0.9364206                 0.9364206
#>                    log2FoldChange_Female_KO_vs_Male_WT
#> ENSMUSG00000004364                          -0.3414172
#>                    pvalue_Female_WT_vs_Male_KO padj_Female_WT_vs_Male_KO
#> ENSMUSG00000004364                   0.4906846                 0.4906846
#>                    log2FoldChange_Female_WT_vs_Male_KO
#> ENSMUSG00000004364                           0.4279003
#>                    pvalue_Female_WT_vs_Male_WT padj_Female_WT_vs_Male_WT
#> ENSMUSG00000004364                   0.9999996                 0.9999996
#>                    log2FoldChange_Female_WT_vs_Male_WT
#> ENSMUSG00000004364                          0.01523491
#>                    pvalue_Male_KO_vs_Male_WT padj_Male_KO_vs_Male_WT
#> ENSMUSG00000004364                 0.4026044               0.4026044
#>                    log2FoldChange_Male_KO_vs_Male_WT Female_KO_mean
#> ENSMUSG00000004364                        -0.4126654       466.9804
#>                    Female_WT_mean Male_KO_mean Male_WT_mean         test
#> ENSMUSG00000004364       597.9456     444.4786     591.6645 dunnett-test
```

To see how each row was sorted:

``` r
head(res$tests)
#>                         shapiro    levene set
#> ENSMUSG00000051285 5.955714e-01 0.4443846  ne
#> ENSMUSG00000057363 3.464768e-01 0.2722540  ne
#> ENSMUSG00000100555 3.976942e-05 0.6284235  nn
#> ENSMUSG00000078184 7.456867e-01 0.5382225  ne
#> ENSMUSG00000061024 3.367137e-01 0.7158967  ne
#> ENSMUSG00000025911 9.854698e-02 0.7137068  ne
table(res$tests$set)
#> 
#> ne nu nn 
#> 62  1 37
```

All the pairwise results can be combined into one table, because they
have the same columns in all three sets:

``` r
pairwise <- lapply(res[c("anova_results", "welch_results", "kw_results")], function(tab) {
    tab[, grep("_vs_|_mean$|^test$", colnames(tab))]
})
all_results <- do.call(rbind, unname(pairwise))
head(all_results)
#>                    pvalue_Female_KO_vs_Female_WT padj_Female_KO_vs_Female_WT
#> ENSMUSG00000051285                     0.9790520                   0.9998409
#> ENSMUSG00000057363                     0.7955281                   0.9998409
#> ENSMUSG00000078184                     0.9699575                   0.9998409
#> ENSMUSG00000061024                     0.8094858                   0.9998409
#> ENSMUSG00000025911                     0.8062227                   0.9998409
#> ENSMUSG00000057335                     0.8020005                   0.9998409
#>                    log2FoldChange_Female_KO_vs_Female_WT
#> ENSMUSG00000051285                            -0.1849281
#> ENSMUSG00000057363                            -0.5285301
#> ENSMUSG00000078184                             0.2579516
#> ENSMUSG00000061024                             0.2521327
#> ENSMUSG00000025911                            -0.4755035
#> ENSMUSG00000057335                            -0.5990134
#>                    pvalue_Female_KO_vs_Male_KO padj_Female_KO_vs_Male_KO
#> ENSMUSG00000051285                   0.8652082                 0.9999797
#> ENSMUSG00000057363                   0.9998882                 0.9999797
#> ENSMUSG00000078184                   0.7280883                 0.9999797
#> ENSMUSG00000061024                   0.1413002                 0.9999797
#> ENSMUSG00000025911                   0.9883949                 0.9999797
#> ENSMUSG00000057335                   0.9992630                 0.9999797
#>                    log2FoldChange_Female_KO_vs_Male_KO
#> ENSMUSG00000051285                           0.4555891
#> ENSMUSG00000057363                           0.0472020
#> ENSMUSG00000078184                           0.7006521
#> ENSMUSG00000061024                           0.7674774
#> ENSMUSG00000025911                           0.2151199
#> ENSMUSG00000057335                           0.1067049
#>                    pvalue_Female_KO_vs_Male_WT padj_Female_KO_vs_Male_WT
#> ENSMUSG00000051285                   0.9999995                         1
#> ENSMUSG00000057363                   0.9899203                         1
#> ENSMUSG00000078184                   0.9050970                         1
#> ENSMUSG00000061024                   0.9469754                         1
#> ENSMUSG00000025911                   0.9925533                         1
#> ENSMUSG00000057335                   0.7267167                         1
#>                    log2FoldChange_Female_KO_vs_Male_WT
#> ENSMUSG00000051285                         0.005584147
#> ENSMUSG00000057363                        -0.196868100
#> ENSMUSG00000078184                         0.413461170
#> ENSMUSG00000061024                         0.148453255
#> ENSMUSG00000025911                        -0.162232725
#> ENSMUSG00000057335                        -0.672042886
#>                    pvalue_Female_WT_vs_Male_KO padj_Female_WT_vs_Male_KO
#> ENSMUSG00000051285                   0.6568296                 0.9999904
#> ENSMUSG00000057363                   0.7597698                 0.9999904
#> ENSMUSG00000078184                   0.9310571                 0.9999904
#> ENSMUSG00000061024                   0.5161829                 0.9999904
#> ENSMUSG00000025911                   0.6260176                 0.9999904
#> ENSMUSG00000057335                   0.7340525                 0.9999904
#>                    log2FoldChange_Female_WT_vs_Male_KO
#> ENSMUSG00000051285                           0.6405172
#> ENSMUSG00000057363                           0.5757321
#> ENSMUSG00000078184                           0.4427005
#> ENSMUSG00000061024                           0.5153447
#> ENSMUSG00000025911                           0.6906234
#> ENSMUSG00000057335                           0.7057183
#>                    pvalue_Female_WT_vs_Male_WT padj_Female_WT_vs_Male_WT
#> ENSMUSG00000051285                   0.9773126                 0.9999998
#> ENSMUSG00000057363                   0.9258447                 0.9999998
#> ENSMUSG00000078184                   0.9954529                 0.9999998
#> ENSMUSG00000061024                   0.9850262                 0.9999998
#> ENSMUSG00000025911                   0.9228793                 0.9999998
#> ENSMUSG00000057335                   0.9990121                 0.9999998
#>                    log2FoldChange_Female_WT_vs_Male_WT
#> ENSMUSG00000051285                          0.19051226
#> ENSMUSG00000057363                          0.33166197
#> ENSMUSG00000078184                          0.15550958
#> ENSMUSG00000061024                         -0.10367949
#> ENSMUSG00000025911                          0.31327074
#> ENSMUSG00000057335                         -0.07302953
#>                    pvalue_Male_KO_vs_Male_WT padj_Male_KO_vs_Male_WT
#> ENSMUSG00000051285                 0.8699999               0.9999869
#> ENSMUSG00000057363                 0.9820468               0.9999869
#> ENSMUSG00000078184                 0.9824965               0.9999869
#> ENSMUSG00000061024                 0.3321566               0.9999869
#> ENSMUSG00000025911                 0.9330056               0.9999869
#> ENSMUSG00000057335                 0.6535879               0.9999869
#>                    log2FoldChange_Male_KO_vs_Male_WT Female_KO_mean
#> ENSMUSG00000051285                        -0.4500049      263.67876
#> ENSMUSG00000057363                        -0.2440701       45.28585
#> ENSMUSG00000078184                        -0.2871910       40.73782
#> ENSMUSG00000061024                        -0.6190242      289.13164
#> ENSMUSG00000025911                        -0.3773527      260.09906
#> ENSMUSG00000057335                        -0.7787478      189.02091
#>                    Female_WT_mean Male_KO_mean Male_WT_mean       test
#> ENSMUSG00000051285      299.73955    192.27781    262.66013 tukey-test
#> ENSMUSG00000057363       65.32297     43.82817     51.90698 tukey-test
#> ENSMUSG00000078184       34.06800     25.06574     30.58677 tukey-test
#> ENSMUSG00000061024      242.77061    169.84857    260.85955 tukey-test
#> ENSMUSG00000025911      361.64261    224.06872    291.05543 tukey-test
#> ENSMUSG00000057335      286.30626    175.54500    301.17223 tukey-test
```

## Part 2: Two-group, one-factor design

First, load the data. For this example we’ll look at gene expression
data from an RNA-seq experiment, normalized by DESeq2. There are two
groups in this data set: samples 1-4 and samples 5-8. For the purposes
of this example, we’ll call them “treated” and “untreated”.

``` r
x <- read.csv("two_group_normalized_counts.csv")
x <- cf(x)
head(x)
#>                      sample_1   sample_2   sample_3   sample_4   sample_5
#> ENSMUSG00000000001 1707.89280 1843.60723 1687.02686 1501.02388 1806.65191
#> ENSMUSG00000000028   27.84608   30.17596   29.61100   69.66777   48.57703
#> ENSMUSG00000000056  468.74232  517.30211  556.20667  404.54810  455.79719
#> ENSMUSG00000000058   60.33317   81.90617   75.22795   51.45915  116.79157
#> ENSMUSG00000000078 2706.87086 3110.99739 2504.93077 2310.91176 2468.12630
#> ENSMUSG00000000085  158.95470  145.13198  226.48416  192.37806  248.05289
#>                      sample_6   sample_7   sample_8
#> ENSMUSG00000000001 1916.63101 1984.70503 2100.84931
#> ENSMUSG00000000028   51.32071   42.77382   52.52123
#> ENSMUSG00000000056  468.79492  497.03174  504.00183
#> ENSMUSG00000000058   88.82430   71.86001   68.68161
#> ENSMUSG00000000078 2762.43574 2557.01869 2697.77332
#> ENSMUSG00000000085  224.03462  203.60336  198.97467
```

Now we make the `sample_info` to map our samples to their group.

``` r
sample_info <- data.frame(sample = colnames(x),
                          group = c(rep("treated", 4), rep("untreated", 4)))
sample_info
#>     sample     group
#> 1 sample_1   treated
#> 2 sample_2   treated
#> 3 sample_3   treated
#> 4 sample_4   treated
#> 5 sample_5 untreated
#> 6 sample_6 untreated
#> 7 sample_7 untreated
#> 8 sample_8 untreated
```

For two groups we use
[`twogroup_main()`](https://ruitreves.github.io/statomatic/reference/twogroup_main.md),
with the name of the column that holds the groups.

``` r
res <- twogroup_main(x, sample_info, factors = "group")
#> Analyzing 100 rows in 2 groups: treated, untreated
#> Sorting rows by normality and equality of variance...
#>   ne: 91  nu: 7  nn: 2
#> Calculating fold changes...
#> Running t-test on 91 rows...
#> Running Welch test on 7 rows...
#> Running Wilcoxon test on 2 rows...
#> Done.
names(res)
#> [1] "ttest_results"  "welch_results"  "wilcox_results" "foldchanges"   
#> [5] "ne"             "nu"             "nn"             "tests"
```

The results are in `ttest_results`, `welch_results` and
`wilcox_results`. Each has the p-value, the adjusted p-value, the log2
fold change, the two group means and the name of the test:

``` r
head(res$ttest_results)
#>                        pvalue      padj log2FoldChange_treated_vs_untreated
#> ENSMUSG00000000001 0.02895439 0.5936437                         -0.21245532
#> ENSMUSG00000000056 0.88359469 0.9950100                          0.01577680
#> ENSMUSG00000000058 0.18825794 0.6814105                         -0.36421697
#> ENSMUSG00000000078 0.84664860 0.9950100                          0.02026956
#> ENSMUSG00000000085 0.12619557 0.5936437                         -0.27483781
#> ENSMUSG00000000088 0.82980304 0.9950100                         -0.02965276
#>                    treated_mean untreated_mean   test
#> ENSMUSG00000000001   1684.88769     1952.20932 t-test
#> ENSMUSG00000000056    486.69980      481.40642 t-test
#> ENSMUSG00000000058     67.23161       86.53937 t-test
#> ENSMUSG00000000078   2658.42770     2621.33851 t-test
#> ENSMUSG00000000085    180.73722      218.66639 t-test
#> ENSMUSG00000000088    593.82628      606.15792 t-test
```

Because all three tables have the same columns, they can be combined
directly:

``` r
all_results <- rbind(res$ttest_results, res$welch_results, res$wilcox_results)
head(all_results)
#>                        pvalue      padj log2FoldChange_treated_vs_untreated
#> ENSMUSG00000000001 0.02895439 0.5936437                         -0.21245532
#> ENSMUSG00000000056 0.88359469 0.9950100                          0.01577680
#> ENSMUSG00000000058 0.18825794 0.6814105                         -0.36421697
#> ENSMUSG00000000078 0.84664860 0.9950100                          0.02026956
#> ENSMUSG00000000085 0.12619557 0.5936437                         -0.27483781
#> ENSMUSG00000000088 0.82980304 0.9950100                         -0.02965276
#>                    treated_mean untreated_mean   test
#> ENSMUSG00000000001   1684.88769     1952.20932 t-test
#> ENSMUSG00000000056    486.69980      481.40642 t-test
#> ENSMUSG00000000058     67.23161       86.53937 t-test
#> ENSMUSG00000000078   2658.42770     2621.33851 t-test
#> ENSMUSG00000000085    180.73722      218.66639 t-test
#> ENSMUSG00000000088    593.82628      606.15792 t-test
```

## Part 3: Multi-group, one-factor design

The data used here is species-level data from a metagenomic analysis.
There are 24 samples in three groups, with 8 samples per group. We’ll
call the groups group1, group2 and group3. This analysis moves a bit
faster, as the steps are the same as before.

``` r
x <- read.csv("species_abundance.csv")
x <- cf(x)
head(x)
#>                                  sample_1 sample_2 sample_3 sample_4 sample_5
#> Lachnospiraceae bacterium              14       33      245      224      167
#> Dubosiella newyorkensis               583      181      118      145      112
#> Parasutterella excrementihominis        4      241      431      327      343
#> Barnesiella sp. CU968                 261      146      168      138      128
#> Bacteroides thetaiotaomicron          356       74       79       91       95
#> Lachnospiraceae bacterium 28-4          2        2     1176        7       92
#>                                  sample_6 sample_7 sample_8 sample_9 sample_10
#> Lachnospiraceae bacterium               9      996       71      312        85
#> Dubosiella newyorkensis               118      665      131      121       117
#> Parasutterella excrementihominis       21      128      241      367       364
#> Barnesiella sp. CU968                 175      251      144      129       129
#> Bacteroides thetaiotaomicron          155       80      129       76        83
#> Lachnospiraceae bacterium 28-4          5     1170       43       19         4
#>                                  sample_11 sample_12 sample_13 sample_14
#> Lachnospiraceae bacterium               91      1731        22        12
#> Dubosiella newyorkensis                113       210       120       460
#> Parasutterella excrementihominis       246       255       393         8
#> Barnesiella sp. CU968                  393       130        25        15
#> Bacteroides thetaiotaomicron            75        82        77       371
#> Lachnospiraceae bacterium 28-4          42       272         0         3
#>                                  sample_15 sample_16 sample_17 sample_18
#> Lachnospiraceae bacterium              221       289        11        35
#> Dubosiella newyorkensis                117       105       114       737
#> Parasutterella excrementihominis       386       120       385         7
#> Barnesiella sp. CU968                  417       191       138       132
#> Bacteroides thetaiotaomicron           119        91        78        76
#> Lachnospiraceae bacterium 28-4         671        25         3         3
#>                                  sample_19 sample_20 sample_21 sample_22
#> Lachnospiraceae bacterium              168      1907        42       103
#> Dubosiella newyorkensis                728       123       245       220
#> Parasutterella excrementihominis         4       433         4       117
#> Barnesiella sp. CU968                  123       116       213       284
#> Bacteroides thetaiotaomicron            75       116       513       582
#> Lachnospiraceae bacterium 28-4          14       461         3         5
#>                                  sample_23 sample_24
#> Lachnospiraceae bacterium              136        28
#> Dubosiella newyorkensis                402        57
#> Parasutterella excrementihominis        10         1
#> Barnesiella sp. CU968                  123       420
#> Bacteroides thetaiotaomicron            81       721
#> Lachnospiraceae bacterium 28-4           8         4
sample_info <- data.frame(sample = colnames(x),
                          group = c(rep("group1", 8), rep("group2", 8), rep("group3", 8)))
```

``` r
set.seed(1)
res <- multigroup_main(x, sample_info, factors = "group")
#> Analyzing 100 rows in 3 groups: group1, group2, group3
#> Sorting rows by normality and equality of variance...
#>   ne: 18  nu: 3  nn: 79
#> Calculating fold changes...
#> Running ANOVA + Tukey on 18 rows...
#> Running Welch + Dunnett T3 on 3 rows...
#> Running Kruskal-Wallis + Dunn on 79 rows...
#> Done.
```

With a single factor, statomatic uses a one-way ANOVA for the ne rows,
and the ANOVA results start with a single `pvalue` and `padj`, like the
Welch and Kruskal-Wallis results.

``` r
head(res$anova_results)
#>                                        pvalue        padj
#> Parasutterella excrementihominis 0.1991824896 0.298773734
#> Muricaecibacterium torontonense  0.0681265365 0.122627766
#> Muribaculaceae bacterium         0.0010328865 0.006281909
#> Sangeribacter muris              0.3985184169 0.448333219
#> Anaerotruncus colihominis        0.6083336357 0.644117967
#> Bacteroides sp. CAG:927          0.0009232984 0.006281909
#>                                  pvalue_group1_vs_group2 padj_group1_vs_group2
#> Parasutterella excrementihominis               0.8062039                     1
#> Muricaecibacterium torontonense                0.5464864                     1
#> Muribaculaceae bacterium                       0.4923547                     1
#> Sangeribacter muris                            0.8117231                     1
#> Anaerotruncus colihominis                      0.7787006                     1
#> Bacteroides sp. CAG:927                        1.0000000                     1
#>                                  log2FoldChange_group1_vs_group2
#> Parasutterella excrementihominis                      -0.3011695
#> Muricaecibacterium torontonense                       -0.5635121
#> Muribaculaceae bacterium                               0.6818240
#> Sangeribacter muris                                    0.2115041
#> Anaerotruncus colihominis                              0.4340958
#> Bacteroides sp. CAG:927                                0.0000000
#>                                  pvalue_group1_vs_group3 padj_group1_vs_group3
#> Parasutterella excrementihominis             0.461211401            0.63860040
#> Muricaecibacterium torontonense              0.357348922            0.53602338
#> Muribaculaceae bacterium                     0.014194068            0.06145687
#> Sangeribacter muris                          0.727643346            0.87317201
#> Anaerotruncus colihominis                    0.949127148            0.99734110
#> Bacteroides sp. CAG:927                      0.002511996            0.04521594
#>                                  log2FoldChange_group1_vs_group3
#> Parasutterella excrementihominis                       0.8531586
#> Muricaecibacterium torontonense                        1.4364879
#> Muribaculaceae bacterium                              -1.0093379
#> Sangeribacter muris                                   -0.2250666
#> Anaerotruncus colihominis                             -0.1610624
#> Bacteroides sp. CAG:927                               -1.6455040
#>                                  pvalue_group2_vs_group3 padj_group2_vs_group3
#> Parasutterella excrementihominis             0.182048409           0.297897397
#> Muricaecibacterium torontonense              0.055825056           0.100485100
#> Muribaculaceae bacterium                     0.000975412           0.007926544
#> Sangeribacter muris                          0.366905054           0.440286065
#> Anaerotruncus colihominis                    0.593457381           0.628366638
#> Bacteroides sp. CAG:927                      0.002511996           0.009043187
#>                                  log2FoldChange_group2_vs_group3 group1_mean
#> Parasutterella excrementihominis                       1.1543281     217.000
#> Muricaecibacterium torontonense                        2.0000000     163.750
#> Muribaculaceae bacterium                              -1.6911619      38.500
#> Sangeribacter muris                                   -0.4365707      19.250
#> Anaerotruncus colihominis                             -0.5951583      15.875
#> Bacteroides sp. CAG:927                               -1.6455040       8.750
#>                                  group2_mean group3_mean       test
#> Parasutterella excrementihominis     267.375     120.125 tukey-test
#> Muricaecibacterium torontonense      242.000      60.500 tukey-test
#> Muribaculaceae bacterium              24.000      77.500 tukey-test
#> Sangeribacter muris                   16.625      22.500 tukey-test
#> Anaerotruncus colihominis             11.750      17.750 tukey-test
#> Bacteroides sp. CAG:927                8.750      27.375 tukey-test
```

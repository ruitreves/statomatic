# Additional usage of statomatic

This vignette demonstrates additional usage of the statomatic package
beyond that shown in the basic workflow vignette. It is recommended to
read the basic workflow vignette before reading this one.

### Begin by loading the statomatic package

``` r
library(statomatic)
```

### A general overview

[`multigroup_main()`](https://ruitreves.github.io/statomatic/reference/multigroup_main.md)
and
[`twogroup_main()`](https://ruitreves.github.io/statomatic/reference/twogroup_main.md)
run the whole statomatic analysis. Each of their steps is also available
as a separate function, which can be used on its own. This section shows
how.

We’ll use the same normalized counts that were used in the multi-group,
two-factor example of the basic workflow vignette.

``` r
x <- read.csv("normalized_counts.csv")
x <- cf(x)
head(x)
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

As before,
[`cf()`](https://ruitreves.github.io/statomatic/reference/cf.md) moves
the gene ids from the first column into the row names, so the data.frame
contains only the numbers we want to analyze. Reading the file with
`read.csv("normalized_counts.csv", row.names = 1)` does the same.

We also need a `sample_info` data.frame that maps our samples (the
columns of `x`) to the experimental factors.

``` r
sample_info <- data.frame(sample = colnames(x),
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

### The steps inside multigroup_main

`multigroup_main(x, sample_info, factors = c("sex", "genotype"))` does
the following:

1.  Matches the rows of `sample_info` to the columns of `x` by sample
    name, checks the design, and combines the factors into groups
    (Female_KO, Female_WT, Male_KO and Male_WT).
2.  Sorts the rows of `x` into the ne, nu and nn sets with
    [`test_norm()`](https://ruitreves.github.io/statomatic/reference/test_norm.md).
3.  Calculates fold changes between every pair of groups with
    [`fold_change()`](https://ruitreves.github.io/statomatic/reference/fold_change.md).
4.  Tests each set:
    - ne:
      [`run_anova()`](https://ruitreves.github.io/statomatic/reference/run_anova.md)
      with all factors, then
      [`run_tukey()`](https://ruitreves.github.io/statomatic/reference/run_tukey.md)
      on the groups
    - nu:
      [`run_welch()`](https://ruitreves.github.io/statomatic/reference/run_welch.md)
      and
      [`run_dunnett()`](https://ruitreves.github.io/statomatic/reference/run_dunnett.md)
      on the groups
    - nn:
      [`run_kruskal()`](https://ruitreves.github.io/statomatic/reference/run_kruskal.md)
      and
      [`run_dunn()`](https://ruitreves.github.io/statomatic/reference/run_dunn.md)
      on the groups
5.  Combines each set’s results with the fold changes.

[`twogroup_main()`](https://ruitreves.github.io/statomatic/reference/twogroup_main.md)
follows the same steps, with
[`run_ttest()`](https://ruitreves.github.io/statomatic/reference/run_ttest.md),
[`run_welch()`](https://ruitreves.github.io/statomatic/reference/run_welch.md)
and
[`run_wilcox()`](https://ruitreves.github.io/statomatic/reference/run_wilcox.md).

### Using the functions individually

There is one important difference between the main functions and the
individual functions. The main functions match samples to their groups
**by name**. The individual functions don’t: they take the group labels
as a plain vector, and **the labels must be in the same order as the
columns of `x`**. If they are not, the results will be wrong without any
error.

A safe way to get the labels in the right order is to look them up by
sample name:

``` r
idx <- match(colnames(x), sample_info$sample)
sex <- sample_info$sex[idx]
genotype <- sample_info$genotype[idx]
```

For the pairwise tests and fold changes, the groups are the combinations
of the factors, joined with “\_“, as
[`multigroup_main()`](https://ruitreves.github.io/statomatic/reference/multigroup_main.md)
does:

``` r
group <- paste(sex, genotype, sep = "_")
group
#>  [1] "Male_WT"   "Male_WT"   "Male_WT"   "Male_WT"   "Male_WT"   "Male_KO"  
#>  [7] "Male_KO"   "Male_KO"   "Male_KO"   "Male_KO"   "Female_WT" "Female_WT"
#> [13] "Female_WT" "Female_WT" "Female_WT" "Female_KO" "Female_KO" "Female_KO"
#> [19] "Female_KO" "Female_KO"
```

[`test_norm()`](https://ruitreves.github.io/statomatic/reference/test_norm.md)
sorts the rows of the data by the distribution of their residuals and
the equality of variances between groups. With two or more factors, give
it a data.frame with one column per factor:

``` r
sorted <- test_norm(x, data.frame(sex = sex, genotype = genotype))
names(sorted)
#> [1] "ne"    "nu"    "nn"    "tests"
```

The result is a list. `ne` has the rows with normally distributed
residuals and equal variances, `nu` the rows with normally distributed
residuals but unequal variances, and `nn` the rows with residuals that
are not normally distributed. `tests` shows the Shapiro-Wilk and Levene
p-values of every row and the set it was sorted into:

``` r
head(sorted$tests)
#>                         shapiro    levene set
#> ENSMUSG00000051285 5.955714e-01 0.4443846  ne
#> ENSMUSG00000057363 3.464768e-01 0.2722540  ne
#> ENSMUSG00000100555 3.976942e-05 0.6284235  nn
#> ENSMUSG00000078184 7.456867e-01 0.5382225  ne
#> ENSMUSG00000061024 3.367137e-01 0.7158967  ne
#> ENSMUSG00000025911 9.854698e-02 0.7137068  ne
```

To use
[`test_norm()`](https://ruitreves.github.io/statomatic/reference/test_norm.md)
with only one factor, give it a vector instead:

``` r
sorted_one_factor <- test_norm(x, group)
table(sorted_one_factor$tests$set)
#> 
#> ne nu nn 
#> 62  1 37
```

The usage of the rest of the functions is similar, and their names
describe what they do:
[`run_anova()`](https://ruitreves.github.io/statomatic/reference/run_anova.md)
runs an ANOVA,
[`run_tukey()`](https://ruitreves.github.io/statomatic/reference/run_tukey.md)
runs a Tukey test, and so on.

[`run_anova()`](https://ruitreves.github.io/statomatic/reference/run_anova.md)
takes one factor (a vector) or several (a data.frame), and runs a
one-way, two-way or three-way ANOVA accordingly. It uses Type II sums of
squares, so the order of the factors doesn’t change the results.

``` r
anova_res <- run_anova(sorted$ne, data.frame(sex = sex, genotype = genotype))
head(anova_res)
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
```

The other test functions take a single grouping:
[`run_tukey()`](https://ruitreves.github.io/statomatic/reference/run_tukey.md),
[`run_welch()`](https://ruitreves.github.io/statomatic/reference/run_welch.md),
[`run_dunnett()`](https://ruitreves.github.io/statomatic/reference/run_dunnett.md),
[`run_kruskal()`](https://ruitreves.github.io/statomatic/reference/run_kruskal.md)
and
[`run_dunn()`](https://ruitreves.github.io/statomatic/reference/run_dunn.md)
for more than two groups, and
[`run_ttest()`](https://ruitreves.github.io/statomatic/reference/run_ttest.md),
[`run_welch()`](https://ruitreves.github.io/statomatic/reference/run_welch.md)
and
[`run_wilcox()`](https://ruitreves.github.io/statomatic/reference/run_wilcox.md)
for two groups.

``` r
tukey_res <- run_tukey(sorted$ne, group)
head(tukey_res)
#>                    pvalue_Female_KO_vs_Female_WT padj_Female_KO_vs_Female_WT
#> ENSMUSG00000051285                     0.9790520                   0.9998409
#> ENSMUSG00000057363                     0.7955281                   0.9998409
#> ENSMUSG00000078184                     0.9699575                   0.9998409
#> ENSMUSG00000061024                     0.8094858                   0.9998409
#> ENSMUSG00000025911                     0.8062227                   0.9998409
#> ENSMUSG00000057335                     0.8020005                   0.9998409
#>                    pvalue_Female_KO_vs_Male_KO padj_Female_KO_vs_Male_KO
#> ENSMUSG00000051285                   0.8652082                 0.9999797
#> ENSMUSG00000057363                   0.9998882                 0.9999797
#> ENSMUSG00000078184                   0.7280883                 0.9999797
#> ENSMUSG00000061024                   0.1413002                 0.9999797
#> ENSMUSG00000025911                   0.9883949                 0.9999797
#> ENSMUSG00000057335                   0.9992630                 0.9999797
#>                    pvalue_Female_KO_vs_Male_WT padj_Female_KO_vs_Male_WT
#> ENSMUSG00000051285                   0.9999995                         1
#> ENSMUSG00000057363                   0.9899203                         1
#> ENSMUSG00000078184                   0.9050970                         1
#> ENSMUSG00000061024                   0.9469754                         1
#> ENSMUSG00000025911                   0.9925533                         1
#> ENSMUSG00000057335                   0.7267167                         1
#>                    pvalue_Female_WT_vs_Male_KO padj_Female_WT_vs_Male_KO
#> ENSMUSG00000051285                   0.6568296                 0.9999904
#> ENSMUSG00000057363                   0.7597698                 0.9999904
#> ENSMUSG00000078184                   0.9310571                 0.9999904
#> ENSMUSG00000061024                   0.5161829                 0.9999904
#> ENSMUSG00000025911                   0.6260176                 0.9999904
#> ENSMUSG00000057335                   0.7340525                 0.9999904
#>                    pvalue_Female_WT_vs_Male_WT padj_Female_WT_vs_Male_WT
#> ENSMUSG00000051285                   0.9773126                 0.9999998
#> ENSMUSG00000057363                   0.9258447                 0.9999998
#> ENSMUSG00000078184                   0.9954529                 0.9999998
#> ENSMUSG00000061024                   0.9850262                 0.9999998
#> ENSMUSG00000025911                   0.9228793                 0.9999998
#> ENSMUSG00000057335                   0.9990121                 0.9999998
#>                    pvalue_Male_KO_vs_Male_WT padj_Male_KO_vs_Male_WT       test
#> ENSMUSG00000051285                 0.8699999               0.9999869 tukey-test
#> ENSMUSG00000057363                 0.9820468               0.9999869 tukey-test
#> ENSMUSG00000078184                 0.9824965               0.9999869 tukey-test
#> ENSMUSG00000061024                 0.3321566               0.9999869 tukey-test
#> ENSMUSG00000025911                 0.9330056               0.9999869 tukey-test
#> ENSMUSG00000057335                 0.6535879               0.9999869 tukey-test
```

[`fold_change()`](https://ruitreves.github.io/statomatic/reference/fold_change.md)
calculates the log2 fold change between every pair of groups, and the
mean of each group:

``` r
fc <- fold_change(x, group)
head(fc)
#>                    log2FoldChange_Female_KO_vs_Female_WT
#> ENSMUSG00000051285                            -0.1849281
#> ENSMUSG00000057363                            -0.5285301
#> ENSMUSG00000100555                             1.3719378
#> ENSMUSG00000078184                             0.2579516
#> ENSMUSG00000061024                             0.2521327
#> ENSMUSG00000025911                            -0.4755035
#>                    log2FoldChange_Female_KO_vs_Male_KO
#> ENSMUSG00000051285                           0.4555891
#> ENSMUSG00000057363                           0.0472020
#> ENSMUSG00000100555                           0.4793357
#> ENSMUSG00000078184                           0.7006521
#> ENSMUSG00000061024                           0.7674774
#> ENSMUSG00000025911                           0.2151199
#>                    log2FoldChange_Female_KO_vs_Male_WT
#> ENSMUSG00000051285                         0.005584147
#> ENSMUSG00000057363                        -0.196868100
#> ENSMUSG00000100555                         1.443360593
#> ENSMUSG00000078184                         0.413461170
#> ENSMUSG00000061024                         0.148453255
#> ENSMUSG00000025911                        -0.162232725
#>                    log2FoldChange_Female_WT_vs_Male_KO
#> ENSMUSG00000051285                           0.6405172
#> ENSMUSG00000057363                           0.5757321
#> ENSMUSG00000100555                          -0.8926021
#> ENSMUSG00000078184                           0.4427005
#> ENSMUSG00000061024                           0.5153447
#> ENSMUSG00000025911                           0.6906234
#>                    log2FoldChange_Female_WT_vs_Male_WT
#> ENSMUSG00000051285                          0.19051226
#> ENSMUSG00000057363                          0.33166197
#> ENSMUSG00000100555                          0.07142283
#> ENSMUSG00000078184                          0.15550958
#> ENSMUSG00000061024                         -0.10367949
#> ENSMUSG00000025911                          0.31327074
#>                    log2FoldChange_Male_KO_vs_Male_WT Female_KO_mean
#> ENSMUSG00000051285                        -0.4500049      263.67876
#> ENSMUSG00000057363                        -0.2440701       45.28585
#> ENSMUSG00000100555                         0.9640249       79.36147
#> ENSMUSG00000078184                        -0.2871910       40.73782
#> ENSMUSG00000061024                        -0.6190242      289.13164
#> ENSMUSG00000025911                        -0.3773527      260.09906
#>                    Female_WT_mean Male_KO_mean Male_WT_mean
#> ENSMUSG00000051285      299.73955    192.27781    262.66013
#> ENSMUSG00000057363       65.32297     43.82817     51.90698
#> ENSMUSG00000100555       30.66305     56.92661     29.18199
#> ENSMUSG00000078184       34.06800     25.06574     30.58677
#> ENSMUSG00000061024      242.77061    169.84857    260.85955
#> ENSMUSG00000025911      361.64261    224.06872    291.05543
```

If the data is already log-transformed, use
`fold_change(x, group, log_input = TRUE)`. The fold change is then the
difference between the group means. The main functions have the same
`log_input` argument.

### Data in a SummarizedExperiment

If your data is stored in a SummarizedExperiment (for example from
Bioconductor), pass its assay and its sample information to the main
functions:

``` r
library(SummarizedExperiment)
res <- multigroup_main(assay(se), as.data.frame(colData(se)), factors = c("sex", "genotype"))
```

### A more advanced example

In this example we have six different tables we want to analyze. These
tables are from a metagenomic analysis, and contain abundance
information at the phylum, class, order, family, genus and species
level. We looked at the species data before, in the basic workflow
vignette.

There are 24 samples in this experiment in 3 groups, with 8 samples per
group. We’ll call the groups group1, group2 and group3, and analyze them
with
[`multigroup_main()`](https://ruitreves.github.io/statomatic/reference/multigroup_main.md)
with one factor.

We’re going to apply the statomatic method to all six tables, to
identify differences between groups at each taxonomic level.

We read in the data all together and store it in a list.

``` r
# get all files in the current directory with abundance.csv in the name
a <- list.files(".", pattern = "abundance.csv")
a
#> [1] "class_abundance.csv"   "family_abundance.csv"  "genus_abundance.csv"  
#> [4] "order_abundance.csv"   "phylum_abundance.csv"  "species_abundance.csv"
```

Notice that the file names are in alphabetical order.

``` r
# read them in with read.csv
read_in <- lapply(a, read.csv)
# look at the first three rows of each table
lapply(read_in, head, 3)
#> [[1]]
#>              class sample_1 sample_2 sample_3 sample_4 sample_5 sample_6
#> 1       Clostridia      102      475     2910     1188     1278      324
#> 2      Bacteroidia     1066     1005      383      801      623      628
#> 3 Erysipelotrichia      588      186      141      157      127      123
#>   sample_7 sample_8 sample_9 sample_10 sample_11 sample_12 sample_13 sample_14
#> 1     4350      641     1176       883       426      5880        82       158
#> 2      596      549      259       264       551       270       207       582
#> 3      752      139      129       132       229       249       131       464
#>   sample_15 sample_16 sample_17 sample_18 sample_19 sample_20 sample_21
#> 1      1632       531       384       829      1708      4864       799
#> 2       972       700       477       724      1227       661      1754
#> 3       144       135       121       746       752       171       251
#>   sample_22 sample_23 sample_24
#> 1       841      1262       524
#> 2      1656       688      1523
#> 3       230       410        62
#> 
#> [[2]]
#>             family sample_1 sample_2 sample_3 sample_4 sample_5 sample_6
#> 1  Lachnospiraceae       36      123     1888      426      586       69
#> 2 Oscillospiraceae       39      215      501      487      344      145
#> 3     unclassified       76      204      384      342      421       83
#>   sample_7 sample_8 sample_9 sample_10 sample_11 sample_12 sample_13 sample_14
#> 1     2786      258      808       266       217      3684        38        56
#> 2     1066      255      175       324        66      1301        24        54
#> 3      585      159      310       229        56      1308        45        67
#>   sample_15 sample_16 sample_17 sample_18 sample_19 sample_20 sample_21
#> 1      1099       395        39       136       764      3133       217
#> 2       173        71        32       258       491       896       293
#> 3       369       124        71       286       514      1188       273
#>   sample_22 sample_23 sample_24
#> 1       290       508        98
#> 2       292       480       284
#> 3       288       262       197
#> 
#> [[3]]
#>          genus sample_1 sample_2 sample_3 sample_4 sample_5 sample_6 sample_7
#> 1 unclassified      281      304     1902      564      688      235     2828
#> 2   Dubosiella      583      181      118      145      112      118      665
#> 3  Bacteroides      406       92       91      106      111      185      106
#>   sample_8 sample_9 sample_10 sample_11 sample_12 sample_13 sample_14 sample_15
#> 1      370      758       322       228      3085       114       152      1477
#> 2      131      121       117       113       210       120       460       117
#> 3      156       91        93        78        95        97       421       165
#>   sample_16 sample_17 sample_18 sample_19 sample_20 sample_21 sample_22
#> 1       611       185       510      1002      3544       611       665
#> 2       105       114       737       728       123       245       220
#> 3       129       102       129       128       170       617       677
#>   sample_23 sample_24
#> 1       717       323
#> 2       402        57
#> 3       140       774
#> 
#> [[4]]
#>                order sample_1 sample_2 sample_3 sample_4 sample_5 sample_6
#> 1      Eubacteriales      100      468     2888     1184     1261      318
#> 2      Bacteroidales     1065     1005      380      800      621      626
#> 3 Erysipelotrichales      588      186      141      157      127      123
#>   sample_7 sample_8 sample_9 sample_10 sample_11 sample_12 sample_13 sample_14
#> 1     4338      632     1167       874       426      5858        82       158
#> 2      593      549      258       263       551       270       206       582
#> 3      752      139      129       132       229       249       131       464
#>   sample_15 sample_16 sample_17 sample_18 sample_19 sample_20 sample_21
#> 1      1625       529       382       812      1696      4845       783
#> 2       967       696       477       716      1224       656      1745
#> 3       144       135       121       746       752       171       251
#>   sample_22 sample_23 sample_24
#> 1       824      1245       521
#> 2      1650       687      1520
#> 3       230       410        62
#> 
#> [[5]]
#>           phylum sample_1 sample_2 sample_3 sample_4 sample_5 sample_6 sample_7
#> 1      Bacillota      710      790     3214     1471     1573      495     5279
#> 2   Bacteroidota     1071     1007      385      801      627      631      603
#> 3 Pseudomonadota       14      269      475      372      396       29      197
#>   sample_8 sample_9 sample_10 sample_11 sample_12 sample_13 sample_14 sample_15
#> 1      836     1506      1117      1109      6689       234       668      1948
#> 2      551      261       267       551       274       207       584       978
#> 3      266      402       395       263       352       413        14       462
#>   sample_16 sample_17 sample_18 sample_19 sample_20 sample_21 sample_22
#> 1       712       999      1707      2770      5492      1231      1214
#> 2       708       477       736      1235       670      1764      1671
#> 3       176       409        25        35       516        35       163
#>   sample_23 sample_24
#> 1      1862       630
#> 2       704      1532
#> 3        29        18
#> 
#> [[6]]
#>                            species sample_1 sample_2 sample_3 sample_4 sample_5
#> 1        Lachnospiraceae bacterium       14       33      245      224      167
#> 2          Dubosiella newyorkensis      583      181      118      145      112
#> 3 Parasutterella excrementihominis        4      241      431      327      343
#>   sample_6 sample_7 sample_8 sample_9 sample_10 sample_11 sample_12 sample_13
#> 1        9      996       71      312        85        91      1731        22
#> 2      118      665      131      121       117       113       210       120
#> 3       21      128      241      367       364       246       255       393
#>   sample_14 sample_15 sample_16 sample_17 sample_18 sample_19 sample_20
#> 1        12       221       289        11        35       168      1907
#> 2       460       117       105       114       737       728       123
#> 3         8       386       120       385         7         4       433
#>   sample_21 sample_22 sample_23 sample_24
#> 1        42       103       136        28
#> 2       245       220       402        57
#> 3         4       117        10         1
```

Each of these tables needs the cf function.

``` r
# apply cf to each table
tables <- lapply(read_in, cf)
```

We can use the same sample_info for all of these tables, since they all
contain the same samples.

``` r
# tables[[1]] is the first table in the tables list. all the tables have the same column names
sample_info <- data.frame(sample = colnames(tables[[1]]),
                          group = c(rep("group1", 8), rep("group2", 8), rep("group3", 8)))
sample_info
#>       sample  group
#> 1   sample_1 group1
#> 2   sample_2 group1
#> 3   sample_3 group1
#> 4   sample_4 group1
#> 5   sample_5 group1
#> 6   sample_6 group1
#> 7   sample_7 group1
#> 8   sample_8 group1
#> 9   sample_9 group2
#> 10 sample_10 group2
#> 11 sample_11 group2
#> 12 sample_12 group2
#> 13 sample_13 group2
#> 14 sample_14 group2
#> 15 sample_15 group2
#> 16 sample_16 group2
#> 17 sample_17 group3
#> 18 sample_18 group3
#> 19 sample_19 group3
#> 20 sample_20 group3
#> 21 sample_21 group3
#> 22 sample_22 group3
#> 23 sample_23 group3
#> 24 sample_24 group3
```

Now we’ll apply the multigroup_main function to each of these tables,
again using lapply. `verbose = FALSE` turns off the progress messages,
which would otherwise be printed for every table, and
[`set.seed()`](https://rdrr.io/r/base/Random.html) makes the Dunnett T3
p-values reproducible.

``` r
set.seed(1)
res <- lapply(tables, multigroup_main, sample_info = sample_info, factors = "group", verbose = FALSE)
```

We know from before that multigroup_main returns a list, and we just
applied this function to a list of data.frames. So what is res? It’s a
list of lists.

``` r
# res has 6 elements, one per table we analyzed
length(res)
#> [1] 6
# each element of res has 8 elements
length(res[[1]])
#> [1] 8
names(res[[1]])
#> [1] "anova_results" "welch_results" "kw_results"    "foldchanges"  
#> [5] "ne"            "nu"            "nn"            "tests"
```

Each element of res is the same as it was when we were only analyzing
one table. The only difference is that now we’re working with several of
them at once.

We’ll look at the results of the first one just to orient ourselves, and
we’re only going to look at the first three elements: the ANOVA, Welch
and Kruskal-Wallis results.

``` r
# get the first element of res (a list of the results of the first analysis).
# This is the class abundances, since it is first alphabetically
x <- res[[1]]
lapply(x[1:3], head)
#> $anova_results
#>                         pvalue      padj pvalue_group1_vs_group2
#> Betaproteobacteria  0.20000097 0.2317104               0.8263342
#> Coriobacteriia      0.06770511 0.2016507               0.5440510
#> Negativicutes       0.04907064 0.2016507               0.6193168
#> Tissierellia        0.15297673 0.2317104               0.1656226
#> Alphaproteobacteria 0.08962253 0.2016507               0.9493171
#> Flavobacteriia      0.16490962 0.2317104               0.5901316
#>                     padj_group1_vs_group2 log2FoldChange_group1_vs_group2
#> Betaproteobacteria              0.9493171                      -0.2803858
#> Coriobacteriia                  0.9493171                      -0.5572681
#> Negativicutes                   0.9493171                       0.6147098
#> Tissierellia                    0.7453016                       1.2537566
#> Alphaproteobacteria             0.9493171                       0.2630344
#> Flavobacteriia                  0.9493171                      -0.8744691
#>                     pvalue_group1_vs_group3 padj_group1_vs_group3
#> Betaproteobacteria                0.4456629             0.5896895
#> Coriobacteriia                    0.3576415             0.5896895
#> Negativicutes                     0.2415234             0.5434278
#> Tissierellia                      0.9467810             0.9467810
#> Alphaproteobacteria               0.1795866             0.5387598
#> Flavobacteriia                    0.1411272             0.5387598
#>                     log2FoldChange_group1_vs_group3 pvalue_group2_vs_group3
#> Betaproteobacteria                        0.8590642              0.18512474
#> Coriobacteriia                            1.3940111              0.05543407
#> Negativicutes                            -0.6890709              0.04158533
#> Tissierellia                              0.1468414              0.27729400
#> Alphaproteobacteria                      -1.0000000              0.10326409
#> Flavobacteriia                           -1.4150375              0.59013161
#>                     padj_group2_vs_group3 log2FoldChange_group2_vs_group3
#> Betaproteobacteria              0.2776871                       1.1394500
#> Coriobacteriia                  0.1663022                       1.9512792
#> Negativicutes                   0.1663022                      -1.3037807
#> Tissierellia                    0.3565209                      -1.1069152
#> Alphaproteobacteria             0.2323442                      -1.2630344
#> Flavobacteriia                  0.6638981                      -0.5405684
#>                     group1_mean group2_mean group3_mean       test
#> Betaproteobacteria      229.000     278.125     126.250 tukey-test
#> Coriobacteriia          173.125     254.750      65.875 tukey-test
#> Negativicutes             6.125       4.000       9.875 tukey-test
#> Tissierellia              3.875       1.625       3.500 tukey-test
#> Alphaproteobacteria       1.500       1.250       3.000 tukey-test
#> Flavobacteriia            0.750       1.375       2.000 tukey-test
#> 
#> $welch_results
#>                   pvalue       padj pvalue_group1_vs_group2
#> Bacteroidia   0.03093208 0.04254664               0.2372377
#> Chitinophagia 0.04254664 0.04254664               0.6030916
#>               padj_group1_vs_group2 log2FoldChange_group1_vs_group2
#> Bacteroidia               0.4744753                       0.5706097
#> Chitinophagia             0.6030916                       1.2223924
#>               pvalue_group1_vs_group3 padj_group1_vs_group3
#> Bacteroidia                 0.2150517             0.2604673
#> Chitinophagia               0.2604673             0.2604673
#>               log2FoldChange_group1_vs_group3 pvalue_group2_vs_group3
#> Bacteroidia                        -0.6241665              0.03432600
#> Chitinophagia                      -1.2801079              0.05631293
#>               padj_group2_vs_group3 log2FoldChange_group2_vs_group3 group1_mean
#> Bacteroidia              0.05631293                       -1.194776     706.375
#> Chitinophagia            0.05631293                       -2.502500       0.875
#>               group2_mean group3_mean         test
#> Bacteroidia       475.625    1088.750 dunnett-test
#> Chitinophagia       0.375       2.125 dunnett-test
#> 
#> $kw_results
#>                         pvalue      padj pvalue_group1_vs_group2
#> Clostridia          0.81668648 0.8538086                       1
#> Erysipelotrichia    0.65079423 0.7878035                       1
#> unclassified        0.73334806 0.8215781                       1
#> Verrucomicrobiae    0.12045536 0.5288267                       1
#> Bacilli             0.08206044 0.5288267                       1
#> Gammaproteobacteria 0.79132952 0.8465386                       1
#>                     padj_group1_vs_group2 log2FoldChange_group1_vs_group2
#> Clostridia                              1                      0.06548115
#> Erysipelotrichia                        1                      0.45625701
#> unclassified                            1                     -0.28237212
#> Verrucomicrobiae                        1                     -0.47375524
#> Bacilli                                 1                     -2.15324626
#> Gammaproteobacteria                     1                     -0.46247849
#>                     pvalue_group1_vs_group3 padj_group1_vs_group3
#> Clostridia                       1.00000000                     1
#> Erysipelotrichia                 1.00000000                     1
#> unclassified                     1.00000000                     1
#> Verrucomicrobiae                 0.19709373                     1
#> Bacilli                          0.08117608                     1
#> Gammaproteobacteria              1.00000000                     1
#>                     log2FoldChange_group1_vs_group3 pvalue_group2_vs_group3
#> Clostridia                              0.007316499               1.0000000
#> Erysipelotrichia                       -0.309751171               1.0000000
#> unclassified                           -0.561999137               1.0000000
#> Verrucomicrobiae                       -0.652276613               0.2581564
#> Bacilli                                -2.106153552               0.4871113
#> Gammaproteobacteria                    -0.173648087               1.0000000
#>                     padj_group2_vs_group3 log2FoldChange_group2_vs_group3
#> Clostridia                              1                     -0.05816465
#> Erysipelotrichia                        1                     -0.76600818
#> unclassified                            1                     -0.27962702
#> Verrucomicrobiae                        1                     -0.17852137
#> Bacilli                                 1                      0.04709271
#> Gammaproteobacteria                     1                      0.28883041
#>                     group1_mean group2_mean group3_mean      test
#> Clostridia             1408.500    1346.000    1401.375 dunn-test
#> Erysipelotrichia        276.625     201.625     342.875 dunn-test
#> unclassified            145.125     176.500     214.250 dunn-test
#> Verrucomicrobiae         82.000     113.875     128.875 dunn-test
#> Bacilli                  25.375     112.875     109.250 dunn-test
#> Gammaproteobacteria      21.500      29.625      24.250 dunn-test
```

With one factor, the results of all three sets have the same columns, so
we can combine each table’s results into a single data.frame.

For each list of results in res, we’ll rbind (row bind) the first three
elements, save this combined table to a new list and then give it a
name.

For the names, we’ll use the result of the list.files function above,
but we’ll remove the .csv extension first.

``` r
# remove .csv from the file names
name_list <- sub("\\.csv$", "", a)
name_list
#> [1] "class_abundance"   "family_abundance"  "genus_abundance"  
#> [4] "order_abundance"   "phylum_abundance"  "species_abundance"
```

Now we’ll combine and save the tables.

``` r
# define an empty list for our new tables
collated_tables <- list()
# use a for loop to combine and save our tables
for (i in seq_along(res)) {
  # the current result we're processing
  temp <- res[[i]]
  # get just the p-values and fold change results
  first_three <- temp[1:3]
  # rbind into one new table
  new_table <- rbind(first_three[[1]], first_three[[2]], first_three[[3]])
  # save the new table to the collated_tables list
  collated_tables[[i]] <- new_table
}
# set the names of collated_tables to be the names from name_list
names(collated_tables) <- name_list
# look at some of the results
lapply(collated_tables, head, 3)
#> $class_abundance
#>                        pvalue      padj pvalue_group1_vs_group2
#> Betaproteobacteria 0.20000097 0.2317104               0.8263342
#> Coriobacteriia     0.06770511 0.2016507               0.5440510
#> Negativicutes      0.04907064 0.2016507               0.6193168
#>                    padj_group1_vs_group2 log2FoldChange_group1_vs_group2
#> Betaproteobacteria             0.9493171                      -0.2803858
#> Coriobacteriia                 0.9493171                      -0.5572681
#> Negativicutes                  0.9493171                       0.6147098
#>                    pvalue_group1_vs_group3 padj_group1_vs_group3
#> Betaproteobacteria               0.4456629             0.5896895
#> Coriobacteriia                   0.3576415             0.5896895
#> Negativicutes                    0.2415234             0.5434278
#>                    log2FoldChange_group1_vs_group3 pvalue_group2_vs_group3
#> Betaproteobacteria                       0.8590642              0.18512474
#> Coriobacteriia                           1.3940111              0.05543407
#> Negativicutes                           -0.6890709              0.04158533
#>                    padj_group2_vs_group3 log2FoldChange_group2_vs_group3
#> Betaproteobacteria             0.2776871                        1.139450
#> Coriobacteriia                 0.1663022                        1.951279
#> Negativicutes                  0.1663022                       -1.303781
#>                    group1_mean group2_mean group3_mean       test
#> Betaproteobacteria     229.000     278.125     126.250 tukey-test
#> Coriobacteriia         173.125     254.750      65.875 tukey-test
#> Negativicutes            6.125       4.000       9.875 tukey-test
#> 
#> $family_abundance
#>                     pvalue       padj pvalue_group1_vs_group2
#> Sutterellaceae 0.197457096 0.26213702               0.8104024
#> Atopobiaceae   0.066157338 0.12138411               0.5388124
#> Prevotellaceae 0.001577777 0.02208887               0.6487139
#>                padj_group1_vs_group2 log2FoldChange_group1_vs_group2
#> Sutterellaceae             0.8727411                      -0.2966753
#> Atopobiaceae               0.8727411                      -0.5663191
#> Prevotellaceae             0.8727411                       0.7670228
#>                pvalue_group1_vs_group3 padj_group1_vs_group3
#> Sutterellaceae              0.45486902             0.7740066
#> Atopobiaceae                0.35585660             0.7740066
#> Prevotellaceae              0.01360419             0.1904586
#>                log2FoldChange_group1_vs_group3 pvalue_group2_vs_group3
#> Sutterellaceae                       0.8589077              0.18096687
#> Atopobiaceae                         1.4272001              0.05406228
#> Prevotellaceae                      -1.2857835              0.00171722
#>                padj_group2_vs_group3 log2FoldChange_group2_vs_group3
#> Sutterellaceae            0.23284092                        1.155583
#> Atopobiaceae              0.13291457                        1.993519
#> Prevotellaceae            0.02404108                       -2.052806
#>                group1_mean group2_mean group3_mean       test
#> Sutterellaceae      219.00     269.000     120.750 tukey-test
#> Atopobiaceae        168.75     249.875      62.750 tukey-test
#> Prevotellaceae       24.25      14.250      59.125 tukey-test
#> 
#> $genus_abundance
#>                        pvalue      padj pvalue_group1_vs_group2
#> Parasutterella     0.19869563 0.3462475               0.8113967
#> Muricaecibacterium 0.06812654 0.2240055               0.5464864
#> Anaerotruncus      0.51580712 0.5731190               0.8242216
#>                    padj_group1_vs_group2 log2FoldChange_group1_vs_group2
#> Parasutterella                 0.8881362                      -0.2961576
#> Muricaecibacterium             0.8881362                      -0.5635121
#> Anaerotruncus                  0.8881362                       0.3156143
#>                    pvalue_group1_vs_group3 padj_group1_vs_group3
#> Parasutterella                   0.4560574                     1
#> Muricaecibacterium               0.3573489                     1
#> Anaerotruncus                    0.8353710                     1
#>                    log2FoldChange_group1_vs_group3 pvalue_group2_vs_group3
#> Parasutterella                           0.8580840              0.18218500
#> Muricaecibacterium                       1.4364879              0.05582506
#> Anaerotruncus                           -0.2503234              0.48405009
#>                    padj_group2_vs_group3 log2FoldChange_group2_vs_group3
#> Parasutterella                 0.4186984                       1.1542416
#> Muricaecibacterium             0.2469730                       2.0000000
#> Anaerotruncus                  0.6181444                      -0.5659377
#>                    group1_mean group2_mean group3_mean       test
#> Parasutterella         218.875     268.750     120.750 tukey-test
#> Muricaecibacterium     163.750     242.000      60.500 tukey-test
#> Anaerotruncus           35.625      28.625      42.375 tukey-test
#> 
#> $order_abundance
#>                      pvalue      padj pvalue_group1_vs_group2
#> Burkholderiales  0.19886795 0.2471578               0.8240061
#> Coriobacteriales 0.06626659 0.1664696               0.5465790
#> Bacillales       0.04799237 0.1664696               0.1124440
#>                  padj_group1_vs_group2 log2FoldChange_group1_vs_group2
#> Burkholderiales              0.9888073                      -0.2828907
#> Coriobacteriales             0.9888073                      -0.5561823
#> Bacillales                   0.8523957                       1.1020982
#>                  pvalue_group1_vs_group3 padj_group1_vs_group3
#> Burkholderiales                0.4457946             0.7371492
#> Coriobacteriales               0.3504407             0.7371492
#> Bacillales                     0.9373723             0.9373723
#>                  log2FoldChange_group1_vs_group3 pvalue_group2_vs_group3
#> Burkholderiales                        0.8616358              0.18382079
#> Coriobacteriales                       1.4224441              0.05427840
#> Bacillales                            -0.1202942              0.05767428
#>                  padj_group2_vs_group3 log2FoldChange_group2_vs_group3
#> Burkholderiales              0.2545714                        1.144526
#> Coriobacteriales             0.1730229                        1.978626
#> Bacillales                   0.1730229                       -1.222392
#>                  group1_mean group2_mean group3_mean       test
#> Burkholderiales      228.500     278.000     125.750 tukey-test
#> Coriobacteriales     170.875     251.250      63.750 tukey-test
#> Bacillales            20.125       9.375      21.875 tukey-test
#> 
#> $phylum_abundance
#>                             pvalue      padj pvalue_group1_vs_group2
#> Pseudomonadota          0.21409442 0.3211416               0.7870908
#> Actinomycetota          0.07304152 0.2191246               0.5632674
#> Thermodesulfobacteriota 0.80778822 0.8077882               0.8845479
#>                         padj_group1_vs_group2 log2FoldChange_group1_vs_group2
#> Pseudomonadota                      0.9175777                      -0.2956677
#> Actinomycetota                      0.9175777                      -0.5316898
#> Thermodesulfobacteriota             0.9175777                       0.2995603
#>                         pvalue_group1_vs_group3 padj_group1_vs_group3
#> Pseudomonadota                        0.5019108             0.6290021
#> Actinomycetota                        0.3622041             0.6290021
#> Thermodesulfobacteriota               0.9863854             0.9863854
#>                         log2FoldChange_group1_vs_group3 pvalue_group2_vs_group3
#> Pseudomonadota                               0.71426786              0.19346553
#> Actinomycetota                               1.32040063              0.06022461
#> Thermodesulfobacteriota                     -0.08746284              0.80490509
#>                         padj_group2_vs_group3 log2FoldChange_group2_vs_group3
#> Pseudomonadota                      0.2901983                       1.0099356
#> Actinomycetota                      0.1806738                       1.8520905
#> Thermodesulfobacteriota             0.8049051                      -0.3870231
#>                         group1_mean group2_mean group3_mean       test
#> Pseudomonadota               252.25     309.625     153.750 tukey-test
#> Actinomycetota               177.00     255.875      70.875 tukey-test
#> Thermodesulfobacteriota        2.00       1.625       2.125 tukey-test
#> 
#> $species_abundance
#>                                       pvalue        padj
#> Parasutterella excrementihominis 0.199182490 0.298773734
#> Muricaecibacterium torontonense  0.068126536 0.122627766
#> Muribaculaceae bacterium         0.001032886 0.006281909
#>                                  pvalue_group1_vs_group2 padj_group1_vs_group2
#> Parasutterella excrementihominis               0.8062039                     1
#> Muricaecibacterium torontonense                0.5464864                     1
#> Muribaculaceae bacterium                       0.4923547                     1
#>                                  log2FoldChange_group1_vs_group2
#> Parasutterella excrementihominis                      -0.3011695
#> Muricaecibacterium torontonense                       -0.5635121
#> Muribaculaceae bacterium                               0.6818240
#>                                  pvalue_group1_vs_group3 padj_group1_vs_group3
#> Parasutterella excrementihominis              0.46121140            0.63860040
#> Muricaecibacterium torontonense               0.35734892            0.53602338
#> Muribaculaceae bacterium                      0.01419407            0.06145687
#>                                  log2FoldChange_group1_vs_group3
#> Parasutterella excrementihominis                       0.8531586
#> Muricaecibacterium torontonense                        1.4364879
#> Muribaculaceae bacterium                              -1.0093379
#>                                  pvalue_group2_vs_group3 padj_group2_vs_group3
#> Parasutterella excrementihominis             0.182048409           0.297897397
#> Muricaecibacterium torontonense              0.055825056           0.100485100
#> Muribaculaceae bacterium                     0.000975412           0.007926544
#>                                  log2FoldChange_group2_vs_group3 group1_mean
#> Parasutterella excrementihominis                        1.154328      217.00
#> Muricaecibacterium torontonense                         2.000000      163.75
#> Muribaculaceae bacterium                               -1.691162       38.50
#>                                  group2_mean group3_mean       test
#> Parasutterella excrementihominis     267.375     120.125 tukey-test
#> Muricaecibacterium torontonense      242.000      60.500 tukey-test
#> Muribaculaceae bacterium              24.000      77.500 tukey-test
```

From here, we could do some further analysis or filtering of the
results. One way to export them is to write.csv each of the elements of
collated_tables, with the same file name as the original abundance
table, with results\_ at the beginning.

``` r
for (i in seq_along(collated_tables)) {
  write.csv(collated_tables[[i]], paste0("results_", a[i]))
}
```

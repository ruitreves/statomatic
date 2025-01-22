# statomatic
This package was developed to provide a straight-forward and high-throughput method to perform statistical analyses on numeric data. 

## Installation Instructions
Statomatic is currently only available for download from this github site, and some of it's dependencies must be installed manually. The following code can be copied
and pasted into R. Please follow the order of the steps, and do them one at a time.

### Step 1: Install BiocManager
```{r install BiocManager}
#install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.20")
```

### Step 2: Install SummarizedExperiment 
If you are prompted by the following step to choose whether or not to install newer versions of some packages "from source", say no. This could lead to compilation errors and cause statomatic to not be installed.

```{r install SummarizedExperiment}
#Use BiocManager to install the dependency "SummarizedExperiment"
BiocManager::install("SummarizedExperiment")
```

### Step 3: Install remotes and statomatic
```{r install remotes and statomatic}
#install remotes
install.packages("remotes")
#use remotes to install statomatic
remotes::install_github("ruitreves/statomatic")
```

## Introduction
Certain experimental designs are often suited to the preplanned statistical analysis that will take place after data is collected. It is vital, however, to apply the appropriate statistical techniques in order to ensure the validity of the results. The purpose of this package is to provide a high-throughput way to generate statistical results based on characteristics of the data itself, rather than on assumptions about how certains types of data "usually" are. 

As an example, consider a researcher who hopes to perform a two way anova test on a dataset. The researcher may plug their data in and generate statistical results, but without consideration for the mathematical assumptions of the two way anova test, their results could be erroneous. To be suitable for a two way anova, the residuals of the data must be normally distributed and the data must be homoscedastic. The statomatic package is equipped to address these assumptions and sort the data based on which tests can be appropriately applied, and apply them. 

### Background
The statomatic package is specifically designed to handle many seperate hypothesis tests where the aim is to look for statistical difference between
the observations (or groups of observations) for a given variable. For instance, for some gene X (the row), one could be interested in how it is expressed by individuals (the columns) subject to different experimental factors. This is just one hypothesis test, and in this scenario one would likely want to test for significance for several thousand genes. This is the purpose of statomatic, to be a high-throughput hypothesis testing package. 

What sets statomatic apart from other similar programs is it's ability to sort data based on which statistical test is mathematically appropriate for the variable being tested. This means that a different statistical test may be applied to different variables in the dataset. Statomatic does not apply a "one-size-fits-all" approach to an entire data set, but rather handles each individual variable separately. 

See the details section for in depth information on how statomatic sorts and analyzes data.

## Details

The aim of the program is to apply the appropriate statistical test to each variable in a data set. Each row is assumed to be a variable, and columns are assumed to be observations of that variable.

1. For each row in a data set, a generalized linear model is created based on the specified factor(s) (groups, etc.). 
From this model, the residuals are tested using a Shapiro Normality Test with alpha = 0.05. 
Rows that are deemed by this criteria to have normally distributed residuals are then tested by Levene's Test to assess the equality of variances (homoscedacity) amongst the specified factors with alpha = 0.05.
 Rows that are deemed by these criteria to be both normally distributed and homoscedastic are used to construct a new data set, "ne" (for normal, equal variance). 
 Rows that are deemed to be normally distributed but heteroscedastic (non-homoscedastic) are used to construct the data set "nu" (for normal, unequal variance). 
 Rows that are deemed to be non-normally distributed are used to construct the data set "nn" (for non-normal).
 
2. In the case of multiple groups (> 2 groups):
      ne is analyzed by a either a one-way or a two-way anova test, depending on your experimental factors. A p-value is reported for each factor and the interaction of multiple factos, if applicable. A Tukey multi comparison test is done based on factor3, which is assumed to be some linear combination of factors 1 and 2 (e.g., if factor1 = A, factor2 = B, then factor 3 should be something like AB, either literally or symbolically). P values are reported for each groupwise comparison. 
      nu is analyzed by a Welch's t-test based on factor3, and a p value is reported for this test. A Dunnett T3 multicomparison test is done based also on factor3, and a p value is reported for each groupwise comparison. 
      nn is analyzed by a Kruskal-Wallis rank sums test and a p value is reported. This is followed by a Dunn multicomparison test and a p value is reported for each groupwise comparison.
      If you have only one experimental factor, it will be used by every test. 

   In the case of two groups:
     ne is analyzed by a Student's t-test and a p value is reported.
     nu is analyzed by a Welch's t-test and a p value is reported.
     nn is analyzed by a Wilcoxon/Mann-Whitney test and a p value is reported. 

3. After the statistical tests are performed, log2 fold changes are calculated. Note that the program will sort the group names alphabetically before calculating the fold changes. The fold changes will be reported as log2FoldChange_GroupA_vs_GroupB, which indicates that the following operation was performed: log2(mean(groupA) / mean(groupB))

4. Your results are ready. 

## See vignettes here

Basic workflow vignette: https://ruitreves.github.io/statomatic/articles/basic_workflow.html

Additional workflows vignette: https://ruitreves.github.io/statomatic/articles/additional_workflows.html

## Quick Start Guide

Load the package:

```{r}
library(statomatic)
```

Read in data:

```{r} 
x <- read.csv("normalized_counts.csv")
head(x)
```

Format data correctly:

```{r}
#set column 1 of x to be rownames of x and remove column 1, using statomatic::cf()
x <- cf(x)
head(x)
```

Create sample metadata to define experimental factors:

```{r}
sample_info <- data.frame(sample = colnames(x), sex = c(rep("m", 10), rep("f", 10)), 
                          genotype = c(rep("wt", 5), rep("ko", 5), rep("wt", 5), rep("ko", 5)),
                          group = c(rep("m_wt", 5), rep("m_ko", 5), rep("f_wt", 5), rep("f_ko", 5)))
sample_info
```

Build statomatic data set

```{r}
sds <- make_sds(x, colData = sample_info, design = ~ sex + genotype)
sds
```

Run analysis:

```{r} 
sds <- sds_analyze(sds)
```

Additionally, statomatic includes stand-alone functions which can be used on data.frames. These functions will treat each row as a seperate variable and each
column as an observation. 

```{r all funcs}
#sort data according to statomatic's procedure, one factor or two
res1 <- test_norm(x, sample_info$group)
res2 <- test_norm(x, sample_info$sex, sample_info$genotype)

#anova, can use either one factor or two
anova_res1 <- run_anova(x, sample_info$group)
anova_res2 <- run_anova(x, sample_info$sex, sample_info$genotype)

#welch's t-test and kruskal-wallis tests, one factor only
welch_res <- run_welch(x, sample_info$group)
kw_res <- run_kruskal(x, sample_info$group)

#Multi-comparison tests: TukeyHSD, DunnettT3, and Dunn's test, one factor only
tukey_res <- run_tukey(x, sample_info$group)
dunnett_res <- run_dunnett(x, sample_info$group)
dunn_res <- run_dunn(x, sample_info$group)

#t-tests, one factor only, strictly for two groups
t_res <- run_ttest(x, sample_info$group)
wilcox_res <- run_wilcox(x, sample_info$group)

#calculating fold changed, one factor only 
fc <- fold_change(x, sample_info$group)
```

## Limitations

Statomatic aims to make addressing the mathematical assumptions of hypothesis tests as easy as possible for users, but be aware that not all kinds of data can be treated in the same way. In multiple examples shown here, bulk RNA-sequencing data is used to demonstrate how to use statomatic. However, there are highly specialized tools specifically designed to handle RNA-seq data, such as DESeq2, which perform many operations and checks on the data that statomatic does not. Although this package draws much inspiration from DESeq2, statomatic should not be seen as a replacement for the highly specialized procedures performed by other tools designed to handle specific data challenges.





<!-- README.md is generated from README.Rmd. Please edit that file -->

<https://github.com/JiujingWu/AHD>

# Adaptive Hellinger Distance (AHD) for Testing Independence between two random vectors

<!-- badges: start -->

[![R-CMD-check](https://github.com/wzhong41/semidist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/wzhong41/semidist/actions/workflows/R-CMD-check.yaml)
<!-- badges: end --> The goal of package `AHD` is to provide an easy way
to implement the Adaptive Hellinger Distance.

## Installation

To install `AHD`,

``` r
install.packages("AHD")
```

## Example

Here is a simple example showing how to use `AHD` to compute the test
statistic between two random vectors, and apply the test statistic on
testing the independence and conduct p-values combination.

``` r
library(AHD)

n <- 50
x <- runif(n, 0, 1)
y <- x * (sample(c(-1, 1), size = n, replace = TRUE)) + runif(n, -0.5, 0.5)

AHD <- AHD(x, y)
AHD_4bins_pval <- AHD_test_4bins(x, y, B = 199)
p_AHD <- AHD_4bins_pval$pvalue
p_AHD_11 <- AHD_4bins_pval$pvalue11
p_AHD_12 <- AHD_4bins_pval$pvalue12
p_AHD_21 <- AHD_4bins_pval$pvalue21
p_AHD_22 <- AHD_4bins_pval$pvalue22

p_val_com <- P_VALUE_COMB(p_AHD_11, p_AHD_12, p_AHD_21, p_AHD_22)
```

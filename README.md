
# AHD

<!-- badges: start -->
<!-- badges: end -->

The goal of package `AHD` is to provide an easy way to implement the Adaptive Hellinger Distance (AHD). 

Key features:

- A rank-based independence test using the AHD.

- The AHD-based tests address the potential zeros in the denominators of the HHG test statistic.

- The AHD-based tests inherit the robustness of the HHG test while conservatively improving its performance.

- The AHD tests can be regarded as a crucial extension of the HHG test. 

- By default, the Manhattan ($L_1$) norm is used as the distance metric for both $X$ and $Y$ spaces. 

- Distance ties are handled using the max-rank method with stable sorting, and self-distances are excluded when computing ranks.


## Installation

To install `AHD`, 
``` r
remotes::install_github("JiujingWu/AHD")
```

## Example

Here is a simple example showing how to use `AHD` to compute the test statistic between two random vectors, 
and apply the test statistic on testing the independence and conduct p-values combination. 

``` r
library(AHD)

B <- 300
n <- 50
x <- runif(n, 0, 1)
y <- x * (sample(c(-1, 1), size = n, replace = TRUE)) + runif(n, -0.5, 0.5)

perms <- t(replicate(B, sample.int(n), simplify = "matrix"))

Dx <- as.matrix(dist(x, method = "manhattan", diag = T, upper = T))
Dy <- as.matrix(dist(y, method = "manhattan", diag = T, upper = T))
  
AHD_stat <- AHD(Dx, Dy)
AHD_4bins_pval <- AHD_test(Dx, Dy, perms)
p_AHD <- AHD_4bins_pval$pvalue
p_AHD_11 <- AHD_4bins_pval$pvalue11
p_AHD_12 <- AHD_4bins_pval$pvalue12
p_AHD_21 <- AHD_4bins_pval$pvalue21
p_AHD_22 <- AHD_4bins_pval$pvalue22

p_val_com <- P_VALUE_COMB(p_AHD_11, p_AHD_12, p_AHD_21, p_AHD_22)
```


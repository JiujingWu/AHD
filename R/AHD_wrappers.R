#' Adaptive Hellinger Distance (AHD) statistic
#'
#' Computes the AHD statistic (and its four quadrant components) from two
#' n x n distance matrices.
#'
#' @param Dx numeric n x n distance matrix for X.
#' @param Dy numeric n x n distance matrix for Y.
#' @return A list with components `AHD`, `AHD11`, `AHD12`, `AHD21`, `AHD22`.
#' @details Distances can be any metric; Dx and Dy must share the same n.
#' @examples
#' \donttest{
#' set.seed(1); n <- 20
#' x <- runif(n); y <- rnorm(n)
#' Dx <- as.matrix(dist(x, method = "manhattan", diag = TRUE, upper = TRUE))
#' Dy <- as.matrix(dist(y, method = "manhattan", diag = TRUE, upper = TRUE))
#' AHD(Dx, Dy)
#' }
#' @export
AHD <- function(Dx, Dy) {
  .Call(`_AHD_AHD`, Dx, Dy)
}

#' Permutation test for AHD
#'
#' Performs a permutation test using a precomputed permutation table
#' (rows are permutations of 1..n). The unbiased p-value formula (k+1)/(B+2)
#' is used.
#'
#' @param Dx numeric n x n distance matrix for X.
#' @param Dy numeric n x n distance matrix for Y.
#' @param perms integer matrix B x n; each row a permutation of 1..n.
#' @return A list with p-values: `pvalue`, `pvalue11`, `pvalue12`, `pvalue21`, `pvalue22`.
#' @examples
#' \donttest{
#' set.seed(1); n <- 20; B <- 50
#' x <- runif(n); y <- rnorm(n)
#' Dx <- as.matrix(dist(x, method = "manhattan", diag = TRUE, upper = TRUE))
#' Dy <- as.matrix(dist(y, method = "manhattan", diag = TRUE, upper = TRUE))
#' perms <- t(replicate(B, sample.int(n)))
#' AHD_test(Dx, Dy, perms)
#' }
#' @export
AHD_test <- function(Dx, Dy, perms) {
  .Call(`_AHD_AHD_test`, Dx, Dy, perms)
}



#' Combine four AHD component p-values
#'
#' @description
#' Combine \eqn{(p_{11}, p_{12}, p_{21}, p_{22})} into a single p-value using
#' the Cauchy Combination Test (CCT) and the Benjamini–Hochberg (BH) minimum.
#'
#' @param pvalue11,pvalue12,pvalue21,pvalue22 Numeric vectors of equal length:
#' permutation p-values of the four quadrant components.
#'
#' @return A list with:
#' * \code{CCT}: combined p-values via the Cauchy Combination Test;
#' * \code{BH}:  combined p-values via the BH-adjusted minimum (Simes-type).
#'
#' @references
#' Liu, Y. and Xie, J. (2019). Cauchy Combination Test: A Powerful Test With Analytic p-Value Calculation Under Arbitrary Dependency Structures.
#' \emph{Journal of the American Statistical Association}, 115(529), 393–402.
#'
#' @examples
#' set.seed(1)
#' p11 <- runif(5); p12 <- runif(5); p21 <- runif(5); p22 <- runif(5)
#' P_VALUE_COMB(p11, p12, p21, p22)
#'
#' @name P_VALUE_COMB
NULL

#' @rawNamespace export(P_VALUE_COMB)
NULL

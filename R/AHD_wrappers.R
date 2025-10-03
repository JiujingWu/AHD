#' Adaptive Hellinger Distance (AHD): statistic and permutation test
#'
#' @description
#' \code{AHD()} computes the AHD statistic (and its four quadrant
#' components) from two \eqn{n \times n} distance matrices.
#' \code{AHD_test()} performs a permutation test using a precomputed
#' permutation table (rows are permutations of \code{1:n}).
#'
#' @param Dx Numeric \eqn{n \times n} distance matrix for \eqn{X}.
#' @param Dy Numeric \eqn{n \times n} distance matrix for \eqn{Y}.
#' @param perms Integer matrix of size \eqn{B \times n}; each row is a
#'   permutation of \code{1:n}. (**Only for** \code{AHD_test()}).
#'
#' @return
#' * \code{AHD()}: list with \code{AHD}, \code{AHD11}, \code{AHD12},
#'   \code{AHD21}, \code{AHD22}.
#' * \code{AHD_test()}: list with permutation p-values \code{pvalue},
#'   \code{pvalue11}, \code{pvalue12}, \code{pvalue21}, \code{pvalue22}.
#'
#' @details
#' Permutation p-values use the unbiased formula \eqn{(k+1)/(B+2)}.
#' Ties in ranks are handled by the max-rank rule.
#' Distances can be any metric, but \code{Dx} and \code{Dy} must both be
#' \eqn{n \times n} with the same \eqn{n}.
#'
#' @usage AHD(Dx, Dy)
#' @usage AHD_test(Dx, Dy, perms)
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 20; B <- 50
#' x <- runif(n); y <- rnorm(n)
#' Dx <- as.matrix(dist(x, method = "manhattan", diag = T, upper = T))
#' Dy <- as.matrix(dist(y, method = "manhattan", diag = T, upper = T))
#' # observed statistic
#' AHD_stat = AHD(Dx, Dy)
#' # permutation test (rows are permutations of 1..n)
#' perms <- t(replicate(B, sample.int(n)))
#' AHD_perm_test = AHD_test(Dx, Dy, perms)
#' }
#'
#' @name AHD
#' @aliases AHD AHD_test
NULL



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

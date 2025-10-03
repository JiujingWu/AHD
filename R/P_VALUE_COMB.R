
#' Adaptive Hellinger Distance (AHD) independence test
#'
#' @description Implement the AHD independence test via the combined p-values.
#'
#' @param pvalue11 The p-value computed from the test statistic \eqn{T_{11}}.
#' @param pvalue12 The p-value computed from the test statistic \eqn{T_{12}}.
#' @param pvalue21 The p-value computed from the test statistic \eqn{T_{21}}.
#' @param pvalue22 The p-value computed from the test statistic \eqn{T_{22}}.
#'
#' @returns A list of the objects about the combined p-values:
#' * `BH`: Combined p-values are adjusted via the Benjamini–Hochberg (BH) method;
#' * `CCT`: Combined p-values are adjusted via the Cauchy Combination Test (CCT) method.
#' @export
P_VALUE_COMB = function(pvalue11, pvalue12, pvalue21, pvalue22){
  p_values = cbind(pvalue11, pvalue12, pvalue21, pvalue22)
  CCT = apply(
    p_values, 1, function(x) {
      cauchy_stat = mean(tan((0.5 - x) * pi))
      0.5 - atan(cauchy_stat) / pi
    }
  )

  Simes = apply(
    p_values, 1, function(x) {min( p.adjust(x, method = "BH", n = length(x)) )}
  )

  list(CCT = CCT, Simes = Simes)
}

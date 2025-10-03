#' Pairwise angular-distance matrix around a reference row
#'
#' @description Computes the pairwise angular-distance matrix (in radians) for rows of \code{A} centered at a reference row.
#'
#' @param A A numeric matrix of size \eqn{n \times p} with sample size \eqn{n} and dimension \eqn{p}.
#' @param r An integer in \code{1:n} specifying the reference row about which the row-wise centering is performed.
#'
#' @returns A numeric matrix of angular distances.
#'
#' @seealso \link[=pcov2]{pcov2}, \link[=pcov2_test]{pcov2_test}
#' @export
arccos = function(A, r){
  A = as.matrix(A)
  n = dim(A)[1]; p = dim(A)[2]
  A1 = A - matrix(1,n,1)%*%matrix(A[r, ],1,p)
  A1 = A1[-r, ]
  if(p > 1){
    A1_std = A1 / (matrix(sqrt(rowSums(A1^2)),n-1,1)%*%matrix(1,1,p))
  } else if(p == 1){
    A1_std = A1 / (matrix(sqrt(A1^2),n-1,1)%*%matrix(1,1,p))
  }
  A1_euc = suppressWarnings(acos(A1_std%*%t(A1_std)))
  A1_euc[is.nan(A1_euc)] = 0
  return(A1_euc)
}

#' Projection Covariance
#'
#' @description Calculate the projection covariance of two random vectors. Two vectors can be with different dimensions, but with equal sample sizes.
#'
#' @param x A numeric matrix of size n * p.
#' @param y A numeric matrix of size n * q.
#'
#' @returns A list with two elements:
#' \item{pcov}{The projection covariance of x and y.}
#' \item{S2}{A component of the projection covariance, which is required for calculating the subsequent test statistic.}
#'
#' @references L.Zhu, K.Xu, R.Li, W.Zhong (2017). Projection correlation between two random vectors. Biometrika, Volume 104, Pages 829-843. https://doi.org/10.1093/biomet/asx043
#' @export
pcov2 = function(x, y){
  x = as.matrix(x); y = as.matrix(y)
  n = dim(x)[1]; p = dim(x)[2]; q = dim(y)[2]
  pcov = S2 = c()
  for(r in 1:n){
    X = arccos(x, r); Y = arccos(y, r)
    XY1 = X%*%Y
    XY2 = X*Y
    S1 = sum(XY2)/n^2
    S2[r] = sum(X)*sum(Y)/n^4
    S3 = sum(XY1)/n^3
    pcov[r] = S1 + S2[r] - 2*S3
  }
  return(list(pcov=pcov, S2 = S2))
}


#' Projection Correlatoin Permutation Test
#'
#' @description Test whether two vectors are independent. Calculate the test result of projection correlation test.
#'
#' @param x A numeric matrix of size n * p.
#' @param y A numeric matrix of size n * q.
#' @param B Integer, the number of permutations. Defaults to \code{299}.
#'
#' @returns A list with:
#' \item{stat}{The observed test statistic \code{stat}.}
#' \item{pvalue}{The permutation p-value \code{pvalue}.}
#'
#' @references L.Zhu, K.Xu, R.Li, W.Zhong (2017). Projection correlation between two random vectors. Biometrika, Volume 104, Pages 829-843. https://doi.org/10.1093/biomet/asx043
#' @export
pcov2_test = function(x, y, B=299){
  x = as.matrix(x); y = as.matrix(y)
  n = dim(x)[1]
  ppcov = pcov2(x,y)
  pcov_n = mean(ppcov$pcov) / (pi^2-mean(ppcov$S2))
  num_rej = 0
  for (b in 1:B) {
    ppcov_star = pcov2(x[sample(1:n), ], y)
    pcov_star = mean(ppcov_star$pcov) / (pi^2-mean(ppcov_star$S2))
    if (pcov_star >= pcov_n) {
      num_rej = num_rej + 1
    }
  }
  pvalue = (num_rej+1) / (B+1)
  return(list(stat=pcov_n, pvalue=pvalue))
}


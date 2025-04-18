#' Adaptive Hellinger Distance (AHD)
#'
#' @description This function computes the AHD between two random vectors.
#'
#' @param x A matrix containing the first random vectors.
#' @param y A matrix containing the second random vectors.
#' @param method The method for calculating the distance between two random vectors.
#' This must be one of `"euclidean"`, `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
#' Any unambiguous substring can be given. Defaults to `"manhattan"`.
#'
#' @returns A list of the objects about the test statistics:
#' * `AHD`: The computed test statistic \eqn{T} value;
#' * `AHD11`: The computed test statistic \eqn{T_{11}} value;
#' * `AHD12`: The computed test statistic \eqn{T_{12}} value;
#' * `AHD21`: The computed test statistic \eqn{T_{21}} value;
#' * `AHD22`: The computed test statistic \eqn{T_{22}} value.
#' @export
AHD = function(x, y, method = "manhattan"){
  x = as.matrix(x); y = as.matrix(y)
  n = dim(x)[1]
  Dx = as.matrix(dist(x, method = method, diag = T, upper = T))
  ORD.x = t(apply(Dx, 1, order))
  A0 = matrix(rep(1:n,n),n,n,byrow=T)
  A1_ =  (A0-2)*(A0>2)/(n-2)

  Dy = as.matrix(dist(y, method = method, diag = T, upper = T))
  indx = matrix(c(rep(1:n,n),as.vector(ORD.x)), , 2)
  Dy = Dy[indx]
  Dy = matrix(Dy,n,n)
  A0 = t(apply(Dy, 1, rank))
  A_1 = (A0-2)*(A0>2)/(n-2)
  A11 = matrix(0,n,n)
  A22 = matrix(0,n,n)
  for(j in 1:n){
    A11[,j] = rowSums(A_1[,1:j] <= matrix(rep(A_1[,j],j),n,j))-2
    A22[,j] = rowSums(A_1[,j:n] > matrix(rep(A_1[,j],n-j+1),n,n-j+1))
  }
  A11 = A11*(A11 > 0)/(n-2)
  A22 = A22/(n-2)
  A2_ = 1 - A1_; A_2 = 1 - A_1
  A12 = A1_ - A11; A21 = A_1 - A11
  A11[A11<0] = 0; A12[A12<0] = 0; A21[A21<0] = 0; A22[A22<0] = 0
  AHD = 4*(n-2)*(sum((sqrt(A11)-sqrt(A_1*A1_))^2 + (sqrt(A22)-sqrt(A_2*A2_))^2 +
                       (sqrt(A12)-sqrt(A1_*A_2))^2 + (sqrt(A21)-sqrt(A2_*A_1))^2))^(1/2)

  AHD11 = 4*(n-2)*(sum( (sqrt(A11)-sqrt(A_1*A1_))^2 ))^(1/2)
  AHD12 = 4*(n-2)*(sum( (sqrt(A12)-sqrt(A1_*A_2))^2 ))^(1/2)
  AHD21 = 4*(n-2)*(sum( (sqrt(A21)-sqrt(A2_*A_1))^2 ))^(1/2)
  AHD22 = 4*(n-2)*(sum( (sqrt(A22)-sqrt(A_2*A2_))^2 ))^(1/2)

  return(list(AHD = AHD, AHD11 = AHD11, AHD12 = AHD12, AHD21 = AHD21, AHD22 = AHD22 ))
}


#' Adaptive Hellinger Distance (AHD) independence test
#'
#' @description Implement the AHD independence test via permutation test.
#'
#' @param x A matrix containing the first random vectors.
#' @param y A matrix containing the second random vectors.
#' @param method The method for calculating the distance between two random vectors.
#' This must be one of `"euclidean"`, `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
#' Any unambiguous substring can be given. Defaults to `"manhattan"`.
#'
#' @returns A list of the objects about the p-values:
#' * `pvalue`: The p-value computed from the test statistic \eqn{T};
#' * `pvalue11`: The p-value computed from the test statistic \eqn{T_{11}};
#' * `pvalue12`: The p-value computed from the test statistic \eqn{T_{12}};
#' * `pvalue21`: The p-value computed from the test statistic \eqn{T_{21}};
#' * `pvalue22`: The p-value computed from the test statistic \eqn{T_{22}}.
#' @export
AHD_test_4bins = function(x, y, method = "manhattan", B = 299) {
  x = as.matrix(x); y = as.matrix(y)
  n = dim(as.matrix(x))[1]
  AHD_out = AHD(x, y, method)
  AHD_n = AHD_out$AHD
  AHD_n_11 = AHD_out$AHD11
  AHD_n_12 = AHD_out$AHD12
  AHD_n_21 = AHD_out$AHD21
  AHD_n_22 = AHD_out$AHD22
  num_rej = num_rej11 = num_rej12 = num_rej21 = num_rej22 = 0
  for (b in 1:B) {
    per_index = sample(1:n)
    AHD_out_star = AHD(x[per_index, ], y, method)
    AHD_star = AHD_out_star$AHD
    AHD_star11 = AHD_out_star$AHD11
    AHD_star12 = AHD_out_star$AHD12
    AHD_star21 = AHD_out_star$AHD21
    AHD_star22 = AHD_out_star$AHD22
    if (AHD_star >= AHD_n) { num_rej = num_rej + 1 }
    if (AHD_star11 >= AHD_n_11) { num_rej11 = num_rej11 + 1 }
    if (AHD_star12 >= AHD_n_12) { num_rej12 = num_rej12 + 1 }
    if (AHD_star21 >= AHD_n_21) { num_rej21 = num_rej21 + 1 }
    if (AHD_star22 >= AHD_n_22) { num_rej22 = num_rej22 + 1 }
  }
  pvalue = (num_rej+1) / (B+1)
  pvalue11 = (num_rej11+1) / (B+1)
  pvalue12 = (num_rej12+1) / (B+1)
  pvalue21 = (num_rej21+1) / (B+1)
  pvalue22 = (num_rej22+1) / (B+1)
  return(list(pvalue = pvalue, pvalue11 = pvalue11, pvalue12 = pvalue12, pvalue21 = pvalue21, pvalue22 = pvalue22))
}


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
      return(0.5 - atan(cauchy_stat) / pi)
    }
  )
  BH = apply(
    p_values, 1, function(x) {min( p.adjust(x, method = "BH", n = length(x)) )}
  )

  return(list(CCT = CCT, BH = BH))
}

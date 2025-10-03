library(RcppParallel)
library(parallel)
library(mvtnorm)
library(MASS)
library(HHG)
library(energy)
library(dHSIC)
library(IndepTest)
library(AHD)

rm(list = ls())


gen2_indep_norm = function(n, p) {
  X = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)

  Y1 = rnorm(n, sd=5)
  Y2 = rnorm(n, sd=5)
  Y3 = rnorm(n, sd=5)
  Y4 = rnorm(n, sd=5)
  Y5 = rnorm(n, sd=5)
  Y6 = rnorm(n, sd=5)
  Y7 = rnorm(n, sd=5)
  Y8 = rnorm(n, sd=5)

  Y_cont = cbind(Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8)

  lin  = 0.5*X[,1]^2 + 3*cos(X[,2])
  prob = 1 / (1 + exp(-lin))
  Y9 = sample(c(0,1), size=n, replace = TRUE)  # rbinom(n, size = 1, prob = prob)

  Y = cbind(Y_cont, Y9)

  return(list(X = X, Y = Y))
}

gen2_indep_abs_t = function(n, p) {
  X = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)

  Y1 = abs(rt(n, df=1))
  Y2 = abs(rt(n, df=1))
  Y3 = abs(rt(n, df=1))
  Y4 = abs(rt(n, df=1))
  Y5 = abs(rt(n, df=1))
  Y6 = abs(rt(n, df=1))
  Y7 = abs(rt(n, df=1))
  Y8 = abs(rt(n, df=1))

  Y_cont = cbind(Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8)

  lin = 0.5*X[,1]^2 + 3*cos(X[,2])
  prob = 1 / (1 + exp(-lin))
  Y9 = sample(c(0,1), size=n, replace = TRUE)  # rbinom(n, size = 1, prob = prob)

  Y = cbind(Y_cont, Y9)

  return(list(X = X, Y = Y))
}


gen2_mixed_depend_norm = function(n, p) {
  X = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)

  Y1 = 5*sin(X[,1]^2) + rnorm(n, sd = 5)
  Y2 = 3*X[,4] * X[,5] + rnorm(n, sd = 5)
  Y3 = 2*log(abs(X[,6]) + 1) + rnorm(n, sd = 5)
  Y4 = 2 * X[,8] - X[,9]^3 + rnorm(n, sd = 5)
  Y5 = X[,10] * (1 + X[,1]) + rnorm(n, sd = 5)
  Y6 = 4 * X[,2] * X[,3] + sqrt(abs(X[,4])) + rnorm(n, sd = 5)
  Y7 = X[,5]^3 / (1 + X[,6]^2) + rnorm(n, sd = 5)
  Y8 = 3*cos(X[,7]^2) + 2*X[,8] * X[,9] + rnorm(n, sd = 5)

  Y_cont = cbind(Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8)

  lin  = 0.5*X[,1]^2 + 3*cos(X[,2])
  prob = 1 / (1 + exp(-lin))
  Y9 = rbinom(n, size = 1, prob = prob)

  Y = cbind(Y_cont, Y9)

  return(list(X = X, Y = Y))
}

gen2_mixed_depend_t = function(n, p) {
  X = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)

  Y1 = 5*sin(X[,1]^2) + abs(rt(n, df = 1))
  Y2 = 3*X[,4] * X[,5] + abs(rt(n, df = 1))
  Y3 = 2*log(abs(X[,6]) + 1) + abs(rt(n, df = 1))
  Y4 = 2 * X[,8] - X[,9]^3 + abs(rt(n, df = 1))
  Y5 = X[,10] * (1 + X[,1]) + abs(rt(n, df = 1))
  Y6 = 4 * X[,2] * X[,3] + sqrt(abs(X[,4])) + abs(rt(n, df = 1))
  Y7 = X[,5]^3 / (1 + X[,6]^2) + abs(rt(n, df = 1))
  Y8 = 3*cos(X[,7]^2) + 2*X[,8] * X[,9] + abs(rt(n, df = 1))

  Y_cont = cbind(Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8)

  lin  = 0.5*X[,1]^2 + 3*cos(X[,2])
  prob = 1 / (1 + exp(-lin))
  Y9 = rbinom(n, size = 1, prob = prob)

  Y = cbind(Y_cont, Y9)

  return(list(X = X, Y = Y))
}


n = 150
p = 30

Sigma = 0.5^(abs(outer(1:p, 1:p, "-")))
eig = eigen(Sigma); sqA = sqrt(eig$values)
sigma = eig$vectors%*%diag(sqA)%*%t(eig$vectors)
B = 300; R = 1000

p_ahd = p_ahd_11 = p_ahd_12 = p_ahd_21 = p_ahd_22 = p_ahd_sum = c()
p_hhg = p_dc = p_pcov = p_hsic = p_mint = numeric(R)

perms = t(replicate(B, sample.int(n), simplify = "matrix"))

set.seed(928)
for(i in 1:R){
  data = gen2_mixed_depend_t(n, p)  # gen2_mixed_depend_norm(n, p) # gen2_indep_abs_t(n, p)  # gen2_indep_norm(n, p)
  x = data$X; y = data$Y

  Dx = as.matrix(dist(x, method = "manhattan", diag = T, upper = T))
  Dy = as.matrix(dist(y, method = "manhattan", diag = T, upper = T))

  AHD_4bins_pval = AHD_test_cpp_fast_parallel(Dx, Dy, perms)
  p_ahd[i] = AHD_4bins_pval$pvalue
  p_ahd_11[i] = AHD_4bins_pval$pvalue11
  p_ahd_12[i] = AHD_4bins_pval$pvalue12
  p_ahd_21[i] = AHD_4bins_pval$pvalue21
  p_ahd_22[i] = AHD_4bins_pval$pvalue22

  p_hhg[i] = hhg.test(Dx, Dy, nr.perm = B)$perm.pval.hhg.sc
  p_dc[i] = dcorT.test(x,y)$p.value
  p_mint[i] = MINTperm(x, y, k=20, B=B)
  p_hsic[i] = dhsic.test(x, y, alpha = 0.05, B=B, kernel = "gaussian")$p.value
  p_pcov[i] = pcov2_test(x, y, B=B)$pvalue
}
p_values = P_VALUE_COMB(p_ahd_11, p_ahd_12, p_ahd_21, p_ahd_22)

out = matrix(c(p_ahd, p_values$Simes, p_values$CCT, p_hhg, p_dc, p_pcov, p_hsic, p_mint), ncol = R, byrow = T)
pval = function(x){mean(x<=0.05)}
power = apply(out, 1, pval)
names(power) = c("AHD", "Simes", "CCT", "HHG", "DCor", "PCOV", "HSIC", "MINT")
power

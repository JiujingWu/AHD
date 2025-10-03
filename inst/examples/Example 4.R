library(RcppParallel)
library(parallel)
library(mvtnorm)
library(MASS)
library(HHG)
library(energy)
library(dHSIC)
library(semidist)
library(AHD)

rm(list = ls())

gen4 = function(n, p, prob, xtype){

  if(prob == "equal"){
    p_1 = 1 / 2
    p_2 = 1 - p_1
    sample_size1 = ceiling(p_1 * n)
    sample_size2 = n - sample_size1
  } else if(prob == "nonequal"){
    p_1 = 1 / 3
    p_2 = 1 - p_1
    sample_size1 = ceiling(p_1 * n)
    sample_size2 = n - sample_size1
  }

  Y = c(rep(1, sample_size1), rep(2, sample_size2))
  data = matrix(0, nrow = n, ncol = p + 1)
  X = matrix(abs(rt(n = n * p, df = 1)), nrow = n)
  data[ ,1:p] = X
  data[ ,p+1] = Y
  w = as.integer(runif(sample_size1) > 0.5)
  if(xtype == "mixnorm"){
    data[which(data[ ,p+1]==1), c(2*1-1, 2*1)] = 5*matrix(w*rnorm(n=sample_size1, mean = -5, sd=0.1) +
                                                            (1-w)*rnorm(n=sample_size1, mean=5, sd=0.1), ncol=2)
    data[which(data[ ,p+1]==2), c(2*2-1, 2*2)] = 5*matrix(w*rnorm(n=sample_size2, mean = -5, sd=0.1) +
                                                            (1-w)*rnorm(n=sample_size2, mean=5, sd=0.1), ncol=2)
  } else if(xtype == "mixcauchy"){
    data[which(data[ ,p+1] == 1), c(2*1-1, 2*1)] = 5*matrix(w*rcauchy(n=sample_size1, location=-1, scale=10) +
                                                              (1-w)*rcauchy(n=sample_size1, location=1, scale=10), ncol=2)
    data[which(data[ ,p+1] == 2), c(2*2-1, 2*2)] = 5*matrix(w*rcauchy(n=sample_size2, location=-1, scale=10) +
                                                              (1-w)*rcauchy(n=sample_size2, location=1, scale=10), ncol=2)
  }

  y = data[ ,p+1]
  x = data[ ,1:p]
  return(list(x=x, y=y))
}


gen4_indep = function(n, p, prob){
  if(prob == "equal"){
    p_1 = 1 / 2
    p_2 = 1 - p_1
    sample_size1 = ceiling(p_1 * n)
    sample_size2 = n - sample_size1
  } else if(prob == "nonequal"){
    p_1 = 1 / 3
    p_2 = 1 - p_1
    sample_size1 = ceiling(p_1 * n)
    sample_size2 = n - sample_size1
  }

  y = c(rep(1, sample_size1), rep(2, sample_size2))
  data = matrix(0, nrow = n, ncol = p + 1)
  x = matrix(abs(rt(n = n * p, df = 1)), nrow = n)

  return(list(x=x, y=y))

}

n = 60; p = 100
B = 300; R = 1000

p_ahd = p_ahd_11 = p_ahd_12 = p_ahd_21 = p_ahd_22 = c()
p_hhg = p_dc = p_mint = p_hsic = p_pcov = p_semidc = c()

perms = t(replicate(B, sample.int(n), simplify = "matrix"))

set.seed(928)
for(i in 1:R){
  data = gen4(n, p, prob = "equal", xtype = "mixnorm")
  # data = gen4_indep(n, p, prob = "equal")
  x = data$x; y = data$y

  Dx = as.matrix(dist(x, method = "manhattan", diag = T, upper = T))
  Dy = as.matrix(dist(y, method = "manhattan", diag = T, upper = T))

  AHD_4bins_pval = AHD_test(Dx, Dy, perms)
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
  p_semidc[i] = sd_test(x, as.factor(y), "perm", num_perm=B)$pvalue
}

p_values = P_VALUE_COMB(p_ahd_11, p_ahd_12, p_ahd_21, p_ahd_22)

out = matrix(c(p_ahd, p_values$Simes, p_values$CCT, p_hhg, p_dc, p_pcov, p_hsic, p_semidc), ncol = R, byrow = T)
pval = function(x){ mean( x <= 0.05 ) }
power = apply(out, 1, pval)
names(power) = c("AHD", "Simes", "CCT", "HHG", "DCor", "PCOV", "HSIC", "SemiDC")
power


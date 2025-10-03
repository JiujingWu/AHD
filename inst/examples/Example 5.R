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

gen2 = function(model){

  ##############################################
  #################### size ####################
  ##############################################
  if(model == "4.1.2"){
    x = matrix(rbinom(n*p, size=10, prob=0.5), n, p);
    y = matrix(rbinom(n*p, size=4, prob=0.5), n, p)
  } else if (model == "4.1.3"){
    z1 = matrix(rbinom(n*p, size=10, prob=0.5), n, p);
    z2 = matrix(rbinom(n*p, size=4, prob=0.5), n, p)
    x = z1^2+z1; y = z2^2
  } else if (model == "4.1.4"){
    x = mvrnorm(n, mu=rep(0,p), Sigma = Sigma1)
    y = mvrnorm(n, mu=rep(0,p), Sigma = Sigma1)
  } else if(

    ##############################################
    #################### power ###################
    ##############################################
    model == "4.1.5"){  # x discrete; y continuous
    x1 = matrix(rbinom(n*(p-1), 10, 0.5), n, p-1); z = matrix(rbinom(n, 20, 0.5), n, 1)
    y = sin(x1[ ,1])*(x1[ ,2]+z)^2/(x1[ ,3]+1) + log(x1[ ,4]+z+1)
    x = cbind(x1, z)
  } else if(model == "4.1.6"){  # x continuous; y continuous
    x = matrix(rt(n*p, df=1), n, p)
    y1 = 3*(4-x[ ,1])^2*x[ ,3]*x[ ,4] + x[ ,2]^2
    y2 = 3*(4-x[ ,2])^2*x[ ,3]*x[ ,4] + x[ ,1]^2
    y = cbind(y1, y2)
  } else if(model == "4.1.9"){   # x continuous; y continuous
    x = matrix(rt(n*p, df=1), n, p)
    y1 = 0.5*x[ ,1]*x[ ,2]
    y2 = 0.5*(x[ ,3]+x[ ,4]+x[ ,5])^2
    y3 = 0.5*sin(x[ ,1]^2) + 0.5*log(x[ ,5]^2)
    y = cbind(y1, y2, y3)
  }

  return(list(x=x, y=y))
}


n = 50; p = 10
Sigma1 = matrix(0.5, ncol=p, nrow = p)
diag(Sigma1) = 1
B = 300; R = 1000

p_ahd = p_ahd_11 = p_ahd_12 = p_ahd_21 = p_ahd_22 = p_ahd_sum = c()
p_hhg = p_dc = p_pcov = p_hsic = p_mint = numeric(R)

perms = t(replicate(B, sample.int(n), simplify = "matrix"))

set.seed(928)
for(i in 1:R){
  data = gen2(model = "4.1.5")
  x = data$x; y = data$y

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

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


gen1 = function(p, sigma, m, xtype, model, beta){

  if (m < p){


    y1 = matrix(rnorm(n*p), n, p)
    if(xtype == "norm"){
      x = matrix(rnorm(n*p), n, p) %*% sigma + matrix(1, n, 1) %*% matrix(0, 1, p)
      epsilon = matrix(rnorm(n*p), n, p)

      if(model == "erci"){
        y = cbind(beta*x[ ,1:m]^2, y1[ ,(m+1):p]) + epsilon
      } else if(model == "log"){
        y = cbind(beta*log(x[ ,1:m]^2), y1[ ,(m+1):p]) + epsilon
      } else if(model == "sin"){
        y = cbind(beta*sin(x[ ,1:m]^2), y1[ ,(m+1):p]) + epsilon
      } else if(model == "linear"){
        y = cbind(beta*x[ ,1:m], y1[ ,(m+1):p]) + epsilon
      }
    } else if(xtype == "t"){
      x = (matrix(rnorm(n*p), n, p) %*% sigma + matrix(1, n, 1) %*% matrix(0, 1, p)) / (sqrt(rchisq(n, df=2))/2)
      epsilon = matrix(rnorm(n*p), n, p)
      if(model == "erci"){
        y = cbind(beta*x[ ,1:m]^2, y1[ ,(m+1):p]) + epsilon
      } else if(model == "log"){
        y = cbind(beta*log(x[ ,1:m]^2), y1[ ,(m+1):p]) + epsilon
      } else if(model == "sin"){
        y = cbind(beta*sin(x[ ,1:m]), y1[ ,(m+1):p]) + epsilon
      } else if(model == "linear"){
        y = cbind(beta*x[ ,1:m], y1[ ,(m+1):p]) + epsilon
      }
    }


  }
  else if (m == p){


    if(xtype == "norm"){
      x = matrix(rnorm(n*p), n, p) %*% sigma + matrix(1, n, 1) %*% matrix(0, 1, p)
      epsilon = matrix(rnorm(n*p), n, p)

      if(model == "erci"){
        y = beta*x^2 + epsilon
      } else if(model == "log"){
        y = beta*log(x^2) + epsilon
      } else if(model == "sin"){
        y = beta*sin(x[ ,1:m]^2)+ epsilon
      }
    }  else if(xtype == "t"){
      x = (matrix(rnorm(n*p), n, p) %*% sigma + matrix(1, n, 1) %*% matrix(0, 1, p)) / (sqrt(rchisq(n, df=2))/2)
      epsilon = matrix(rnorm(n*p), n, p)
      if(model == "erci"){
        y = beta*x^2 + epsilon
      } else if(model == "log"){
        y = beta*log(x^2) + epsilon
      } else if(model == "sin"){
        y = beta*sin(x) + epsilon
      }
    }


  }
  return(list(x=x, y=y, epsilon=epsilon))
}

gen1_size = function(p, sigma, xtype){
  if(xtype == "norm"){
    x = matrix(rnorm(n*p), n, p) %*% sigma + matrix(1, n, 1) %*% matrix(0, 1, p)
  } else if (xtype == "t"){
    x = (matrix(rnorm(n*p), n, p) %*% sigma + matrix(1, n, 1) %*% matrix(0, 1, p)) / (sqrt(rchisq(n, df=2))/2)
  }
  y = matrix(rnorm(n*p), n, p)
  return(list(x=x, y=y))
}
pval = function(x){ mean(x <= 0.05) }

n = 30; p = 5
rho = 0.5
Sigma1 = rho^(abs(outer(1:p, 1:p, "-")))
eig = eigen(Sigma1); sqA = sqrt(eig$values)
sigma = eig$vectors%*%diag(sqA)%*%t(eig$vectors)

B = 300; R = 1000
p_ahd = p_ahd_11 = p_ahd_12 = p_ahd_21 = p_ahd_22 = p_ahd_sum = c()
p_hhg = p_dc = p_pcov = p_hsic = p_mint = numeric(R)

perms = t(replicate(B, sample.int(n), simplify = "matrix"))

set.seed(928)
for(m in 1:5){

  for(i in 1:R){
    set.seed(i)

    data = gen1(p, Sigma1, m, xtype = "t", model = "sin", beta = 6)
    # data = gen1_size(p, Sigma1, xtype = "t")
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

  p_values = cbind(p_ahd_11, p_ahd_12, p_ahd_21, p_ahd_22)
  cct = apply(
    p_values, 1, function(x) {
      cauchy_stat = mean(tan((0.5 - x) * pi))
      return(0.5 - atan(cauchy_stat) / pi)
    }
  )
  Simes = apply(
    p_values, 1, function(x) {min( p.adjust(x, method = "BH", n = length(x)) )}
  )

  pvals = data.frame(p_rhd, Simes, cct, p_hhg, p_dc, p_pcov1, p_hsic, p_mint)
  cat("******************************", "\n",
      apply(pvals, 2, pval), "\n")
}


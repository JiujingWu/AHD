library(RcppParallel)
library(parallel)
library(mvtnorm)
library(AHD)

rm(list = ls())

datagenSpiral = function(n, noise = 0) {
  theta = runif(n, 0, 6 * pi)
  r = theta + runif(n,0,noise)
  x = r * cos(theta) 
  y = r * sin(theta)  
  return (cbind(x, y))
}

datagenSine = function(n, noise = 1) {
  x = runif(n, 0, 1)
  y = sin(8*pi * x) + runif(n, -noise, noise)
  return (cbind(x, y))
}

datagencubic = function(n, noise = 2){
  x = rnorm(n)
  y = x^3 + noise*rt(n, df=2)
  return (cbind(x, y))
}

datagenCircle = function(n, r1 = 1, r2 = 3) {
  theta = runif(n, 0, 2*pi)
  r = runif(n, r1, r2)
  x = r * cos(theta)
  y = r * sin(theta)
  return (cbind(x, y))
}


n = 100
B = 300; R = 1000

p_ahd = p_ahd_11 = p_ahd_12 = p_ahd_21 = p_ahd_22 = c()

perms = t(replicate(B, sample.int(n), simplify = "matrix"))

set.seed(928)
for(i in 1:R){
  data = datagenSpiral(n, noise = 0)   # datagenCircle(n, r1 = 1, r2 = 3) # datagencubic(n, noise = 2) # datagenSine(n, noise = 1)
  x = data[ ,1]; y = data[ ,2]
  
  Dx = as.matrix(dist(x, method = "manhattan", diag = T, upper = T))
  Dy = as.matrix(dist(y, method = "manhattan", diag = T, upper = T))
  
  AHD_4bins_pval = AHD_test(Dx, Dy, perms)
  p_ahd[i] = AHD_4bins_pval$pvalue
  p_ahd_11[i] = AHD_4bins_pval$pvalue11
  p_ahd_12[i] = AHD_4bins_pval$pvalue12
  p_ahd_21[i] = AHD_4bins_pval$pvalue21
  p_ahd_22[i] = AHD_4bins_pval$pvalue22
}

out = matrix(c(p_ahd, p_ahd_11, p_ahd_12, p_ahd_21, p_ahd_22), ncol = R, byrow = T)
pval = function(x){mean(x<=0.05)}
power = apply(out, 1, pval)
names(power) = c("AHD", "AHD11", "AHD12", "AHD21", "AHD22")
power

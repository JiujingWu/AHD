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


##example ***
datagenWedge = function(n, noise) {
  x = runif(n, 0, 1)
  y = x * (sample(c(-1, 1), size = n, replace = TRUE)) + runif(n, -noise, noise)
  return (cbind(x, y))
}

#example***
datagenW = function(n, noise){
  x = runif(n, -1, 1)
  u = x + runif(n,-noise,noise)
  v = 4 * ((x ^ 2 - 0.5) ^ 2 + runif(n, -noise, noise) / 500)
  return (cbind(u, v))
}

##example ***
datagenParabola.corrected = function(n,noise) {
  x = runif(n, -1, 1); y = (x^2)/2+runif(n,-noise,noise)
  return (cbind(x, y))
}

##example ***
datagenCircle = function(n, r1, r2) {
  theta = runif(n, 0, 2*pi)
  r = runif(n, r1, r2)
  x = r * cos(theta)
  y = r * sin(theta)
  return (cbind(x, y))
}

##example ***
datagenCross = function(n, noise) {
  x = runif(n, 0, 1)
  s = sample(c(-1, 1), size = n, replace = TRUE)
  y = 2 * s * x - (s-1) + runif(n, -noise, noise)
  return (cbind(x, y))
}

##example ***
datagenCross = function(n, noise) {
  x = runif(n, 0, 1)
  s = sample(c(-1, 1), size = n, replace = TRUE)
  y = 2 * s * x - (s-1) + runif(n, -noise, noise)
  return (cbind(x, y))
}

##example ***
datagenDoppler = function(n, noise) {
  x = runif(n, 0, 1)
  y = sqrt(x * (1-x)) * sin(2.1*pi/(x+0.05)) + runif(n, -noise, noise)
  return (cbind(x, y))
}

##example ***
datagenSine = function(n, noise) {
  x = runif(n, 0, 1)
  y = sin(8 * pi * x) + runif(n, -noise, noise)
  return (cbind(x, y))
}

##example ***
datagenHeavisine = function(n, noise) {
  x = runif(n, 0, 1)
  y = (4 * sin(4 * pi * x) - sign(x - 0.3) - sign (0.72 - x)) + runif(n, -noise, noise)
  return (cbind(x, y))
}

##example ***
datagenindep = function(n){
  x = runif(n, 0, 1)
  y = runif(n, -1, 1)
  return (cbind(x, y))
}


n = 50
B = 300; R = 1000

p_ahd = p_ahd_11 = p_ahd_12 = p_ahd_21 = p_ahd_22 = c()
p_hhg = p_dc = p_pcov = p_hsic = p_mint = c()

perms = t(replicate(B, sample.int(n), simplify = "matrix"))

set.seed(928)
for(i in 1:R){
  data = datagenWedge(n, 0.6)
  # data = datagenW(n, 0.4)
  # data = datagenParabola.corrected(n, 0.3)
  # data = datagenCircle(n, 1, 1.8)
  # data = datagenCross(n, 0.4)
  # data = datagenDoppler(n, 0.3)
  # data = datagenSine(n, 0.4)
  # data = datagenHeavisine(n, 4)
  # data = datagenindep(n)
  x = data[ ,1]; y = data[ ,2]

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
}
p_values = P_VALUE_COMB(p_ahd_11, p_ahd_12, p_ahd_21, p_ahd_22)

out = matrix(c(p_ahd, p_values$Simes, p_values$CCT, p_hhg, p_dc, p_pcov, p_hsic, p_mint), ncol = R, byrow = T)
pval = function(x){mean(x<=0.05)}
power = apply(out, 1, pval)
names(power) = c("AHD", "Simes", "CCT", "HHG", "DCor", "PCOV", "HSIC", "MINT")
power

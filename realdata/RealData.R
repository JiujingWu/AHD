
# Load packages
library(HellCor)
library(HHG)
library(mvtnorm)
library(energy)
library(dHSIC)
library(minerva)
library(MASS)
library(IndepTest)
library(AHD)

# calculate for PCOV
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

pcov2 = function(x,y){
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

# p-values combination
p_val_combin = function(x, y, K){
  ahd_4bins_pval = AHD_test_4bins(x, y, B = K)
  p_ahd_11 = ahd_4bins_pval$pvalue11
  p_ahd_12 = ahd_4bins_pval$pvalue12
  p_ahd_21 = ahd_4bins_pval$pvalue21
  p_ahd_22 = ahd_4bins_pval$pvalue22
  p_values = cbind(p_ahd_11, p_ahd_12, p_ahd_21, p_ahd_22)

  cct = apply(p_values, 1, function(x) {
    cauchy_stat = mean(tan((0.5 - x) * pi))
    0.5 - atan(cauchy_stat) / pi
  })
  BH = apply(p_values, 1, function(x) { min(p.adjust(x, method = "BH")) })

  return(list(
    p_ahd=ahd_4bins_pval$pvalue,
    p_ahd_11=p_ahd_11, p_ahd_12=p_ahd_12, p_ahd_21=p_ahd_21, p_ahd_22=p_ahd_22,
    cct=cct, BH=BH))
}

calc_pvals = function(x, y, K = 300) {

  Dx = as.matrix(dist(x, method = "manhattan", diag = TRUE, upper = TRUE))
  Dy = as.matrix(dist(y, method = "manhattan", diag = TRUE, upper = TRUE))

  p_hhg = hhg.test(Dx, Dy, nr.perm = K)$perm.pval.hhg.sc
  dcor_res = dcor.test(x, y, R = K)
  p_dc = dcor_res$p.value
  p_mint = MINTperm(x, y, k=5, B=K)
  p_hsic = dhsic.test(x, y, alpha = 0.05, B=K, kernel = "gaussian")$p.value
  p_pcov = pcov2_test(x, y, B=K)$pvalue
  p_combination = p_val_combin(x, y, K)
  p_ahd = p_combination$p_ahd
  cct = p_combination$cct
  bh = p_combination$BH

  return(c(AHD = p_ahd, CCT = cct, BH = bh, HHG = p_hhg, DC = p_dc, PCOV = p_pcov, HSIC = p_hsic, MINT = p_mint))
}


run_all_tests = function(df, K = 300) {
  vars = colnames(df)
  res_list = list()

  for (i in 1:(length(vars)-1)) {
    for (j in (i+1):length(vars)) {
      x = df[[i]]
      y = df[[j]]

      pvals = calc_pvals(x, y, K)
      res_list[[paste(vars[i], "vs", vars[j])]] = pvals
    }
  }

  res_df = do.call(rbind, res_list)
  rownames(res_df) = names(res_list)
  return(res_df)
}

data(Chagos)

df = Chagos[, c("Treatment", "Seabirds_ha", "kg_N_ha_yr", "Number_of_fishes", "biomass")]
df$Treatment = as.numeric(Chagos$Treatment)
x = as.numeric(Treatment); y = cbind(Fish, Biomass)
x = Nitrogen; y = cbind(Fish, Biomass)
x = Seabirds; y = cbind(Fish, Biomass)

K = 300
results_df = run_all_tests(df, K)
print(round(results_df, 3))

x = Chagos$kg_N_ha_yr
y = Chagos[, c("Number_of_fishes", "biomass")]
calc_pvals(x, y, K)
p_val_combin(x, y, K)


#######################################################
###################### SRBCT ##########################
library(plsgenomics)
data("SRBCT")
x = SRBCT[["X"]]
y = SRBCT[["Y"]]
standard_func = function(x){ (x - mean(x)) / sd(x) }
x_std = t(apply(x, 1, standard_func))
x = x_std
n = dim(x)[1]; p = dim(x)[2]
R = length(unique(y))


prob_values = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
K = 300
results = matrix(0, nrow = length(prob_values), ncol = 8)
colnames(results) = c("HHG", "DCor", "HSIC", "MINT", "PCOV", "AHD", "BH", "CCT")

for (i in 1:length(prob_values)) {
  set.seed(i)
  prob = prob_values[i]
  index = sample(1:n, size = floor(prob * n))
  y[index] = sample(y[index])

  Dx = as.matrix(dist(x, method = "manhattan", diag = T, upper = T))
  Dy = as.matrix(dist(y, method = "manhattan", diag = T, upper = T))
  out = p_val_combin(x, y, K)

  results[i, 1] = hhg.test(Dx, Dy, nr.perm = K)$perm.pval.hhg.sc
  results[i, 2] = dcor.test(x, y, R = K)$p.value
  results[i, 3] = dhsic.test(x, y, alpha = 0.05, B = K, kernel = "gaussian")$p.value
  results[i, 4] = MINTperm(x, y, k = 5, B = K)
  results[i, 5] = pcov2_test(x, y, B = K)$pvalue
  results[i, 6] = AHD_test_4bins(x, y, B = K)$pvalue
  results[i, 7] = out$BH
  results[i, 8] = out$cct
}

rownames(results) = prob_values
round(t(results), 3)

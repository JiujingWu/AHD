
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

# p-values combination
p_val_combin = function(Dx, Dy, B){
  perms = t(replicate(B, sample.int(n), simplify = "matrix"))

  ahd_4bins_pval = AHD_test(Dx, Dy, perms)
  p_ahd_11 = ahd_4bins_pval$pvalue11
  p_ahd_12 = ahd_4bins_pval$pvalue12
  p_ahd_21 = ahd_4bins_pval$pvalue21
  p_ahd_22 = ahd_4bins_pval$pvalue22
  p_values = cbind(p_ahd_11, p_ahd_12, p_ahd_21, p_ahd_22)

  cct = apply(p_values, 1, function(x) {
    cauchy_stat = mean(tan((0.5 - x) * pi))
    0.5 - atan(cauchy_stat) / pi
  })
  Simes = apply(p_values, 1, function(x) { min(p.adjust(x, method = "BH")) })

  return(list(
    p_ahd=ahd_4bins_pval$pvalue,
    p_ahd_11=p_ahd_11, p_ahd_12=p_ahd_12, p_ahd_21=p_ahd_21, p_ahd_22=p_ahd_22,
    cct=cct, Simes=Simes))
}

calc_pvals = function(x, y, B = 300) {

  Dx = as.matrix(dist(x, method = "manhattan", diag = TRUE, upper = TRUE))
  Dy = as.matrix(dist(y, method = "manhattan", diag = TRUE, upper = TRUE))

  p_hhg = hhg.test(Dx, Dy, nr.perm = B)$perm.pval.hhg.sc
  dcor_res = dcor.test(x, y, R = B)
  p_dc = dcor_res$p.value
  p_mint = MINTperm(x, y, k=5, B=B)
  p_hsic = dhsic.test(x, y, alpha = 0.05, B=B, kernel = "gaussian")$p.value
  p_pcov = pcov2_test(x, y, B=B)$pvalue
  p_combination = p_val_combin(Dx, Dy, B)
  p_ahd = p_combination$p_ahd
  cct = p_combination$cct
  simes = p_combination$Simes

  return(c(AHD = p_ahd, CCT = cct, Simes = simes, HHG = p_hhg, DC = p_dc, PCOV = p_pcov, HSIC = p_hsic, MINT = p_mint))
}


run_all_tests = function(df, B = 300) {
  vars = colnames(df)
  res_list = list()

  for (i in 1:(length(vars)-1)) {
    for (j in (i+1):length(vars)) {
      x = df[[i]]
      y = df[[j]]

      pvals = calc_pvals(x, y, B)
      res_list[[paste(vars[i], "vs", vars[j])]] = pvals
    }
  }

  res_df = do.call(rbind, res_list)
  rownames(res_df) = names(res_list)
  return(res_df)
}

###################### Chagos ##########################

data(Chagos)

df = Chagos[, c("Treatment", "Seabirds_ha", "kg_N_ha_yr", "Number_of_fishes", "biomass")]
df$Treatment = as.numeric(df$Treatment)
n = length(df$Treatment)

B = 300
results_df = run_all_tests(df, B)
print(round(results_df, 3))


###################### SRBCT ##########################

rm(list = ls())
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
B = 200
results = matrix(0, nrow = length(prob_values), ncol = 8)
colnames(results) = c("HHG", "DCor", "HSIC", "MINT", "PCOV", "AHD", "Simes", "CCT")

for (i in 1:length(prob_values)) {
  set.seed(i)
  prob = prob_values[i]
  index = sample(1:n, size = floor(prob * n))
  y[index] = sample(y[index])

  Dx = as.matrix(dist(x, method = "manhattan", diag = T, upper = T))
  Dy = as.matrix(dist(y, method = "manhattan", diag = T, upper = T))
  out = p_val_combin(Dx, Dy, B)

  results[i, 1] = hhg.test(Dx, Dy, nr.perm = B)$perm.pval.hhg.sc
  results[i, 2] = dcor.test(x, y, R = B)$p.value
  results[i, 3] = dhsic.test(x, y, alpha = 0.05, B = B, kernel = "gaussian")$p.value
  results[i, 4] = MINTperm(x, y, k = 5, B = B)
  results[i, 5] = pcov2_test(x, y, B = B)$pvalue
  results[i, 6] = out$p_ahd
  results[i, 7] = out$Simes
  results[i, 8] = out$cct
}

rownames(results) = prob_values
round(t(results), 3)

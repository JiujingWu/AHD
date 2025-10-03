
library(HHG)

rm(list = ls())

datagenSine  = function(n, sig, noise) {
  x = runif(n, 0, 1)
  y = sig * sin(8 * pi * x) + runif(n, -noise, noise)
  return(cbind(x, y))
}


n = 200
prob = c(0.1, 0.5, 0.9)
sig.num = 8
signal = seq(0, 0.5, length = sig.num)
rep = 500
num_perm = 2000  

result = matrix(0, length(prob), 2 * sig.num)


permutation_test = function(data, x0, y0, n, num_perm) {
  s1_permuted = numeric(num_perm)
  ss1_permuted = numeric(num_perm)
  
  for (m in 1:num_perm) {
    x = data[, 1]
    y_permut = data[sample(1:n), 2]
    s1_tmp = ss1_tmp = 0
    
    for (j in 1:n) {
      Rx = abs(x0 - data[j, 1])
      Ry = abs(y0 - y_permut[j])
      
      A11 = sum((abs(x-x0) <= Rx) * (abs(y_permut-y0) <= Ry))
      A12 = sum((abs(x-x0) <= Rx) * (abs(y_permut-y0) > Ry))
      A21 = sum((abs(x-x0) > Rx) * (abs(y_permut-y0) <= Ry))
      A22 = sum((abs(x-x0) > Rx) * (abs(y_permut-y0) > Ry))
      
      A1. = max(A11+A12, 1)
      A.1 = max(A11+A21, 1)
      A2. = max(A21+A22, 1)
      A.2 = max(A12+A22, 1)
      
      s1_tmp = s1_tmp + (n * (A11*A22 - A12*A21)^2 / A1./A.1/A2./A.2)
      ss1_tmp = ss1_tmp + (4 * (sqrt(A11) - sqrt(A1.*A.1 / n))^2 +
                             4 * (sqrt(A12) - sqrt(A1.*A.2 / n))^2 +
                             4 * (sqrt(A21) - sqrt(A2.*A.1 / n))^2 +
                             4 * (sqrt(A22) - sqrt(A2.*A.2 / n))^2)
    }
    
    s1_permuted[m] = s1_tmp
    ss1_permuted[m] = ss1_tmp
  }
  
  return(list(s1_threshold = quantile(s1_permuted, 0.95),
              ss1_threshold = quantile(ss1_permuted, 0.95)))
}

for (jjj in 1:length(prob)) {
  power1 = power2 = NULL
  for (kkk in 1:length(signal)) {
    s1 = ss1 = numeric(rep)
    
    for (m in 1:rep) {
      data = datagenSine(n, signal[kkk], 0.4)
      x = data[, 1]
      y = data[, 2]
      x0 = quantile(x, prob[jjj])
      y0 =quantile(y, prob[jjj])  # quantile(y, 1 - prob[jjj])
      
      thresholds = permutation_test(data, x0, y0, n, num_perm)
      
      for (j in 1:n) {
        Rx = abs(x0 - x[j])
        Ry = abs(y0 - y[j])
        
        A11 = sum((abs(x-x0) <= Rx) * (abs(y-y0) <= Ry))
        A12 = sum((abs(x-x0) <= Rx) * (abs(y-y0) > Ry))
        A21 = sum((abs(x-x0) > Rx) * (abs(y-y0) <= Ry))
        A22 = sum((abs(x-x0) > Rx) * (abs(y-y0) > Ry))
        
        A1. = max(A11+A12, 1)
        A.1 = max(A11+A21, 1)
        A2. = max(A21+A22, 1)
        A.2 = max(A12+A22, 1)
        
        s1[m] = s1[m] + (n * (A11 * A22 - A12 * A21)^2 / A1./ A.1/ A2./ A.2)
        ss1[m] = ss1[m] + (4 * (sqrt(A11) - sqrt(A1.*A.1 / n))^2 +
                             4 * (sqrt(A12) - sqrt(A1.*A.2 / n))^2 +
                             4 * (sqrt(A21) - sqrt(A2.*A.1 / n))^2 +
                             4 * (sqrt(A22) - sqrt(A2.*A.2 / n))^2)
      }
    }
    
    power1 = c(power1, mean(s1 >= thresholds$s1_threshold))
    power2 = c(power2, mean(ss1 >= thresholds$ss1_threshold))
  }
  result[jjj, ] = c(power1, power2)
}

result

HHG = result[ ,1:sig.num]
AHD = result[ ,(sig.num+1):(2*sig.num)]

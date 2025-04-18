test_that("multiplication works", {

  n = 50
  x <- runif(n, 0, 1)
  y <- x * (sample(c(-1, 1), size = n, replace = TRUE)) + runif(n, -0.5, 0.5)

  AHD = AHD(x, y)
  AHD_4bins_pval = AHD_test_4bins(x, y, B = 199)
  p_AHD = AHD_4bins_pval$pvalue
  p_AHD_11 = AHD_4bins_pval$pvalue11
  p_AHD_12 = AHD_4bins_pval$pvalue12
  p_AHD_21 = AHD_4bins_pval$pvalue21
  p_AHD_22 = AHD_4bins_pval$pvalue22

  p_val_com = P_VALUE_COMB(p_AHD_11, p_AHD_12, p_AHD_21, p_AHD_22)

  expect_type(AHD_4bins_pval, "list")
  expect_type(p_val_com, "list")
  expect_length(p_AHD_11, 1)

})

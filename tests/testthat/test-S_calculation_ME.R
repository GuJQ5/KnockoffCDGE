test_that("S_calculation_ME", {
  data("Sigma")
  S_calculation_ME(Sigma,M=5,verbose=TRUE,tol=0.001)
})

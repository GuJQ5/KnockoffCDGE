test_that("Low rank approximation works", {
  data("Sigma")
  S_calculation_PCA(Sigma,r=5)
})

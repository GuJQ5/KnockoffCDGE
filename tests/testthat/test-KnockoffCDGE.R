test_that("KnockoffCDGE works", {
  data("DGE_result")
  data("Sigma")
  Z<-Z_calculation(DGE_result$pvalue,DGE_result$log2FoldChange)
  KnockoffCDGE(Z,Sigma,M=5,n=23010,method="ME",Gene_info=DGE_result)
})

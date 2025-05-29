test_that("KO filter", {
  set.seed(433)
  tau<-c(runif(20)+0.5,runif(80))
  kappa<-c(sample(0:5,20,prob=c(0.8,rep(0.04,5)),replace=TRUE),sample(0:5,80,replace=TRUE))
  KO_Filter(tau,kappa,M=5)
})

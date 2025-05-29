#' Knockoff q-values computation
#'
#' This function computes knockoff q-values from feature statistics of all features.
#'
#' @param tau Feature statistics (magnitude)
#' @param kappa Feature statistics (sign)
#' @param M The number of knockoff copies (under the multiple knockoff framework)
#'
#' @return Knockoff q-values of all features (the smallest target FDR level for each feature to be selected).
#' @export
#'
#' @examples
#' set.seed(433)
#' tau<-c(runif(20)+0.5,runif(80))
#' kappa<-c(sample(0:5,20,prob=c(0.8,rep(0.04,5)),replace=TRUE),sample(0:5,80,replace=TRUE))
#' KO_Filter(tau,kappa,M=5)
KO_Filter<-function(tau,kappa,M)
{
  ord<-order(tau,decreasing = TRUE)
  tau_ord<-tau[ord]
  kappa_ord<-kappa[ord]
  kappa_ord[tau_ord<1e-6]<-1

  keep<-sum(tau_ord>1e-6)
  tau_ord<-tau_ord[1:keep]
  kappa_ord<-kappa_ord[1:keep]

  RR<-(kappa_ord==0)

  aux_term<-0*RR+1

  RR<-cumsum(kappa_ord==0)
  RR<-RR+(RR==0)
  V<-1/M*aux_term+1/M*cumsum(kappa_ord!=0)

  V_RR<-V/RR
  q_V_RR<-rep(Inf,length(ord))
  V_RR_order<-order(V_RR,decreasing = TRUE)
  for(i in V_RR_order)
  {
    q_V_RR[1:i]<-V_RR[i]
  }
  q_V_RR[kappa_ord!=0]<-Inf

  q<-rep(Inf,length(ord))
  q[ord]<-q_V_RR
  return(q)
}

#' Knockoff-based conditional differential gene expression (CDGE) analysis
#'
#' This function performs the knockoff-based CDGE analysis with Z-scores of DGE analysis and correlation matrix of gene expressions.
#'
#' @param Z DGE Z-scores.
#' @param Sigma Correlation matrix.
#' @param M The number of knockoff copies (under the multiple knockoff framework).
#' @param n Sample size.
#' @param method The knockoff method used.
#' @param Gene_info Gene information.
#' @param verbose Indicator of whether progress of computation is plotted.
#' @param tol Convergence tolerance
#'
#' @return Knockoff q-values of all features (the smallest target FDR level for each feature to be selected).
#' @importFrom ghostbasil GhostMatrix
#' @importFrom ghostbasil ghostbasil
#' @importFrom stats rnorm
#' @importFrom stats qbeta
#' @importFrom Rfast spdinv
#' @export
#'
#' @examples
#' data("DGE_result")
#' data("Sigma")
#' set.seed(433)
#' Z<-Z_calculation(DGE_result$pvalue,DGE_result$log2FoldChange)
#' KnockoffCDGE(Z,Sigma,M=5,n=23010,method="ME",Gene_info=DGE_result)
KnockoffCDGE<-function(Z,Sigma,M,n,method=c("ME","SDP","SDP_no_julia","PCA"),Gene_info=NULL,verbose=TRUE,tol=0.001)
{
  is_null_gene_info<-is.null(Gene_info)
  p<-nrow(Sigma)
  if(method=="ME")
  {
    Inv_Sigma<-spdinv(Sigma)
    S<-S_calculation_ME(Sigma,M,verbose,tol)

    EZ_KO<-Z-colSums(Inv_Sigma*Z)*S
    V1<-diag((M+1)/M*S)-Inv_Sigma*(S%*%t(S))
    V1.right<-chol(V1)
    Z1<-colSums(V1.right*rnorm(p))
    Z2<-matrix(rnorm(M*p),p,M)
    Z2<-Z2-rowMeans(Z2)
    Z2<-Z2*sqrt(S)
    Z_KO<-Z2+EZ_KO+Z1

    Sigma_minus_S<-Sigma
    diag(Sigma_minus_S)<-diag(Sigma_minus_S)-S
    A<-GhostMatrix(Sigma = Sigma_minus_S,D=S,n.groups = M+1)
  }
  if(method=="SDP")
  {
    Inv_Sigma<-spdinv(Sigma)
    S<-S_calculation_SDP(Sigma,M,verbose,tol)

    EZ_KO<-Z-colSums(Inv_Sigma*Z)*S
    V1<-diag((M+1)/M*S)-Inv_Sigma*(S%*%t(S))
    V1.right<-chol(V1)
    Z1<-colSums(V1.right*rnorm(p))
    Z2<-matrix(rnorm(M*p),p,M)
    Z2<-Z2-rowMeans(Z2)
    Z2<-Z2*sqrt(S)
    Z_KO<-Z2+EZ_KO+Z1

    Sigma_minus_S<-Sigma
    diag(Sigma_minus_S)<-diag(Sigma_minus_S)-S
    A<-GhostMatrix(Sigma = Sigma_minus_S,D=S,n.groups = M+1)
  }
  if(method=="SDP_no_julia")
  {
    Inv_Sigma<-spdinv(Sigma)
    S<-S_calculation_SDP_no_julia(Sigma,M,verbose,tol)

    EZ_KO<-Z-colSums(Inv_Sigma*Z)*S
    V1<-diag((M+1)/M*S)-Inv_Sigma*(S%*%t(S))
    V1.right<-chol(V1)
    Z1<-colSums(V1.right*rnorm(p))
    Z2<-matrix(rnorm(M*p),p,M)
    Z2<-Z2-rowMeans(Z2)
    Z2<-Z2*sqrt(S)
    Z_KO<-Z2+EZ_KO+Z1

    Sigma_minus_S<-Sigma
    diag(Sigma_minus_S)<-diag(Sigma_minus_S)-S
    A<-GhostMatrix(Sigma = Sigma_minus_S,D=S,n.groups = M+1)
  }
  if(method=="PCA")
  {
    print(paste0("What is the number of principal components to be conditioned on? Please input integer from 1 to ",p,": "))
    while(TRUE)
    {
      r<-readline()
      if(r%in%as.character(1:p))
      {
        print("Valid input.")
        r<-as.integer(r)
        break
      }
      else
      {
        print(paste0("Invalid input. Please input integer from 1 to ",p,": "))
      }
    }
    N1_list<-S_calculation_PCA(Sigma,r,verbose,tol)
    S<-N1_list$S
    Z_KO<-((N1_list$eigenvectors)%*%(t(N1_list$eigenvectors)%*%Z))[,1]+matrix(rnorm(p*M),p,M)*sqrt(S)

    Sigma_minus_S<-N1_list$eigenvectors%*%(t(N1_list$eigenvectors)*N1_list$eigenvalues)
    A<-GhostMatrix(Sigma = Sigma_minus_S,D=S,n.groups = M+1)
  }

  # Run ghostbasil
  {
    rr<-c(Z,as.vector(Z_KO))/sqrt(n)
    gamma_max<-max(abs(rr))
    B<-1000
    q_b<-qbeta((1:B-0.5)/B,(M+1)*p,1)
    u_b<-qnorm((q_b+1)/2)
    gamma_min = 0.6*mean(u_b)/sqrt(n)
    if(gamma_max>gamma_min)
    {
      gamma_seq<-seq(gamma_max,gamma_min,(gamma_min-gamma_max)/500)
      fitted_model<-ghostbasil(A,rr,user.lambdas = gamma_seq)
      beta_lasso<-matrix(fitted_model$betas[,ncol(fitted_model$betas)],nrow = p,ncol = M+1)
    }else{
      beta_lasso<-matrix(0,nrow = p,ncol = M+1)
    }
  }

  #### Knockoffs Inference
  {
    kappa_ko<-apply(abs(beta_lasso),1,which.max)-1 # Magnitude of gene signal
    tau_ko<-apply(abs(beta_lasso),1,tau_calculation) # Indicator of whether the signal comes from the true gene or the knockoff copy
    q_value_ko<-KO_Filter(tau_ko,kappa_ko,M=M) # Compute the q-value of KO for CDGE analysis (q-value is of the same use as padj)
    if(is_null_gene_info)
    {
      Gene_info<-data.frame(Z=Z,S=S,q_value_ko=q_value_ko)
    }else
    {
      Gene_info<-as.data.frame(Gene_info)
      Gene_info$Z<-Z
      Gene_info$S<-S
      Gene_info$q_value_ko<-q_value_ko
    }
  }
  return(Gene_info)
}

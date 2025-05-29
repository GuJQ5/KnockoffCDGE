#' S-matrix computation (PCA)
#'
#' This function computes the S-matrix from correlation matrix under the PCA knockoff framework.
#'
#' @param Sigma Correlation matrix.
#' @param r The number of principal components conditioned on.
#' @param verbose Indicator of whether progress of computation is plotted.
#' @param tol Convergence tolerance
#'
#' @return A list of objects for knockoff generation.
#' @importFrom knockoffsr knockoff_setup
#' @export
S_calculation_PCA<-function(Sigma,r,verbose=TRUE,tol=0.001)
{
  p<-nrow(Sigma)
  NN1<-Sigma
  N_0<-0*Sigma
  diag(NN1)<-0

  eigs_k<-eigs(NN1,r)
  N_1<-(eigs_k$vectors)%*%(t(eigs_k$vectors)*eigs_k$values[1:r])
  while(mean(abs(N_1-N_0))>tol)
  {
    N_0<-N_1
    diag(NN1)<-diag(N_0)
    eigs_k<-eigs(NN1,r)
    N_1<-(eigs_k$vectors)%*%(t(eigs_k$vectors)*eigs_k$values[1:r])
  }
  S<-diag(Sigma)-diag(N_1)
  N_1_list<-list(eigenvalues=eigs_k$values[1:r],eigenvectors=eigs_k$vectors,S=S)
  return(N_1_list)
}

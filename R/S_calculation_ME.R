#' S-matrix computation (maximizing the entropy)
#'
#' This function computes the S-matrix from correlation matrix under the maximizing the entropy criterion.
#'
#' @param Sigma Correlation matrix.
#' @param M The number of knockoff copies (under the multiple knockoff framework).
#' @param verbose Indicator of whether progress of computation is plotted.
#' @param tol Convergence tolerance
#'
#' @return Diagnoal elements of the S-matrix for knockoff generation.
#' @importFrom knockoffsr knockoff_setup
#' @export
S_calculation_ME<-function(Sigma,M,verbose=TRUE,tol=0.001)
{
  ko<-knockoff_setup()
  p<-nrow(Sigma)
  solver <- "maxent"
  result <- ko$solve_s_graphical_group(Sigma, 1:p, 1:p, solver, m=M, verbose=verbose,tol=tol)
  S<-diag(result[[2]])
  return(S)
}

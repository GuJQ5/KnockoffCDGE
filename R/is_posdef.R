#' Positive semidefiniteness detection
#'
#' This function detects whether a input matrix is positive semidefinite or not
#'
#' @param A The input matrix
#'
#' @return An indicator whether the input matrixis positive semidefinite or not.
#' @export
#' @importFrom RSpectra eigs
#'
is_posdef = function(A) {
  A<-(A+t(A))/2
  if(nrow(A)<=2)
  {
    lambda_min<-min(eigen(A,only.values = TRUE)$values)
  }
  else
  {
    lambda_min<-eigs(A, 1, which = "SR")$values
  }
  return (lambda_min>1e-6)
}

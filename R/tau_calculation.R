#' Feature Statistics (Magnitude)
#'
#' This function computes feature statistics (magnitude) from importance scores of original feature and its knockoff copies.
#'
#' @param x Feature importance scores of original feature and its knockoff copies.
#'
#' @return Feature statistics (magnitude)
#' @export
#' @importFrom stats median
#'
#' @examples
#' tau_calculation(1:6)
tau_calculation<-function(x){
  l<-sort(x,decreasing = T)
  tau<-l[1]-median(l[-1])
  return(tau)
}

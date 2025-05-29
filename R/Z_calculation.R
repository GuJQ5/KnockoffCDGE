#' Z-scores computation
#'
#' This function computes Z-scores from DGE p-values (unadjusted) and log2 fold change.
#'
#' @param p_value DGE p-values.
#' @param log2FoldChange Log2 fold change.
#'
#' @return DGE Z-scores.
#' @export
#' @importFrom stats qnorm
#'
Z_calculation<-function(p_value,log2FoldChange)
{
  Z<-sign(log2FoldChange)*qnorm(1-p_value/2)
  return(Z)
}

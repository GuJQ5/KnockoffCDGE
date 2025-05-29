#' DGE results from psoriasis dataset.
#'
#' The DGE analysis results on 682 genes in 23010 T-cells of 3 psoriasis patients. The original records of cells are derived from the single-cell RNA-seq dataset collected by Reynolds et al. (2021) in the study of inflammatory skin diseases.
#'
#' @format A data frame with 682 rows and 6 variables:
#' \describe{
#'   \item{Gene}{The name of genes.}
#'   \item{baseMean}{Mean of normalized counts for all samples}
#'   \item{log2FoldChange}{Log2 fold change.}
#'   \item{lfcSE}{Standard error of log2 fold change.}
#'   \item{pvalue}{Wald test p-value.}
#'   \item{padj}{BH adjusted p-values.}
#' }
#'
#' @examples
#' data(DGE_result)
#' head(DGE_result)
"DGE_result"

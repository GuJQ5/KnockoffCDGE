% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sigma.R
\docType{data}
\name{Sigma}
\alias{Sigma}
\title{Correlation matrix from psoriasis dataset.}
\format{
A matrix with 682 rows and 682 column.
}
\usage{
Sigma
}
\description{
The estimated correlation matrix of expression levels of 682 genes in 23010 T-cells of 3 psoriasis patients. The original records of cells are derived from the single-cell RNA-seq dataset collected by Reynolds et al. (2021) in the study of inflammatory skin diseases.
}
\examples{
data(Sigma)
Sigma[1:6,1:6]
}
\keyword{datasets}

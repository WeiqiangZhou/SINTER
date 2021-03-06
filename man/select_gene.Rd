% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SINTER_functions.R
\name{select_gene}
\alias{select_gene}
\title{Select variable gene from single-cell RNA-seq data}
\usage{
select_gene(expr, RNA_train_all, filter_th = 0.1, var_th = 0.2,
  dataplot = FALSE)
}
\arguments{
\item{expr}{Gene expression matrix of scRNA-seq (normalized counts).}

\item{RNA_train_all}{Gene expression matrix of bulk RNA-seq from ENCODE.}

\item{filter_th}{Threshold for filtering out genes with no expression in the given percentage of cells.}

\item{var_th}{Threshold for filtering out genes with low variability in the given percentage of genes.}

\item{dataplot}{Plot the mean and variance of all genes and the selected genes are highlighted.}
}
\value{
\item{expr_select}{Log2-transformed gene expression matrix of scRNA-seq for the selected genes.}
 \item{RNA_train_select}{Gene expression matrix of bulk RNA-seq for the selected genes.}
}
\description{
This function is used for selecting variable gene from scRNA-seq data.
}
\examples{
\dontrun{
expr_data_select <- select_gene(expr_in,RNA_train_all,filter_th=0.1,var_th=0.2,dataplot=TRUE)
}
}
\keyword{gene}
\keyword{scRNA-seq}
\keyword{variable}

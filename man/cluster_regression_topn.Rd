% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SINTER_functions.R
\name{cluster_regression_topn}
\alias{cluster_regression_topn}
\title{Linear regression function based on sure independence screening}
\usage{
cluster_regression_topn(DNase_train, RNA_train_mean, RNA_test_mean, top_n)
}
\arguments{
\item{DNase_train}{ENCODE cluster features from DNase-seq data for building the regression model.}

\item{RNA_train_mean}{Gene cluster mean from ENCODE RNA-seq data for building the regression model.}

\item{RNA_test_mean}{Gene cluster mean from scRNA-seq data for making predictions.}

\item{top_n}{Number of predictors used in the regression model.}
}
\value{
\item{y_pre}{A vector of predicted ENCODE cluster features.}
}
\description{
This function is used for building linear regression models based on the top N predictors that are most correlated with reponse.
}
\keyword{prediction}

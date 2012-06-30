\docType{methods}
\name{getCorrelation}
\alias{getCorrelation}
\alias{getCorrelation,-methods}
\alias{getCorrelation,methylBase-method}
\alias{getCorrelation-method}
\title{get correlation between samples in methylBase object}
\usage{
  getCorrelation(.Object,method="pearson",plot=FALSE)
}
\arguments{
  \item{.Object}{a methylBase object}

  \item{method}{a character string indicating which
  correlation coefficient (or covariance) is to be computed
  (default:"pearson", other options are "kendall" and
  "spearman")}

  \item{plot}{scatterPlot if TRUE (default:FALSE)}
}
\value{
  a correlation matrix object and plot scatterPlot
}
\description{
  get correlation between samples in methylBase object
}


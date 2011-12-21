\docType{methods}
\name{get.methylDiff}
\alias{get.methylDiff}
\alias{get.methylDiff,methylDiff-method}
\title{gets differentially methylated regions/bases based on cutoffs}
\usage{
  get.methylDiff(.Object,difference=25,qvalue=0.01)
}
\arguments{
  \item{.Object}{a methylDiff object}

  \item{difference}{cutoff for absolute value of %
  methylation change between test and control (default:25)}

  \item{qvalue}{cutoff for qvalue of differential
  methylation statistic (default:0.01)}
}
\value{
  a methylDiff object containing the differential
  methylated locations satisfying the criteria
}
\description{
  gets differentially methylated regions/bases based on
  cutoffs
}


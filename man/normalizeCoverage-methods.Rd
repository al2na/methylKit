\docType{methods}
\name{normalizeCoverage}
\alias{normalizeCoverage}
\alias{normalizeCoverage,methylRawList-method}
\title{normalize read coverage between samples}
\arguments{
  \item{obj}{\code{methylRawList} object}

  \item{method}{a string "mean" or "median" which denotes
  median or mean should be used to calculate scaling
  factor. (Default:median)}
}
\value{
  a \code{methylRawList} object
}
\description{
  The function normalizes coverage values between samples
  using a scaling factor derived from differences between
  mean or median of coverage distributions
}
\examples{
library(methylKit)
data(methylKit)
newObj=normalizeCoverage(methylRawList.obj)
}
\author{
  Altuna Akalin
}


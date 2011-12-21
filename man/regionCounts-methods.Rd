\docType{methods}
\name{regionCounts}
\alias{regionCounts}
\alias{regionCounts,methylRaw,GRanges-method}
\alias{regionCounts,methylRaw,GRangesList-method}
\alias{regionCounts,methylRawList,GRanges-method}
\alias{regionCounts,methylRawList,GRangesList-method}
\title{GETs regional counts for given GRanges or GRangesList object}
\usage{
  regionCounts(methylObj,regions,cov.bases=0)
}
\arguments{
  \item{methylObj}{a \code{methylRaw} or
  \code{methlRawList} object}

  \item{regions}{a GRanges or GRangesList object.}

  \item{cov.bases}{number minimum bases covered per region
  (Default:0). Only regions with base coverage above this
  threshold are returned.}
}
\value{
  RETURNS a new methylRaw or methylRawList object
}
\description{
  Convert \code{methylRaw} or \code{methylRawList} object
  into regional counts for a given GRanges or GRangesList
  object.
}


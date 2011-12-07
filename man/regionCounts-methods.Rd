\docType{methods}
\name{regionCounts}
\alias{regionCounts}
\alias{regionCounts,-methods}
\alias{regionCounts,methylRaw,methyRawList,ANY-method}
\alias{regionCounts,methylRawList-method}
\title{GETs regional counts for given GRanges or GRangesList object}
\arguments{
  \item{methylObj}{a \code{methylRaw} or
  \code{methlRawList} object}

  \item{regions}{a GRanges or GRangesList object.}
}
\value{
  RETURNS a new methylRaw or methylRawList object
}
\description{
  Convert CpG dinucleotide methylation \code{methylRaw} or
  \code{methylRawList} object into regional counts for a
  given GRanges or GRangesList object.
}


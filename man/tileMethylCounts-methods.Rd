\docType{methods}
\name{tileMethylCounts}
\alias{tileMethylCounts}
\alias{tileMethylCounts,methylRaw-method}
\alias{tileMethylCounts,methylRawList-method}
\title{Get methylated/unmethylated base counts for tilling windows}
\arguments{
  \item{methylObj}{\code{methylRaw} or \code{methylRawList}
  object containing base pair resolution methylation
  information}

  \item{win.size}{an integer for the size of the tiling
  windows}

  \item{step.size}{an integer for the step size of tiling
  windows}

  \item{cov.bases}{minimum number of bases to be covered in
  a given window}
}
\value{
  \code{methylRaw} or \code{methylRawList} object
}
\description{
  The function summarizes methylated/unmethylated base
  counts over tilling windows accross genome. This function
  can be used when differential methylated analysis is
  preferable to tilling windows instead of base pairs.
}


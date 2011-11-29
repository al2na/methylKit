\docType{methods}
\name{unite}
\alias{unite}
\alias{unite,-methods}
\alias{unite,ANY-method}
\alias{unite,methylRawList-method}
\title{unites methylRawList to a single table}
\arguments{
  \item{destrand}{if TRUE, reads covering both strands of a
  CpG dinucleotide will be merged, do not set to TRUE if
  not only interested in CpGs (default: FALSE)}
}
\value{
  a methylBase object
}
\description{
  This functions unites \code{methylRawList} object that
  only bases with coverage from all samples are retained.
  The resulting object is a class of \code{methylBase}
}


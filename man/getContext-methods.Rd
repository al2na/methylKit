\docType{methods}
\name{getContext}
\alias{getContext}
\alias{getContext,methylBase-method}
\alias{getContext,methylDiff-method}
\alias{getContext,methylRaw-method}
\title{get the context of methylation}
\usage{
  getContext(x)
}
\arguments{
  \item{x}{an
  \code{\link{methylBase}},\code{\link{methylRaw}} or an
  \code{\link{methylDiff}} object}
}
\value{
  a string for the context methylation
}
\description{
  The function returns the context of methylation. For
  example: "CpG","CHH" or "CHG"
}
\examples{
data(methylKit)

getContext(methylBase.obj)
getContext(methylDiff.obj)
getContext(methylRawList.obj[[1]])
}


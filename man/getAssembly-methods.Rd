\docType{methods}
\name{getAssembly}
\alias{getAssembly}
\alias{getAssembly,methylBase-method}
\alias{getAssembly,methylDiff-method}
\alias{getAssembly,methylRaw-method}
\title{get assembly of the genome}
\usage{
  getAssembly(x)
}
\arguments{
  \item{x}{an
  \code{\link{methylBase}},\code{\link{methylRaw}} or
  \code{\link{methylDiff}} object}
}
\value{
  the assembly string for the object
}
\description{
  The function returns the genome assembly stored in any of
  the
  \code{\link{methylBase}},\code{\link{methylRaw}},\code{\link{methylDiff}}
  objects
}
\examples{
data(methylKit)

getAssembly(methylBase.obj)
getAssembly(methylDiff.obj)
getAssembly(methylRawList.obj[[1]])
}


\docType{methods}
\name{extract}
\alias{[}
\alias{[,methylBase-method}
\alias{[,methylDiff-method}
\alias{[,methylRaw-method}
\alias{extract}
\title{extract parts of methylRaw,methylBase and methylDiff data}
\arguments{
  \item{x}{an
  \code{\link{methylBase}},\code{\link{methylRaw}} or
  \code{\link{methylDiff}} object}

  \item{i}{a numeric or logical vector. This vector
  corresponds to bases or regions contained in
  \code{methylKit} objects.The vector is used to subset the
  data.}

  \item{j}{This argument can not be used for the extraction
  of columns. As unintentional extraction of the columns
  will cause an error in the downstream analysis. Using
  this argument will cause an error. Use
  \code{\link[methylKit]{getData}} to access the data part
  of the objects.}
}
\description{
  The function extracts part of the data and returns a new
  object.
}
\examples{
data(methylKit)

# selects first hundred rows, returns a methylRaw object
subset1=methylRawList.obj[[1]][1:100]

# selects first hundred rows, returns a methylBase object
subset2=methylBase.obj[1:100,]

# selects first hundred rows, returns a methylDiff object
subset3=methylDiff.obj[1:100,]

# This will get chromomsomes, will return a factor
# That means the resulting object will ceases to be a methylKit object
chrs=methylDiff.obj[[2]]
}


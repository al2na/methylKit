\docType{methods}
\name{select}
\alias{select}
\alias{select,methylBase-method}
\alias{select,methylDiff-method}
\alias{select,methylRaw-method}
\title{selects rows from of methylKit objects}
\usage{
  select(x,i)
}
\arguments{
  \item{x}{an
  \code{\link{methylBase}},\code{\link{methylRaw}} or
  \code{\link{methylDiff}} object}

  \item{i}{a numeric or logical vector. This vector
  corresponds to bases or regions contained in
  \code{methylKit} objects.The vector is used to subset the
  data.}
}
\value{
  a \code{\link{methylBase}},\code{\link{methylRaw}} or
  \code{\link{methylDiff}} object depending on the input
  object.
}
\description{
  The function returns a subset of data contained in the
  \code{methylKit} objects.
}
\examples{
data(methylKit)
 # selects first hundred rows, returns a methylRaw object
subset1=select(methylRawList.obj[[1]],1:100)

# selects first hundred rows, returns a methylBase object
subset2=select(methylBase.obj,1:100)

# selects first hundred rows, returns a methylDiff object
subset3=select(methylDiff.obj,1:100)
}


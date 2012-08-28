\docType{methods}
\name{getData}
\alias{getData}
\alias{getData-method}
\alias{getData,-methods}
\alias{getData,methylBase-method}
\alias{getData,methylDiff-method}
\alias{getData,methylRaw-method}
\title{get the data slot from the methylKit objects}
\usage{
  getData(x)
}
\arguments{
  \item{x}{an
  \code{\link{methylBase}},\code{\link{methylRaw}} or
  \code{\link{methylDiff}} object}
}
\value{
  data.frame for methylation events
}
\description{
  The functions retrieves the table containing methylation
  information from \code{methylKit} Objects. The data
  retrived from this function is of a
  \code{\link{data.frame}}. This is basically containing
  all relevant methylation information per genomic region
  or base.
}
\examples{
data(methylKit)

# following commands show first lines of returned data.frames from getData() function
head(
getData(methylBase.obj)
)

head( getData(methylDiff.obj))

head(getData(methylRawList.obj[[1]]))
}


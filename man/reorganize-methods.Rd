\docType{methods}
\name{reorganize}
\alias{reorganize}
\alias{reorganize,methylBase-method}
\alias{reorganize,methylRawList-method}
\title{reorganize methylRawList and methylBase objects by creating new objects from subset of samples}
\usage{
  reorganize(methylObj,sample.ids,treatment)
}
\arguments{
  \item{methylObj}{a \code{methylRawList} or
  \code{methylBase} object}

  \item{sample.ids}{a vector for sample.ids to be subset.
  Order is important and the order should be similar to
  treatment. sample.ids should be a subset or reordered
  version of sample ids in the input object.}

  \item{treatment}{treatment vector, should be same length
  as sample.ids vector}
}
\value{
  RETURNS a \code{methylRawList} or \code{methylBase}
  object depending on the input object
}
\description{
  Create a new \code{methylRawList} or \code{methylBase}
  object by selecting a subset of samples from the input
  object, which is a \code{methylRawList} or
  \code{methylBase} object. You can use the function to
  partition a large methylRawList or methylBase object to
  smaller object based on sample ids or when you want to
  reorder samples and treatmet vector.
}


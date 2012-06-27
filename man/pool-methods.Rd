\docType{methods}
\name{pool}
\alias{pool}
\alias{pool,methylBase-method}
\title{function pools replicates within groups to a single sample per group}
\usage{
  pool(obj,sample.ids)
}
\arguments{
  \item{obj}{\code{methylBase} object with two groups or
  more and each group should have multiple samples}

  \item{sample.ids}{a character vector of new sample.ids
  ex:c("test","control"), should follow the same order as
  unique treatment vector, and should be equal to the
  length of the unique treatment vector}
}
\value{
  a \code{methylBase} object
}
\description{
  The function sums up coverage, numCs and numTs values
  within each group so one representative sample for each
  group will be created in a new methylBase object
}
\examples{
library(methylKit)
data(methylKit)
newBase=pool(methylBase.obj,sample.ids=c("test","control"))
}
\author{
  Altuna Akalin
}


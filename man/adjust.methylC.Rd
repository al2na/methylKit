\docType{methods}
\name{adjust.methylC}
\alias{adjust.methylC}
\alias{adjust.methylC,methylRaw,methylRaw-method}
\alias{adjust.methylC,methylRawList,methylRawList-method}
\title{adjust measured 5mC levels using 5hmC levels}
\usage{
  adjust.methylC(mc,hmc)
}
\arguments{
  \item{mc}{a \code{methylRawList} or \code{methylRaw}
  containing 5mC levels of a sample or set of samples}

  \item{hmc}{a \code{methylRawList} or \code{methylRaw}
  containing 5hmC levels of a sample or set of samples. If
  a \code{methylRawList} givenm the sample order should be
  same as "mc" \code{methylRawList} object.}
}
\value{
  returns adjusted 5-methyl cytosine levels in the form of
  \code{methylRawList} or \code{methylRaw} object depending
  on the input object
}
\description{
  Measured 5mC levels via bisulfite sequencing might be a
  combination of 5hmC and 5mC levels since bisulfite
  sequencing can not distinguish between the two. This
  function can adjust 5mC levels of a bisulfite sequencing
  experiment if the user supplies corresponding 5hmC levels
  from the same sample.
}


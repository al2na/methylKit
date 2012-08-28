\docType{methods}
\name{unite}
\alias{unite}
\alias{unite,-methods}
\alias{unite,methylRawList-method}
\title{unites methylRawList to a single table}
\usage{
  unite(.Object,destrand=FALSE,min.per.group=NULL)
}
\arguments{
  \item{.Object}{a methylRawList object to be merged by
  common locations covered by reads}

  \item{destrand}{if TRUE, reads covering both strands of a
  CpG dinucleotide will be merged, do not set to TRUE if
  not only interested in CpGs (default: FALSE). If the
  methylRawList object contains regions rather than bases
  setting destrand to TRUE will have no effect.}

  \item{min.per.group}{an integer denoting minimum number
  of samples per replicate needed to cover a region/base.
  By default only regions/bases that are covered in all
  samples are united as methylBase object, however by
  supplying an integer for this argument users can control
  how many samples needed to cover region/base to be united
  as methylBase object. For example, if min.per.group set
  to 2 and there are 3 replicates per condition, the
  bases/regions that are covered in at least 2 replicates
  will be united and missing data for uncovered
  bases/regions will appear as NAs.}
}
\value{
  a methylBase object
}
\description{
  This functions unites \code{methylRawList} object that
  only bases with coverage from all samples are retained.
  The resulting object is a class of \code{methylBase}
}
\examples{
data(methylKit)
## Following
my.methylBase=unite(methylRawList.obj)
my.methylBase=unite(methylRawList.obj,destrand=TRUE)
}


\docType{methods}
\name{regionCounts}
\alias{regionCounts}
\alias{regionCounts,methylBase,GRanges-method}
\alias{regionCounts,methylRaw,GRanges-method}
\alias{regionCounts,methylRaw,GRangesList-method}
\alias{regionCounts,methylRawList,GRanges-method}
\alias{regionCounts,methylRawList,GRangesList-method}
\title{GETs regional counts for given GRanges or GRangesList object}
\usage{
  regionCounts(object,regions,cov.bases=0,strand.aware=FALSE)
}
\arguments{
  \item{object}{a \code{methylRaw} or \code{methlRawList}
  object}

  \item{regions}{a GRanges or GRangesList object.}

  \item{cov.bases}{number minimum bases covered per region
  (Default:0). Only regions with base coverage above this
  threshold are returned.}

  \item{strand.aware}{if set to TRUE only CpGs that match
  the strand of the region will be summarized.
  (default:FALSE)}
}
\value{
  a new methylRaw,methylBase or methylRawList object. If
  \code{strand.aware} is set to FALSE (default). Even
  though the resulting object will have the strand
  information of \code{regions} it will still contain
  methylation information from both strands.
}
\description{
  Convert \code{\link{methylRaw}},
  \code{\link{methylRawList}} or \code{\link{methylBase}}
  object into regional counts for a given
  \code{\link{GRanges}} or \code{\link{GRangesList}}
  object.
}
\examples{
data(methylKit)

# get the windows of interest as a GRanges object, this can be any set
# of genomic locations
library(GenomicRanges)
my.win=GRanges(seqnames="chr21",
ranges=IRanges(start=seq(from=9764513,by=10000,length.out=20),width=5000) )

# getting counts per region
regional.methylRaw=regionCounts(object=methylRawList.obj, regions=my.win,
cov.bases=0,strand.aware=FALSE)
}


\docType{methods}
\name{read.transcript.features}
\alias{read.transcript.features}
\alias{read.transcript.features,character-method}
\title{function reading exon intron, promoter structure from a given bed file}
\usage{
  read.transcript.features(location,remove.unsual=TRUE,up.flank=1000,down.flank=1000,unique.prom=TRUE)
}
\arguments{
  \item{location}{location of the bed file with 12 or more
  columns}

  \item{remove.unsual}{remove the chromomesomes with unsual
  names, mainly random chromsomes etc}

  \item{up.flank}{up-stream from TSS to detect promoter
  boundaries}

  \item{down.flank}{down-stream from TSS to detect promoter
  boundaries}

  \item{unique.prom}{get only the unique promoters,
  promoter boundaries will not have a gene name if you set
  this option to be TRUE}
}
\value{
  a \code{\link{GRangesList}} containing locations of
  exon/intron/promoter/TSS
}
\description{
  function reading exon intron, promoter structure from a
  given bed file
}
\note{
  one bed track per file is only accepted, the bed files
  with multiple tracks will cause en error
}


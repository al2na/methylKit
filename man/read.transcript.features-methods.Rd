\docType{methods}
\name{read.transcript.features}
\alias{read.transcript.features}
\alias{read.transcript.features,character-method}
\title{Read  transcript features from a BED file}
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
  a \code{\link[GenomicRanges]{GRangesList}} containing
  locations of exon/intron/promoter/TSS
}
\description{
  The function returns a
  \code{\link[GenomicRanges]{GRangesList}} containing exon,
  intron, TSS(transcription start site) and promoter
  locations
}
\note{
  one bed track per file is only accepted, the bed files
  with multiple tracks will cause en error
}
\examples{
gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))
}


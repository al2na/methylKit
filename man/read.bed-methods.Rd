\docType{methods}
\name{read.bed}
\alias{read.bed}
\alias{read.bed,character-method}
\title{read a bed file and convert it to GRanges}
\usage{
  read.bed(location,remove.unsual=TRUE)
}
\arguments{
  \item{location}{location of the file, a character string
  such as: "/home/user/my.bed"}

  \item{remove.unsual}{if TRUE(default) remove the
  chromomesomes with unsual names, mainly random chromsomes
  etc}
}
\value{
  \code{\link{GRanges}} object
}
\description{
  The function reads a BED file and coverts it to a
  \code{\link[GenomicRanges]{GRanges}} object
}
\note{
  one bed track per file is only accepted, the bed files
  with multiple tracks will cause an error
}
\examples{
bed.file=system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit")
bed.gr=read.bed(location=bed.file,remove.unsual=TRUE)
}


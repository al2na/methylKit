\docType{methods}
\name{convert.bed2introns}
\alias{convert.bed2introns}
\alias{convert.bed2introns,data.frame-method}
\title{convert a data frame read-in from a bed file to a GRanges object for introns}
\usage{
  convert.bed2introns(bed.df)
}
\arguments{
  \item{bed.df}{a data.frame where column order and content
  resembles a bed file with 12 columns}
}
\value{
  \code{\link{GRanges}} object
}
\description{
  convert a data frame read-in from a bed file to a GRanges
  object for introns
}
\note{
  one bed track per file is only accepted, the bed files
  with multiple tracks will cause en error
}


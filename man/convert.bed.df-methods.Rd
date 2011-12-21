\docType{methods}
\name{convert.bed.df}
\alias{convert.bed.df}
\alias{convert.bed.df,data.frame-method}
\title{convert a data frame read-in from a bed file to a GRanges object}
\usage{
  convert.bed.df(bed)
}
\arguments{
  \item{bed}{a data.frame where column order and content
  resembles a bed file with 12 columns}
}
\value{
  \code{\link{GRanges}} object
}
\description{
  convert a data frame read-in from a bed file to a GRanges
  object
}
\note{
  one bed track per file is only accepted, the bed files
  with multiple tracks will cause en error
}


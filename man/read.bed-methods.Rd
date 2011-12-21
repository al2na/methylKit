\docType{methods}
\name{read.bed}
\alias{read.bed}
\alias{read.bed,character-method}
\title{read a bed file and convert it to GRanges}
\usage{
  read.bed(location,remove.unsual=T)
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
  read a bed file and convert it to GRanges
}
\note{
  one bed track per file is only accepted, the bed files
  with multiple tracks will cause en error
}


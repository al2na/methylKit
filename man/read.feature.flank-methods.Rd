\docType{methods}
\name{read.feature.flank}
\alias{read.feature.flank}
\alias{read.feature.flank,character-method}
\title{a function to read-in genomic features and their upstream and downstream adjecent regions such as CpG islands and their shores}
\usage{
  read.feature.flank(location,remove.unsual=T,flank=2000,clean=T,feature.flank.name=NULL)
}
\arguments{
  \item{location}{for the bed file of the feature}

  \item{flank}{number of basepairs for the flanking
  regions}

  \item{clean}{If set to TRUE, flanks overlapping with
  other main features will be trimmed}

  \item{remove.unsual}{remove chromsomes with unsual names
  random, Un and antyhing with "_" character}

  \item{feature.flank.name}{the names for feature and flank
  ranges, it should be a character vector of length 2.
  example: c("CpGi","shores")}
}
\value{
  a GenomicRangesList contatining one GRanges object for
  flanks and one for GRanges object for the main feature.
  NOTE:This can not return a GRangesList at the moment
  because flanking regions do not have to have the same
  column name as the feature. GRangesList elements should
  resemble eachother in the column content. We can not
  satisfy that criteria for the flanks
}
\description{
  a function to read-in genomic features and their upstream
  and downstream adjecent regions such as CpG islands and
  their shores
}


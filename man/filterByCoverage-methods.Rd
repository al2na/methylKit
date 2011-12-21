\docType{methods}
\name{filterByCoverage}
\alias{filterByCoverage}
\alias{filterByCoverage,methylRaw-method}
\alias{filterByCoverage,methylRawList-method}
\title{filter methylRaw and methylRawList object based on read coverage}
\usage{
  filterByCoverage(methylObj,lo.count=NULL,lo.perc=NULL,hi.count=NULL,hi.perc=NULL)
}
\arguments{
  \item{methylObj}{a \code{methylRaw} or
  \code{methylRawList} object}

  \item{lo.count}{An integer for read counts.Bases/regions
  having lower coverage than this count is discarded}

  \item{lo.perc}{A double [0-100] for percentile of read
  counts. Bases/regions having lower coverage than this
  percentile is discarded}

  \item{hi.count}{An integer for read counts. Bases/regions
  having higher coverage than this is count discarded}

  \item{hi.perc}{A double [0-100] for percentile of read
  counts. Bases/regions having higher coverage than this
  percentile is discarded}
}
\value{
  \code{methylRaw} or \code{methylRawList} object depending
  on input object
}
\description{
  This function filters \code{methylRaw} and
  \code{methylRawList} objects. You can filter based on
  lower read cutoff or high read cutoff. Higher read cutoff
  is usefull to eliminate PCR effects Lower read cutoff is
  usefull for doing better statistical tests.
}


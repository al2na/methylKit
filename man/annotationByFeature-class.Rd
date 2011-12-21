\name{annotationByFeature-class}
\alias{annotationByFeature-class}
\title{An S4 class that information on overlap of target features with annotation features}
\description{
  This object is desgined to hold statistics and
  information about genomic feature overlaps
}
\section{Slots}{
  \describe{ \item{members}{a matrix showing overlap of
  target features with annotation genomic features}

  \item{annotation}{a named vector of percentages}

  \item{precedence}{a named vector of percentages}

  \item{num.hierarchica}{vector}

  \item{no.of.OlapFeat}{vector}

  \item{perc.of.OlapFeat}{vector} }
}


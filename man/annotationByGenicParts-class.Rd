\name{annotationByGenicParts-class}
\alias{annotationByGenicParts-class}
\title{An S4 class that information on overlap of target features with annotation features}
\description{
  This object is desgined to hold statistics and
  information about genomic feature overlaps
}
\section{Slots}{
  \describe{ \item{\code{members}}{a matrix showing overlap
  of target features with annotation genomic features}

  \item{\code{annotation}}{a named vector of percentages}

  \item{\code{precedence}}{a named vector of percentages}

  \item{\code{num.hierarchica}}{vector}

  \item{\code{no.of.OlapFeat}}{vector}

  \item{\code{perc.of.OlapFeat}}{vector}

  \item{dist.to.TSS}{a data frame showing distances to TSS
  and gene/TSS names and strand} }
}


\docType{methods}
\name{getAssociationWithTSS}
\alias{getAssociationWithTSS}
\alias{getAssociationWithTSS,-methods}
\alias{getAssociationWithTSS,annotationByGenicParts-method}
\title{Get distance to nearest TSS and gene id from annotationByGenicParts}
\usage{
  getAssociationWithTSS(x)
}
\arguments{
  \item{x}{a \code{annotationByGenicParts} object}
}
\value{
  RETURNS a data.frame containing row number of the target
  features,distance of target to nearest TSS, TSS/Gene
  name, TSS strand
}
\description{
  This accessor function gets the nearest TSS, its distance
  to target feature,strand and name of TSS/gene from
  annotationByGenicParts object
}


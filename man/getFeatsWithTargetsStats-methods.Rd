\docType{methods}
\name{getFeatsWithTargetsStats}
\alias{getFeatsWithTargetsStats}
\alias{getFeatsWithTargetsStats,annotationByFeature-method}
\title{Get the percentage/count of annotation features overlapping with target features from annotationByFeature}
\arguments{
  \item{x}{a \code{annotationByFeature} object}

  \item{percentage}{TRUE|FALSE. If TRUE percentage of
  annotation features will be returned. If FALSE, number of
  annotation features will be returned}
}
\value{
  RETURNS a vector of percentages or counts showing
  quantity of annotation features overlapping with target
  features
}
\description{
  This function retrieves percentage/number of annotation
  features overlapping with targets. For example, if
  \code{annotationByFeature} object is containing
  statistics of differentially methylated regions
  overlapping with gene annotation. This function will
  return number/percentage of introns,exons and promoters
  overlapping with differentially methylated regions.
}


\docType{methods}
\name{getTargetAnnotationStats}
\alias{getTargetAnnotationStats}
\alias{getTargetAnnotationStats,annotationByFeature-method}
\title{Get the percentage of target features overlapping with annotation from annotationByFeature}
\arguments{
  \item{x}{a \code{annotationByFeature} object}

  \item{percentage}{TRUE|FALSE. If TRUE percentage of
  target features will be returned. If FALSE, number of
  target features will be returned}

  \item{precedence}{TRUE|FALSE. If TRUE there will be a
  hierachy of annotation features when calculating numbers
  (with promoter>exon>intron precedence) That means if a
  feature overlaps with a promoter it will be counted as
  promoter overlapping only, or if it is overlapping with a
  an exon but not a promoter, it will be counted as exon
  overlapping only whether or not it overlaps with an
  intron.}
}
\value{
  RETURNS a vector of percentages or counts showing
  quantity of target features overlapping with annotation
}
\description{
  This function retrieves percentage/number of target
  features overlapping with annotation
}


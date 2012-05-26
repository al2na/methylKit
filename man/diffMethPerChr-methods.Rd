\docType{methods}
\name{diffMethPerChr}
\alias{diffMethPerChr}
\alias{diffMethPerChr,methylDiff-method}
\title{Gets and plots the number of hyper/hypo methylated regions per chromosome}
\arguments{
  \item{x}{a \code{annotationByFeature} object}

  \item{plot}{TRUE|FALSE. If TRUE horizontal barplots for
  proportion of hypo/hyper methylated bases/regions}

  \item{qvalue.cutoff}{cutoff for q-value}

  \item{meth.cutoff}{cutoff for percent methylation
  difference}

  \item{exclude}{names of chromosomes to be excluded}

  \item{...}{extra graphical parameters to be passed to
  \code{barplot} function}
}
\value{
  plots a piechart or a barplot for percentage of the
  target features overlapping with annotation
}
\description{
  This accessor function gets the nearest TSS, its distance
  to target feature,strand and name of TSS/gene from
  annotationByGenicParts object.
}


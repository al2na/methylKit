\docType{methods}
\name{plotTargetAnnotation}
\alias{plotTargetAnnotation}
\alias{plotTargetAnnotation,annotationByFeature-method}
\title{Plot annotation categories from annotationByGenicParts or annotationByFeature}
\usage{
  plotTargetAnnotation(x,precedence=TRUE,col=rainbow(length(x@annotation)),
    ...)
}
\arguments{
  \item{x}{an \code{\link[methylKit]{annotationByFeature}}
  or an \code{\link[methylKit]{annotationByGenicParts}}
  object}

  \item{precedence}{TRUE|FALSE. If TRUE there will be a
  hierachy of annotation features when calculating numbers
  (with promoter>exon>intron precedence).  This option is
  only valid when x is a \code{annotationByGenicParts}
  object}

  \item{col}{a vector of colors for piechart or the bar
  plot}

  \item{...}{graphical parameters to be passed to
  \code{pie} or \code{barplot} functions}
}
\value{
  plots a piechart or a barplot for percentage of the
  target features overlapping with annotation
}
\description{
  This function plots a pie or bar chart for showing
  percentages of targets annotated by genic parts or other
  query features.
}
\examples{
data(methylKit)
gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))
ann=annotate.WithGenicParts(methylDiff.obj, gene.obj)
plotTargetAnnotation(ann,precedence=FALSE)
}


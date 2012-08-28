\docType{methods}
\name{getAssociationWithTSS}
\alias{getAssociationWithTSS}
\alias{getAssociationWithTSS,annotationByGenicParts-method}
\title{Get distance to nearest TSS and gene id from annotationByGenicParts}
\usage{
  getAssociationWithTSS(x)
}
\arguments{
  \item{x}{an
  \code{\link[methylKit]{annotationByGenicParts}} object}
}
\value{
  a \code{\link{data.frame}} containing row number of the
  target features,distance of target to nearest TSS,
  TSS/Gene name, TSS strand
}
\description{
  This accessor function gets the nearest TSS, its distance
  to target feature,strand and name of TSS/gene from
  annotationByGenicParts object
}
\examples{
data(methylKit)
gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))
ann=annotate.WithGenicParts(methylDiff.obj, gene.obj)
df=getAssociationWithTSS(ann)
head(df)
}


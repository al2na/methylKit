\docType{methods}
\name{getMembers}
\alias{getMembers}
\alias{getMembers,annotationByFeature-method}
\title{Get the membership slot of annotationByFeature}
\usage{
  getMembers(x)
}
\arguments{
  \item{x}{an \code{\link[methylKit]{annotationByFeature}}
  object}
}
\value{
  a matrix showing overlap of target features with
  annotation features. 1 for overlap, 0 for non-overlap.
  Each row in the matrix corresponds to a genomic feature
  that is annoted by one of the following functions:
  \code{\link[methylKit]{annotate.WithFeature}},
  \code{\link[methylKit]{annotate.WithFeature.Flank}},
  \code{\link[methylKit]{annotate.WithGenicParts}}
}
\description{
  Membership slot defines the overlap of target features
  with annotation features as a matrix.
}
\examples{
data(methylKit)
cpg.obj=read.feature.flank(system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit"),feature.flank.name=c("CpGi","shores"))

ann=annotate.WithFeature.Flank(methylDiff.obj,cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="Shores")
mat=getMembers(ann)
head(mat)
}


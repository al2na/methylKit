\docType{methods}
\name{annotate.WithFeature}
\alias{annotate.WithFeature}
\alias{annotate.WithFeature,GRanges,GRanges-method}
\alias{annotate.WithFeature,methylDiff,GRanges-method}
\title{Annotate object with a set of genomic features}
\usage{
  annotate.WithFeature(target,feature,strand=FALSE,extend=0,feature.name="feat1")
}
\arguments{
  \item{target}{a \code{\link[methylKit]{methylDiff}} or a
  \code{\link[GenomicRanges]{GRanges}} object storing
  chromosome locations to be annotated}

  \item{feature}{a \code{\link[GenomicRanges]{GRanges}}
  object storing chromosome locations of a feature (can be
  CpG islands, ChIP-seq peaks, etc)}

  \item{strand}{If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FALSE)}

  \item{extend}{specifiying a positive value will extend
  the feature on both sides as much as \code{extend}}

  \item{feature.name}{name of the annotation feature. For
  example: H3K4me1,CpGisland etc.}
}
\value{
  returns an \code{\link[methylKit]{annotationByFeature}}
  object
}
\description{
  The function annotates given genomic feature or methylKit
  object with a set of annotation. It is primarily useful
  when annotating objects with simple genomic features,
  such as enhancer locations.
}
\examples{
data(methylKit)
cpg.gr=read.bed(system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit"),remove.unsual=TRUE)

annotate.WithFeature(methylDiff.obj,cpg.gr,strand=FALSE,extend=0,feature.name="CpGi")
}
\seealso{
  \code{\link[methylKit]{getMembers}},
  \code{\link[methylKit]{getTargetAnnotationStats}},
  \code{\link[methylKit]{getFeatsWithTargetsStats}},
  \code{\link[methylKit]{plotTargetAnnotation}}
}


\docType{methods}
\name{annotate.WithFeature.Flank}
\alias{annotate.WithFeature.Flank}
\alias{annotate.WithFeature.Flank,GRanges,GRanges,GRanges-method}
\alias{annotate.WithFeature.Flank,methylDiff,GRanges,GRanges-method}
\title{function to annotate given GRanges object with promoter,exon,intron & intergenic values}
\usage{
  annotate.WithFeature.Flank(target,feature,flank,feature.name="feat",flank.name="flank",strand=FALSE)
}
\arguments{
  \item{target}{a methylDiff or a granges object storing
  chromosome locations to be annotated}

  \item{feature}{a granges object storing chromosome
  locations of a feature (can be CpG islands, ChIP-seq
  peaks, etc)}

  \item{flank}{a granges object storing chromosome
  locations of the flanks of the feature}

  \item{feature.name}{string for the name of the feature}

  \item{flank.name}{string for the name o}

  \item{strand}{If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FALSE)}
}
\value{
  returns an \code{annotationByFeature} object
}
\description{
  function to annotate given GRanges object with
  promoter,exon,intron & intergenic values
}
\examples{
data(methylKit)
cpg.obj=read.feature.flank(system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit"),feature.flank.name=c("CpGi","shores"))

annotate.WithFeature.Flank(methylDiff.obj,cpg.obj$CpGi,cpg.obj$shores)
}


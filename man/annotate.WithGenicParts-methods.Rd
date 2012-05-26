\docType{methods}
\name{annotate.WithGenicParts}
\alias{annotate.WithGenicParts}
\alias{annotate.WithGenicParts,GRanges,GRangesList-method}
\alias{annotate.WithGenicParts,methylDiff,GRangesList-method}
\title{function to annotate given GRanges object with promoter,exon,intron & intergenic values}
\arguments{
  \item{target}{: a methylDiff or a granges object storing
  chromosome locations to be annotated}

  \item{GRangesList.obj}{: A GRangesList object containing
  GRanges object for promoter,exons,introns and TSSes, or
  simply output of read.transcript.features function}

  \item{strand}{: If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FALSE)}
}
\value{
  \code{annotationByGenicParts} object
}
\description{
  function to annotate given GRanges object with
  promoter,exon,intron & intergenic values
}
\examples{
data(methylKit)
gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))
annotate.WithGenicParts(methylDiff.obj,gene.obj)
}


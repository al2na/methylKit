\docType{methods}
\name{annotate.WithGenicParts}
\alias{annotate.WithGenicParts}
\alias{annotate.WithGenicParts,GRanges,GRangesList-method}
\alias{annotate.WithGenicParts,methylDiff,GRangesList-method}
\title{function to annotate given GRanges object with promoter,exon,intron & intergenic values}
\usage{
  annotate.WithGenicParts(target,GRangesList.obj,strand=F)
}
\arguments{
  \item{target}{: a granges object storing chromosome
  locations to be annotated}

  \item{GRangesList.obj}{: A GRangesList object containing
  GRanges object for promoter,exons,introns and TSSes, or
  simply output of read.transcript.features function}

  \item{strand}{: If set to TRUE, annotation features and
  target features will be overlapped based on strand
  (def:FAULT)}
}
\value{
  \code{annotationByGenicParts} object
}
\description{
  function to annotate given GRanges object with
  promoter,exon,intron & intergenic values
}


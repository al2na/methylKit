\name{methylRaw-class}
\alias{methylRaw-class}
\title{An S4 class for holding raw methylation data from alignment pipeline.}
\description{
  This object stores the raw mehylation data that is read
  in through read function and extends data.frame
}
\section{Slots}{
  \describe{ \item{\code{sample.id}:}{string for an
  identifier of the sample} \item{\code{assembly}:}{string
  for genome assembly, ex: hg18,hg19,mm9} }
}


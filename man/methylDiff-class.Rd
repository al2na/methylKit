\docType{class}
\name{methylDiff-class}
\alias{methylDiff-class}
\title{An S4 class that holds differential methylation information}
\description{
  This object is desgined to hold statistics and locations
  for differentially methylated regions/bases
}
\section{Slots}{
  \describe{ \item{\code{sample.ids}}{ids/names of samples
  in a vector} \item{\code{assembly}}{a name of genome
  assembly, such as :hg18,mm9, etc}
  \item{\code{context}}{numeric vector identifying which
  samples are which group } \item{\code{treatment}}{numeric
  vector identifying which samples are which group }
  \item{\code{destranded}}{logical denoting if methylation
  inormation is destranded or not}
  \item{\code{.Data}}{data.frame holding the locations and
  statistics}

  }
}


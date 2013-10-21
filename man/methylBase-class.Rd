\docType{class}
\name{methylBase-class}
\alias{methylBase}
\alias{methylBase-class}
\title{An S4 class for methylation events sampled in multiple experiments}
\description{
  This class is designed to contain methylation information
  such as coverage, number of methylated bases, etc.. The
  methylation events contained in the class must be sampled
  in multiple experiments (ex: only CpG bases covered in
  multiple experiments are stored in the object of this
  class). The class extends \code{data.frame} and creates
  an object that holds methylation information and genomic
  location. The object belonging to this class is produced
  by \code{\link{unite}} function.
}
\section{Slots}{
  \describe{ \item{\code{sample.ids}:}{character vector for
  ids of samples in the object}

  \item{\code{assembly}:}{name of the genome assembly}

  \item{\code{context}:}{context of methylation. Ex:
  CpG,CpH,CHH, etc}

  \item{\code{treatment}:}{treatment vector denoting which
  samples are test and control}

  \item{\code{coverage.index}:}{vector denoting which
  columns in the data correspons to coverage values}

  \item{\code{numCs.index}:}{vector denoting which columns
  in the data correspons to number of methylatedCs values}
  \item{\code{numTs.index}:}{vector denoting which columns
  in the data correspons to number of unmethylated Cs
  values} \item{\code{destranded}:}{ logical value. If
  \code{TRUE} object is destranded, if \code{FALSE} it is
  not.} \item{\code{resolution}:}{ resolution of
  methylation information, allowed values: 'base' or
  'region'} }
}

\section{Details}{
  \code{methylBase} class extends \code{\link{data.frame}}
  class therefore providing novice and experienced R users
  with a data structure that is well known and ubiquitous
  in many R packages.
}

\section{Subsetting}{
  In the following code snippets, \code{x} is a
  \code{methylDiff}.  Subsetting by \code{x[i,]} will
  produce a new object if subsetting is done on rows.
  Column subsetting is not directly allowed to prevent
  errors in the downstream analysis. see ?methylKit[ .
}

\section{Accessors}{
  The following functions provides access to data slots of
  methylDiff:
  \code{\link[methylKit]{getData}},\code{\link[methylKit]{getAssembly}},
  \code{\link[methylKit]{getContext}}
}

\section{Coercion}{
  \code{methylBase} object can be coerced to
  \code{\link[GenomicRanges]{GRanges}} object via
  \code{\link{as}} function.
}
\examples{
data(methylKit)
library(GenomicRanges)
my.gr=as(methylBase.obj,"GRanges")
}


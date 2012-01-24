\name{methylBase-class}
\alias{methylBase-class}
\title{An S4 class that holds base-pair resolution methylation information for multiple experiments, only bases that are covered in all experiments are held in this class}
\description{
  extends data.frame and creates an object that holds
  methylation information and genomic location
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
  values} \item{\code{resolution}:}{ resolution of
  methylation information, allowed values: 'base' or
  'region'} }
}


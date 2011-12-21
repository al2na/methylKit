\name{methylBase-class}
\alias{methylBase-class}
\title{An S4 class that holds base-pair resolution methylation information for multiple experiments, only bases that are covered in all experiments are held in this class}
\description{
  extends data.frame and creates an object that holds
  methylation information and genomic location
}
\section{Slots}{
  \describe{ \item{sample.ids}{character vector for ids of
  samples in the object}

  \item{assembly}{name of the genome assembly}

  \item{context}{context of methylation. Ex: CpG,CpH,CHH,
  etc}

  \item{treatment}{treatment vector denoting which samples
  are test and control}

  \item{coverage.index}{vector denoting which columns in
  the data correspons to coverage values}

  \item{numCs.index}{vector denoting which columns in the
  data correspons to number of methylatedCs values}

  \item{numTs.index}{vector denoting which columns in the
  data correspons to number of unmethylated Cs values} }
}


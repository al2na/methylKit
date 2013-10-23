\docType{methods}
\name{reconstruct}
\alias{reconstruct}
\title{Reconstruct methylBase object based on a new methylation percentage matrix}
\usage{
  reconstruct(methMat, mBase)
}
\arguments{
  \item{methMat}{percent methylation matrix, row order and
  order of the samples same as the methylBase object}

  \item{mBase}{\code{\link{methylBase}} object to be
  reconstructed}
}
\value{
  new \code{\link{methylBase}} object where methylation
  percentage matches input \code{methMat} and coverages
  matches input \code{mBase}
}
\description{
  The function reconstructs a new methylBase object from an
  input methylBase object and percent methylation matrix.
  Basically, it uses the read coverages in the input
  methylBase object and deduces new number of methylated Cs
  and unmethylated Cs based on the input percent
  methylation matrix. It is ideally to be used to
  reconstruct methylBase objects after batch correction on
  percent methylation values. The percent methylation
  matrix rows must match methylBase object rows in order
  ,and in addition column order (the order of samples) in
  input methylBase must match the order in percent
  methylation matrix.
}
\note{
  Batch effect correction (if any batch effect exists) is a
  tricky issue. We provide some simple ways to deal with it
  (see \code{\link{assocComp}} and \code{\link{removeComp}}
  ), But if you can find other ways to correct for batch
  effects and want to create a methylBase object with the
  corrected percent methylation values, you can use this
  function.
}
\examples{
data(methylKit)

# get percent methylation
mat=percMethylation(methylBase.obj)

# do some changes in the matrix
# this is just a toy example
# ideally you want to correct the matrix
# for batch effects
mat[mat==100]=80

# reconstruct the methylBase from the corrected matrix
newobj=reconstruct(mat,methylBase.obj)
}
\author{
  Altuna Akalin
}


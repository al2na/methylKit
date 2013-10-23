\docType{methods}
\name{percMethylation}
\alias{percMethylation}
\alias{percMethylation,methylBase-method}
\title{get percent methylation scores from methylBase object}
\usage{
  percMethylation(methylBase.obj, rowids = FALSE)
}
\arguments{
  \item{methylBase.obj}{a methylBase object}

  \item{rowids}{if TRUE, matrix rownames have identifiers
  as base/region location (default:FALSE)}
}
\value{
  matrix with percent methylation values per base/region
  across all samples, row names would be base/region
  identifiers
}
\description{
  get percent methylation scores from methylBase object
}
\examples{
data(methylKit)
mat=percMethylation(methylBase.obj)
head(mat)
}


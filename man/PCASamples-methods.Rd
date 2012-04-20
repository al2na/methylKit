\docType{methods}
\name{PCASamples}
\alias{PCASamples}
\alias{PCASamples,methylBase-method}
\title{CpG Dinucleotide Methylation Principal Components Analysis}
\usage{
  PCASamples(.Object, cor=TRUE, screeplot=FALSE,
  adj.lim=c(0.0004,0.1),scale=TRUE,center=TRUE)
}
\arguments{
  \item{.Object}{a \code{methylBase} object}

  \item{cor}{[Not used anymore] cor a logical value
  indicating whether the calculation should use the
  correlation matrix or the covariance matrix. (default:
  TRUE)}

  \item{screeplot}{a logical value indicating whether to
  plot the variances against the number of the principal
  component. (default: FALSE)}

  \item{adj.lim}{a vector indicating the propotional
  adjustment of xlim (adj.lim[1]) and ylim (adj.lim[2]).
  (default: c(0.0004,0.1))}

  \item{scale}{logical indicating if \code{prcomp} should
  scale the data to have unit variance or not (default:
  TRUE)}

  \item{center}{logical indicating if \code{prcomp} should
  center the data or not (default: TRUE)}
}
\value{
  The form of the value returned by \code{PCASamples} is
  the summary of principal component analysis by
  \code{prcomp}.
}
\description{
  CpG Dinucleotide Methylation Principal Components
  Analysis
}
\note{
  cor option is not in use anymore, since \code{prcomp} is
  used for PCA analysis instead of \code{princomp}
}


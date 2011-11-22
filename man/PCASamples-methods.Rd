\docType{methods}
\name{PCASamples}
\alias{"PCASamples,"}
\alias{"PCASamples,-methods"}
\alias{"PCASamples,ANY-method"}
\alias{"methylBase-method"}
\alias{PCASamples}
\title{CpG Dinucleotide Methylation Principal Components Analysis}
\arguments{
  \item{.Object}{a \code{methylBase} object}

  \item{cor}{a logical value indicating whether the
  calculation should use the correlation matrix or the
  covariance matrix. (default: TRUE)}

  \item{screeplot}{a logical value indicating whether to
  plot the variances against the number of the principal
  component. (default: FALSE)}
}
\value{
  The form of the value returned by \code{PCASamples} is
  the summary of principal component analysis by
  \code{princomp}.
}
\description{
  CpG Dinucleotide Methylation Principal Components
  Analysis
}


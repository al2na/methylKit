\docType{methods}
\name{getCoverageStats}
\alias{getCoverageStats}
\alias{getCoverageStats,-methods}
\alias{getCoverageStats,methylRaw-method}
\title{get coverage stats from methylRaw object}
\arguments{
  \item{.Object}{a \code{methylRaw} object}

  \item{plot}{plot a histogram of coverage if TRUE
  (default:FALSE)}

  \item{both.strands}{do stats and plot for both strands if
  TRUE (default:FALSE)}

  \item{labels}{should the bars of the histrogram have
  labels showing the percentage of values in each bin
  (default:TRUE)}

  \item{...}{options to be passed to \code{hist} function}
}
\value{
  a summary of coverage statistics or plot a histogram of
  coverage
}
\description{
  get coverage stats from methylRaw object
}


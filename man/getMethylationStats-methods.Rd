\docType{methods}
\name{getMethylationStats}
\alias{getMethylationStats}
\alias{getMethylationStats,methylRaw-method}
\title{get Methylation stats from methylRaw object}
\usage{
  getMethylationStats(.Object,plot=FALSE,both.strands=FALSE,labels=TRUE,...)
}
\arguments{
  \item{.Object}{a \code{methylRaw} object}

  \item{plot}{plot a histogram of Methylation if TRUE
  (deafult:FALSE)}

  \item{both.strands}{do plots and stats for both strands
  seperately if TRUE (deafult:FALSE)}

  \item{labels}{should the bars of the histrogram have
  labels showing the percentage of values in each bin
  (default:TRUE)}

  \item{...}{options to be passed to \code{hist} function.}
}
\value{
  a summary of Methylation statistics or plot a histogram
  of coverage
}
\description{
  get Methylation stats from methylRaw object
}


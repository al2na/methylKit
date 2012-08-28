\docType{methods}
\name{getCoverageStats}
\alias{getCoverageStats}
\alias{getCoverageStats,methylRaw}
\alias{getCoverageStats,methylRaw-method}
\title{get coverage stats from methylRaw object}
\usage{
  getCoverageStats(.Object,plot=FALSE,both.strands=FALSE,labels=TRUE,...)
}
\arguments{
  \item{.Object}{a \code{methylRaw} object}

  \item{plot}{plot a histogram of coverage if TRUE
  (default:FALSE)}

  \item{both.strands}{do stats and plot for both strands if
  TRUE (default:FALSE)}

  \item{labels}{should the bars of the histrogram have
  labels showing the percentage of values in each bin
  (default:TRUE)}

  \item{...}{options to be passed to
  \code{\link[graphics]{hist}} function}
}
\value{
  a summary of coverage statistics or plot a histogram of
  coverage
}
\description{
  The function returns basic statistics about read coverage
  per base. It can also plot a histogram of read coverage
  values.
}
\examples{
data(methylKit)

# gets coverage stats for the first sample in methylRawList.obj object
getCoverageStats(methylRawList.obj[[1]],plot=TRUE,both.strands=FALSE,labels=TRUE)
}


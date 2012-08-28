\docType{methods}
\name{diffMethPerChr}
\alias{diffMethPerChr}
\alias{diffMethPerChr,methylDiff-method}
\title{Get and plot the number of hyper/hypo methylated regions/bases per chromosome}
\usage{
  diffMethPerChr(x,plot=T,qvalue.cutoff=0.01,
  meth.cutoff=25,exclude=NULL,...)
}
\arguments{
  \item{x}{a \code{\link{methylDiff}} object}

  \item{plot}{TRUE|FALSE. If TRUE horizontal barplots for
  proportion of hypo/hyper methylated bases/regions}

  \item{qvalue.cutoff}{cutoff for q-value}

  \item{meth.cutoff}{cutoff for percent methylation
  difference}

  \item{exclude}{names of chromosomes to be excluded}

  \item{...}{extra graphical parameters to be passed to
  \code{barplot} function}
}
\value{
  plots a piechart or a barplot for percentage of the
  target features overlapping with annotation
}
\description{
  This function gets number of hyper/hypo methylated
  regions/bases from \code{\link{methylDiff}} object. It
  can also plot percentages of differentially methylated
  bases per chromosome.
}
\examples{
data(methylKit)

# get a list of differentially methylated bases/regions per chromosome and overall
diffMethPerChr(methylDiff.obj, plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25,exclude=NULL)
}


\docType{methods}
\name{clusterSamples}
\alias{clusterSamples}
\alias{clusterSamples,methylBase-method}
\title{Hierarchical Clustering using methylation data}
\usage{
  clusterSamples(.Object, dist="correlation",
  method="ward", plot=TRUE)
}
\arguments{
  \item{.Object}{a \code{methylBase} object}

  \item{dist}{the distance measure to be used. This must be
  one of "\code{correlation}", "\code{euclidean}",
  "\code{maximum}", "\code{manhattan}", "\code{canberra}",
  "\code{binary}" or "\code{minkowski}". Any unambiguous
  substring can be given. (default:"\code{correlation}")}

  \item{method}{the agglomeration method to be used. This
  should be (an unambiguous abbreviation of) one of
  "\code{ward}", "\code{single}", "\code{complete}",
  "\code{average}", "\code{mcquitty}", "\code{median}" or
  "\code{centroid}". (default:"\code{ward}")}

  \item{plot}{a logical value indicating whether to plot
  hierarchical clustering. (default:TRUE)}
}
\value{
  a \code{tree} object of a hierarchical cluster analysis
  using a set of dissimilarities for the n objects being
  clustered.
}
\description{
  The function clusters samples using
  \code{\link[stats]{hclust}} function and various distance
  metrics derived from percent methylation per base or per
  region for each sample.
}
\examples{
data(methylKit)

clusterSamples(methylBase.obj, dist="correlation", method="ward", plot=TRUE)
}


\docType{methods}
\name{clusterSamples}
\alias{"clusterSamples,-methods"}
\alias{"getCorrelation,methylBase-method"}
\alias{clusterSamples}
\title{Hierarchical cluster analysis on samples in methylBase object}
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

  \item{plot}{clustering plot if TRUE (default:TRUE)}
}
\value{
  a \code{tree} object produced by hclust and plot
  hierarchical clustering
}
\description{
  Hierarchical cluster analysis on samples in methylBase
  object
}


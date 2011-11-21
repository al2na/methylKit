\docType{methods}
\name{clusterSamples}
\alias{"clusterSamples,-methods"}
\alias{"getCorrelation,methylBase-method"}
\alias{clusterSamples}
\title{Hierarchical cluster analysis on samples in methylBase object}
\arguments{
  \item{.Object}{a methylBase object}

  \item{dist}{the distance measure to be used. This must be
  one of "correlation", "euclidean", "maximum",
  "manhattan", "canberra", "binary" or "minkowski". Any
  unambiguous substring can be given.
  (default:"correlation")}

  \item{method}{the agglomeration method to be used. This
  should be (an unambiguous abbreviation of) one of "ward",
  "single", "complete", "average", "mcquitty", "median" or
  "centroid". (default:"ward")}

  \item{plot}{clustering plot if TRUE (default:TRUE)}
}
\value{
  a tree object produced by hclust and plot hierarchical
  clustering
}
\description{
  Hierarchical cluster analysis on samples in methylBase
  object
}


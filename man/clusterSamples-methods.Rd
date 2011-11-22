\docType{methods}
\name{clusterSamples}
\alias{"clusterSamples,"}
\alias{"clusterSamples,-methods"}
\alias{"clusterSamples,ANY-method"}
\alias{"methylBase-method"}
\alias{clusterSamples}
\title{CpG Dinucleotide Methylation Hierarchical Cluster Analysis}
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
  CpG Dinucleotide Methylation Hierarchical Cluster
  Analysis
}


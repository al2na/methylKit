\docType{methods}
\name{get.methylDiff}
\alias{get.methylDiff}
\alias{get.methylDiff,methylDiff-method}
\title{gets differentially methylated regions/bases based on cutoffs}
\usage{
  get.methylDiff(.Object,difference=25,qvalue=0.01,type="all")
}
\arguments{
  \item{.Object}{a methylDiff object}

  \item{difference}{cutoff for absolute value of %
  methylation change between test and control (default:25)}

  \item{qvalue}{cutoff for qvalue of differential
  methylation statistic (default:0.01)}

  \item{type}{one of the "hyper","hypo" or "all" strings.
  Specifies what type of differentially menthylated
  bases/regions should be returned. For retrieving
  Hyper-methylated regions/bases type="hyper", for
  hypo-methylated type="hypo" (default:"all")}
}
\value{
  a methylDiff object containing the differential
  methylated locations satisfying the criteria
}
\description{
  gets differentially methylated regions/bases based on
  cutoffs
}


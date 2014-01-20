\docType{methods}
\name{calculateDiffMethDSS}
\alias{calculateDiffMethDSS}
\alias{calculateDiffMethDSS,methylBase-method}
\title{Calculate differential methylation statistics using DSS package}
\usage{
  calculateDiffMeth(.Object,slim=TRUE, num.cores=1)
}
\arguments{
  \item{.Object}{a methylBase object to calculate
  differential methylation}

  \item{slim}{If TRUE(default) SLIM method will be used for
  P-value adjustment.  If FALSE, \code{\link{p.adjust}}
  with method="fdr" option (same as Benjamini-Hochberg "BH"
  option) will be used for P-value
  correction.}



  \item{num.cores}{integer for denoting how many cores
  should be used for differential methylation calculations
  (only can be used in machines with multiple cores)
  (NOT CURRENTLY IMPLEMENTED)
  }
}
\value{
  a methylDiff object containing the differential
  methylation statistics and locations
}
\description{
  The function calculates differential methylation
  statistics between two groups of samples using
  a hierarchical Bayesian model
  
}
\section{Details}{
  The function uses code from the DSS package, which uses
  a hierarchical Bayesian model to estimate methylation
  difference between two groups. Read counts are modeled
  with a beta-binomial distribution,
  with a log-normal prior on the methylation disperions 
  for each group.
  
  Currently doesn't work with NA values - use the unite()
  function on a methylRawList object with the min.per.group
  parameter set to NULL (the default value) to ensure there are
  no NAs.
}
\examples{
data(methylKit)

# Use hierarchicalBayesian model to calculate differential methylation
# between two groups

my.diffMethDSS=calculateDiffMethDSS(methylBase.obj,slim=TRUE,weigthed.mean=TRUE,num.cores=1)


}
\references{
  Altuna Akalin, Matthias Kormaksson, Sheng Li, Francine E.
  Garrett-Bakelman, Maria E. Figueroa, Ari Melnick,
  Christopher E. Mason. (2012).  "methylKit: A
  comprehensive R package for the analysis of genome-wide
  DNA methylation profiles." Genome Biology.
  
  Hao Wu, Chi Wang, Zhijin Wu (2012): A new shrinkage estimator for
  dispersion improves differential expression detection in RNA-seq
  data.  Biostatistics:10.1093/biostatistics/kxs033
  
  
}
\seealso{
  \code{\link[methylKit]{pool}},
  \code{\link[methylKit]{reorganize}}
}
\docType{methods}
\name{calculateDiffMeth}
\alias{calculateDiffMeth}
\alias{calculateDiffMeth,methylBase-method}
\title{Calculate differential methylation statistics}
\usage{
  calculateDiffMeth(.Object,slim=TRUE,weighted.mean=TRUE,num.cores=1)
}
\arguments{
  \item{.Object}{a methylBase object to calculate
  differential methylation}

  \item{slim}{If TRUE(default) SLIM method will be used for
  P-value adjustment.  If FALSE, \code{\link{p.adjust}}
  with method="BH" option will be used for P-value
  correction.}

  \item{weighted.mean}{calculate the mean methylation
  difference between groups using read coverage as weights}

  \item{num.cores}{integer for denoting how many cores
  should be used for differential methylation calculations
  (only can be used in machines with multiple cores)}
}
\value{
  a methylDiff object containing the differential
  methylation statistics and locations
}
\description{
  The function calculates differential methylation
  statistics between two groups of samples. The function
  uses either logistic regression test or Fisher's Exact
  test to calculate differential methylation. See
  references for detailed explanation on statistics.
}
\section{Details}{
  The function either uses a logistic regression (when
  there are multiple samples per group) or fisher's exact
  when there is one sample per group.
}
\examples{
data(methylKit)

# Logistic regression test will be applied since there are multiple samples in each group
# in methylBase.obj object
my.diffMeth=calculateDiffMeth(methylBase.obj,slim=TRUE,weighted.mean=TRUE,num.cores=1)

# P-value adjustment with Benjamini-Hocberg via p.adjust
my.diffMeth=calculateDiffMeth(methylBase.obj,slim=FALSE,weighted.mean=TRUE)


# pool samples in each group
pooled.methylBase=pool(methylBase.obj,sample.ids=c("test","control"))

# After applying pool() function, there is one sample in each group.
# Fisher's exact test will be applied for differential methylation
my.diffMeth2=calculateDiffMeth(pooled.methylBase,slim=TRUE,
                               weighted.mean=TRUE,num.cores=1)
}
\references{
  Altuna Akalin, Matthias Kormaksson, Sheng Li, Francine E.
  Garrett-Bakelman, Maria E. Figueroa, Ari Melnick,
  Christopher E. Mason. (2012).  "methylKit: A
  comprehensive R package for the analysis of genome-wide
  DNA methylation profiles." Genome Biology.
}
\seealso{
  \code{\link[methylKit]{pool}},
  \code{\link[methylKit]{reorganize}}
}


\docType{methods}
\name{reorganize}
\alias{reorganize}
\alias{reorganize,methylBase-method}
\alias{reorganize,methylRawList-method}
\title{reorganize methylRawList and methylBase objects by creating new objects from subset of samples}
\arguments{
  \item{methylObj}{a \code{methylRawList} or
  \code{methylBase} object}

  \item{sample.ids}{a vector for sample.ids to be subset.
  Order is important and the order should be similar to
  treatment. sample.ids should be a subset or reordered
  version of sample ids in the input object.}

  \item{treatment}{treatment vector, should be same length
  as sample.ids vector}
}
\value{
  RETURNS a \code{methylRawList} or \code{methylBase}
  object depending on the input object
}
\description{
  Create a new \code{methylRawList} or \code{methylBase}
  object by selecting a subset of samples from the input
  object, which is a \code{methylRawList} or
  \code{methylBase} object. You can use the function to
  partition a large methylRawList or methylBase object to
  smaller object based on sample ids or when you want to
  reorder samples and treatmet vector.
}
\examples{
# this is a list of example files, ships with the package
file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
               system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
               system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
               system.file("extdata", "control2.myCpG.txt", package = "methylKit") )


# read the files to a methylRawList object: myobj
myobj=read( file.list,
          sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))
meth=unite(myobj,destrand=TRUE)

myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
meth2 =reorganize(meth,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
}


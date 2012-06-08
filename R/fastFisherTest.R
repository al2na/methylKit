dyn.load(paste("R/fastFisherTest", .Platform$dynlib.ext, sep=""))

.First.lib <- function(lib, pkg)
{
    library.dynam("LIBNAME", pkg, lib)
}

#' fisher's exact test
#'
#' The functions calculate the two-sided pvalue of fisher's exact test. 
#' @param x a \code{matrix} with four columns, the numCs_sample1, numTs_sample1, numCs_sample2, numTs_sample2 
#'
#' @usage fast.fisher.test(x)
#'
#' @author Sheng Li and Altuna Akalin
# @examples
# library(methylKit)
# data(methylKit)
setGeneric("fast.fisher.test", function(x) standardGeneric("fast.fisher.test"))

#' @rdname fast.fisher.test
#' @aliases fast.fisher.test
setMethod("fast.fisher.test", "matrix",
    function(x){

    if(is.matrix(x)){
      a=as.integer(x[,1])
      b=as.integer(x[,2])
      c=as.integer(x[,3])
      d=as.integer(x[,4])
      n <- as.integer(length(a))
      r <- .C("fastFisherTest",a,b,c,d,n, pval=double(n))$pval
      return (r);
    } else {
      show("input have to be a matrix with 4 columns")
    }

})

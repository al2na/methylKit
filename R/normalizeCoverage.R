#' normalize read coverage between samples
#'
#' The function normalizes coverage values between samples using a scaling factor derived from differences between mean or median of coverage distributions
#' @param obj  \code{methylRawList} or \code{methylRawListDB} object 
#' @param method  a string "mean" or "median" which denotes median or mean should be used to calculate scaling factor. (Default:median)
#' @param chunk.size Number of rows to be taken as a chunk for processing the \code{methylRawListDB} objects. (Default: 1e6)
#' @return a  \code{methylRawList} or \code{methylRawList} object depending on class of input object
#'
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the \code{chunk.size}.
#' 
#' @usage normalizeCoverage(obj,method="median")
#' @author  Altuna Akalin
#' @export
#' @examples
#' 
#' data(methylKit)
#' newObj=normalizeCoverage(methylRawList.obj,method="median")
#'
#' @docType methods
#' @rdname normalizeCoverage-methods
setGeneric("normalizeCoverage", function(obj,method="median",chunk.size=1e6) standardGeneric("normalizeCoverage"))

#' @rdname normalizeCoverage-methods
#' @aliases normalizeCoverage,methylRawList-method
setMethod("normalizeCoverage", "methylRawList",
                    function(obj,method){


                      if(method=="median"){
                        x=sapply(obj,function(x) median(x$coverage) )
                      }else if(method=="mean"){
                        x=sapply(obj,function(x) mean(x$coverage) )
                      }else{
                        stop("method option should be either 'mean' or 'median'\n")
                      }
                      sc.fac=max(x)/x #get scaling factor

                      for(i in 1:length(obj)){
                        all.cov=obj[[i]]$coverage
                        fCs    =obj[[i]]$numCs/all.cov
                        fTs    =obj[[i]]$numT/all.cov
                        obj[[i]]$coverage=round(sc.fac[i]*obj[[i]]$coverage)
                        obj[[i]]$numCs   =round(obj[[i]]$coverage*fCs)
                        obj[[i]]$numTs   =round(obj[[i]]$coverage*fTs)
                      }
                      obj
})

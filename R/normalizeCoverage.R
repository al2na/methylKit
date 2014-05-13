#' normalize read coverage between samples
#'
#' The function normalizes coverage values between samples using a scaling factor derived from differences between mean or median of coverage distributions
#' @param obj  \code{methylRawList} object 
#' @param method  a string "mean" or "median" which denotes median or mean should be used to calculate scaling factor. (Default:median)
#' @return a  \code{methylRawList} object
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
setGeneric("normalizeCoverage", function(obj,method="median") standardGeneric("normalizeCoverage"))

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
                        fTs    =obj[[i]]$numTs/all.cov
                        obj[[i]]$coverage=round(sc.fac[i]*obj[[i]]$coverage)
                        obj[[i]]$numCs   =round(obj[[i]]$coverage*fCs)
                        obj[[i]]$numTs   =round(obj[[i]]$coverage*fTs)
                      }
                      obj
})

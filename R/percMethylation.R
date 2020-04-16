#' get percent methylation scores from methylBase or methylBaseDB object
#' 
#' @param methylBase.obj a methylBase or methylBaseDB object 
#' @param rowids if TRUE, matrix rownames have identifiers as base/region 
#' location (default:FALSE)
#' @param save.txt if TRUE, the matrix will be written to a text file, 
#' but only for methylBaseDB objects (default: FALSE) 
#' @param chunk.size Number of rows to be taken as a chunk for processing 
#' the \code{methylBaseDB} objects (default: 1e6)
#' @return matrix with percent methylation values per base/region across all 
#'         samples, row names would be base/region identifiers
#' @examples
#' 
#' data(methylKit)
#' mat=percMethylation(methylBase.obj)
#' head(mat)
#' 
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with 
#' \code{methylBaseDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized 
#' objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for
#'  most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the
#'  \code{chunk.size}.
#' 
#' @export
#' @docType methods
#' @rdname percMethylation-methods
setGeneric("percMethylation", function(methylBase.obj,rowids=FALSE,
                                       save.txt=FALSE,chunk.size=1e6) 
  standardGeneric("percMethylation"))

#' @rdname percMethylation-methods
#' @aliases percMethylation,methylBase-method
setMethod("percMethylation", "methylBase",
                    function(methylBase.obj,rowids=FALSE){
    x=getData(methylBase.obj)
    meth.mat = 100 * x[, methylBase.obj@numCs.index]/(
      x[,methylBase.obj@numCs.index] + x[,methylBase.obj@numTs.index] )                                      
    names(meth.mat)=methylBase.obj@sample.ids
    rownames(meth.mat) <- NULL
    if(rowids){
      rownames(meth.mat)=as.character(paste(x[,1],x[,2],x[,3],sep=".") )
    }
    return(as.matrix(meth.mat))
  
})
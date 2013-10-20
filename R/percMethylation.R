#' get percent methylation scores from methylBase object
#' 
#' @param methylBase.obj a methylBase object 
#' @param rowids if TRUE, matrix rownames have identifiers as base/region location (default:FALSE)
#' @return matrix with percent methylation values per base/region across all 
#'         samples, row names would be base/region identifiers
#' @usage percMethylation(methylBase.obj)
#' @examples
#' 
#' data(methylKit)
#' mat=percMethylation(methylBase.obj)
#' head(mat)
#' 
#' 
#' @export
#' @docType methods
#' @rdname percMethylation-methods
setGeneric("percMethylation", function(methylBase.obj,rowids=FALSE) standardGeneric("percMethylation"))

#' @rdname percMethylation-methods
#' @aliases percMethylation,methylBase-method
setMethod("percMethylation", "methylBase",
                    function(methylBase.obj,rowids=FALSE){
                        x=getData(methylBase.obj)
                        meth.mat = 100 * x[, methylBase.obj@numCs.index]/(x[,methylBase.obj@numCs.index] + x[,methylBase.obj@numTs.index] )                                      
                        names(meth.mat)=methylBase.obj@sample.ids
                        if(rowids){
                          rownames(meth.mat)=as.character(paste(x[,1],x[,2],x[,3],sep=".") )
                        }
                        return(as.matrix(meth.mat))
                      
})
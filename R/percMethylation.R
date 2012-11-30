#' get percent methylation scores from methylBase object
#' 
#' @param methylBase.obj a methylBase object 
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
setGeneric("percMethylation", function(methylBase.obj) standardGeneric("percMethylation"))

#' @rdname percMethylation-methods
#' @aliases percMethylation,methylBase-method
setMethod("percMethylation", "methylBase",
                    function(methylBase.obj){
                        x=getData(methylBase.obj)
                        meth.mat = 100 * x[, methylBase.obj@numCs.index]/(x[,methylBase.obj@numCs.index] + x[,methylBase.obj@numTs.index] )                                      
                        names(meth.mat)=methylBase.obj@sample.ids
                        rownames(meth.mat)=as.character(x[,1])
                        
                        return(as.matrix(meth.mat))
                      
})
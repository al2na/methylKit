
#' update methylKit objects 
#'  
#' The method updates object from earlier versions (<v1.0.0) to latest object.
#' 
#' @param object a methylKit object: methylRaw, methylRawList, methylBase or methylDiff
#'
#' @usage updateMethObject(objects)
#' @return \code{\link{methylRaw}},\code{\link{methylDiff}},\code{\link{methylBase}}
#' or \code{\link{methylRawList}} object
#'
#'
#'  @export 
#'  @docType methods
#'  @rdname updateMethObject
setGeneric("updateMethObject",function(object) standardGeneric("updateMethObject"))

#'  @aliases updateMethObject,methylRaw-method
#'  @rdname updateMethObject
setMethod("updateMethObject" ,signature(object = "methylRaw" ),
          function(object){
            
            new("methylRaw",getData(object)[,-1],sample.id=object@sample.id,
                assembly=object@assembly,context=object@context,resolution=object@resolution)
          })

#' @aliases updateMethObject,methylRawList-method
#'  @rdname updateMethObject
setMethod("updateMethObject" ,signature(object = "methylRawList" ),
          function(object){
            
            lst=lapply(object,updateMethObject)
            new("methylRawList",lst,treatment=object@treatment)
          })


#'  @aliases updateMethObject,methylBase-method
#'  @rdname updateMethObject
setMethod("updateMethObject" ,signature(object = "methylBase" ),
          function(object){
            
            new("methylBase",getData(object)[,-1],
                sample.ids=object@sample.ids,
                assembly=object@assembly,
                context=object@context,
                treatment=object@treatment,
                coverage.index=object@coverage.index-1,
                numCs.index=object@numCs.index-1,
                numTs.index=object@numTs.index-1,
                destranded=object@destranded,
                resolution=object@resolution
            )
            
          })

#'  @aliases updateMethObject,methylDiff-method
#'  @rdname updateMethObject
setMethod("updateMethObject" ,signature(object = "methylDiff" ),
          function(object){
            
            new("methylDiff",getData(object)[,-1],
                sample.ids=object@sample.ids,
                assembly=object@assembly,
                context=object@context,
                treatment=object@treatment,
                destranded=object@destranded,
                resolution=object@resolution
            )
          })
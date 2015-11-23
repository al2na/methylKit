#' normalize read coverage between samples
#'
#' The function normalizes coverage values between samples using a scaling factor derived from differences between mean or median of coverage distributions
#' @param obj  \code{methylRawList} or \code{methylRawListDB} object 
#' @param method  a string "mean" or "median" which denotes median or mean should be used to calculate scaling factor. (Default:median)
#' @param chunk.size Number of rows to be taken as a chunk for processing the \code{methylRawListDB} objects. (Default: 1e6)
#' @param save.db A Logical to decide whether the resulting object should be saved as flat file database or not, default: explained in Details sections  
#' @param ... optional Arguments used when save.db is TRUE
#'            
#'            \code{suffix}
#'                  A character string to append to the name of the output flat file database, 
#'                  only used if save.db is true, default actions: append \dQuote{_filtered} to current filename 
#'                  if database already exists or generate new file with filename \dQuote{sampleID_filtered}
#'                  
#'            \code{dbdir} 
#'                  The directory where flat file database(s) should be stored, defaults
#'                  to getwd(), working directory for newly stored databases
#'                  and to same directory for already existing database
#'                  
#            \code{dbtype}
#                  The type of the flat file database, currently only option is "tabix"
#                  (only used for newly stored databases)
#'
#' @return a  \code{methylRawList} or \code{methylRawList} object depending on class of input object
#'
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects as \code{methylRawListDB}, 
#' while being per default FALSE for \code{methylRawList}. If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the database allows the calculation in-memory, 
#' then you might want to change the value of this parameter.
#' 
#' @usage normalizeCoverage(obj,method="median",chunk.size,save.db,...)
#' @author  Altuna Akalin
#' @export
#' @examples
#' 
#' data(methylKit)
#' # normalize by the median coverage
#' newObj = normalizeCoverage(methylRawList.obj,method="median")
#' 
#' # normalize by mean coverage and save to database in folder methylDB
#' newDBObj = normalizeCoverage(methylRawList.obj,method="mean",save.db=TRUE,dbdir="methylDB")
#'
#' @docType methods
#' @rdname normalizeCoverage-methods
setGeneric("normalizeCoverage", function(obj,method="median",chunk.size=1e6,save.db=FALSE,...) standardGeneric("normalizeCoverage"))

#' @rdname normalizeCoverage-methods
#' @aliases normalizeCoverage,methylRawList-method
setMethod("normalizeCoverage", "methylRawList",
                    function(obj,method,save.db=FALSE,...){


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
                      
                      if(!save.db) {
                        obj
                      
                      } else {
                        
                        # catch additional args 
                        args <- list(...)
                        
                        if( !( "dbdir" %in% names(args)) ){
                          dbdir <- .check.dbdir(getwd())
                        } else { dbdir <- .check.dbdir(args$dbdir) }
                        if(!( "suffix" %in% names(args) ) ){
                          suffix <- NULL
                        } else { 
                          suffix <- paste0("_",args$suffix)
                        }
                        
                        # create new methylRawListDB
                        new.list=list()    
                        for(i in 1:length(obj)){
                          # create methylRawDB
                          tmp <- makeMethylRawDB(df=getData(obj[[i]]),dbpath=dbdir,dbtype="tabix",sample.id=paste0(obj[[i]]@sample.id,suffix),
                                                 assembly=obj[[i]]@assembly,context=obj[[i]]@context,resolution=obj[[i]]@resolution)
                          tmp@sample.id <- obj[[i]]@sample.id
                          new.list[[i]]=tmp
                        }
                        
                        new("methylRawListDB",new.list,treatment=obj@treatment)
                        
                      }
                      
                      
})






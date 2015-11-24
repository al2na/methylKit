#' Pool replicates within groups to a single sample per group
#'
#' The function sums up coverage, numCs and numTs values within each group so 
#' one representative sample for each group will be created in a new methylBase object
#' @param obj  \code{methylBase} or \code{methylBaseDB} object with two groups or more and each group 
#' should have multiple samples
#' @param sample.ids  a character vector of new sample.ids ex:c("test","control"), 
#'                     should follow the same order as unique treatment vector,
#'                    and should be equal to the length of the unique treatment
#'                     vector
#' @param chunk.size Number of rows to be taken as a chunk for processing the \code{methylRawListDB} objects, default: 1e6
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
#' @return a  \code{methylBase} or \code{methylBaseDB} object depending on class of input object
#'
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with \code{methylBaseDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects as \code{methylBaseDB}, 
#' while being per default FALSE for \code{methylBase}. If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the database allows the calculation in-memory, 
#' then you might want to change the value of this parameter.
#'
#' @usage pool(obj,sample.ids,chunk.size,save.db,...)
#' @author  Altuna Akalin
#' @export
#' @examples
#' library(methylKit)
#' data(methylKit)
#' 
#' # methylBase.obj has two groups, each group has two samples,
#' # the following function will pool the samples in each group
#' # so that each group will be represented by one pooled sample
#' pooled.methylBase=pool(methylBase.obj,sample.ids=c("test","control"))
#'
#' @docType methods
#' @rdname pool-methods
setGeneric("pool", function(obj,sample.ids,chunk.size=1e6,save.db=FALSE,...) standardGeneric("pool"))

#' @rdname pool-methods
#' @aliases pool,methylBase-method
setMethod("pool", "methylBase",
                    function(obj,sample.ids,save.db,...){
  df=getData(obj)

  treat=unique(obj@treatment)
  res=df[,1:4]
  for(i in 1:length(treat) ){

     # get indices
     setCs=obj@numCs.index[obj@treatment==treat[i]]
     setTs=setCs+1
     set.cov=setCs-1

     if(length(setCs)>1){
       Cs=rowSums(df[,setCs],na.rm=TRUE)
       Ts=rowSums(df[,setTs],na.rm=TRUE)
       covs=rowSums(df[,set.cov],na.rm=TRUE)

     }else{
       Cs  =df[,setCs]
       Ts  =df[,setTs]
       covs=df[,set.cov]
     }
     res=cbind(res,covs,Cs,Ts) # bind new data
     # change names of columns
     names(res)[(2:4)+3*i]=paste( c("coverage" ,"numCs" ,"numTs"),i,sep="") 

  }
  coverage.ind=3*(1:length(treat)) + 2
  
  if(!save.db) {
      obj1=new("methylBase",as.data.frame(res),sample.ids=sample.ids,
         assembly=obj@assembly,context=obj@context,
         treatment=treat,coverage.index=coverage.ind,
         numCs.index=coverage.ind+1,numTs.index=coverage.ind+2,
           destranded=obj@destranded,resolution=obj@resolution )
      obj1
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
    
    # create methylBaseDB
    makeMethylBaseDB(df=as.data.frame(res),dbpath=dbdir,dbtype="tabix",sample.ids=sample.ids,
                     assembly=obj@assembly,context=obj@context,
                     treatment=treat,coverage.index=coverage.ind,
                     numCs.index=coverage.ind+1,numTs.index=coverage.ind+2,
                     destranded=obj@destranded, resolution=obj@resolution,
                     suffix=suffix )
  }
})

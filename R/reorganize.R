# The function to reorganize the methylRawList and methylBase objects
# 
# 

#' Reorganize methylKit objects by creating new objects from subset of samples
#' 
#' The function creates a new  \code{methylRawList}, \code{methylRawListDB}, 
#' \code{methylBase} or \code{methylBaseDB} 
#' object by selecting a subset of samples from the input object, which is 
#' a \code{methylRawList} or \code{methylBase} object. You can use the function
#'  to partition a large methylRawList or methylBase object
#' to smaller object based on sample ids or when you want to reorder samples
#'  and/or give a new treatmet vector.
#'
#' @param methylObj a \code{methylRawList}, \code{methylRawListDB}, 
#' \code{methylBase} or \code{methylBaseDB} object
#' @param sample.ids a vector for sample.ids to be subset. Order is 
#'        important and the order should be similar to treatment. sample.ids 
#'        should be
#'        a subset or reordered version of sample ids in the input object.
#' @param treatment  treatment vector, should be same length as sample.ids vector
#' @param chunk.size Number of rows to be taken as a chunk for processing the 
#' \code{methylBaseDB} or \code{methylRawListDB} objects, default: 1e6
#' @param save.db A Logical to decide whether the resulting object should be 
#' saved as flat file database or not, default: explained in Details sections  
#' @param ... optional Arguments used when save.db is TRUE
#'            
#'            \code{suffix}
#'                  A character string to append to the name of the output 
#'                  flat file database, 
#'                  only used if save.db is true, default actions: 
#'                  append \dQuote{_filtered} to current filename 
#'                  if database already exists or generate new file 
#'                  with filename \dQuote{sampleID_filtered}
#'                  
#'            \code{dbdir} 
#'                  The directory where flat file database(s) should be 
#'                  stored, defaults
#'                  to getwd(), working directory for newly stored databases
#'                  and to same directory for already existing database
#'                  
#'           \code{dbtype}
#'                  The type of the flat file database, currently only option 
#'                 "tabix"
#                  (only used for newly stored databases)
#'
#' @return returns a \code{methylRawList}, \code{methylRawListDB}, 
#' \code{methylBase} or \code{methylBaseDB} object depending on the input object
#' @examples
#' # this is a list of example files, ships with the package
#' file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
#'
#'
#' # read the files to a methylRawList object: myobj
#' myobj=methRead( file.list,
#'           sample.id=list("test1","test2","ctrl1","ctrl2"),
#'           assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))
#' meth=unite(myobj,destrand=TRUE)
#'
#' # get samples named "test1" and "ctrl2" from myobj and create a new methylRawList object
#' myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
#' 
#' # # get samples named "test1" and "ctrl2" from meth and create a new methylBase object
#' meth2 =reorganize(meth,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
#'
#'
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with 
#' \code{methylBaseDB} or \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized 
#' objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the 
#' \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects as 
#' \code{methylBaseDB} and \code{methylRawListDB}, 
#' while being per default FALSE for \code{methylBase} and 
#' \code{methylRawList}. If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the 
#' database allows the calculation in-memory, 
#' then you might want to change the value of this parameter.
#' 
#' 
#' @export
#' @docType methods
#' @rdname reorganize-methods
setGeneric("reorganize", function(methylObj,sample.ids,treatment,chunk.size=1e6,
                                  save.db=FALSE,...) 
  standardGeneric("reorganize") )


#' @rdname reorganize-methods
#' @aliases reorganize,methylBase-method
setMethod("reorganize", signature(methylObj="methylBase"),
                    function(methylObj,sample.ids,treatment,save.db=FALSE,...){
                      
  #sample.ids length and treatment length should be equal
  if(length(sample.ids) != length(treatment) ){
    stop("length of sample.ids should be equal to treatment")
  }
  
  if( ! all(sample.ids %in% methylObj@sample.ids) ){
    stop("provided sample.ids is not a subset of the sample ids of the object")
  }

  temp.id = methylObj@sample.ids # get the subset of ids
  # get the column order in the original matrix
  col.ord = order(match(temp.id,sample.ids))[1:length(sample.ids)] 
  
  # make a matrix indices for easy access 
  ind.mat=rbind(methylObj@coverage.index[col.ord],  
                methylObj@numCs.index[col.ord],
                methylObj@numTs.index[col.ord])
  dat    =getData(methylObj) # get data
  
  #newdat =cbind(dat[,1:5],dat[ind.mat[,1]],dat[ind.mat[,2]],dat[ind.mat[,3]]) 
  # reorder columns using ind.mat
  newdat =dat[,1:4]
  for(i in 1:ncol(ind.mat))
  {
    newdat=cbind(newdat,dat[,ind.mat[,i]])
  }
                          
  # get indices of coverage,numCs and numTs in the data frame 
  coverage.ind=seq(5,by=3,length.out=length(sample.ids))
  numCs.ind   =coverage.ind+1
  numTs.ind   =coverage.ind+2
  
  # update column names
  #colnames(newdat)[6:ncol(newdat)]=paste(c("coverage","numCs","numTs"),
  # rep(1:length(sample.ids),each=length(sample.ids)),sep="")                        

  names(newdat)[coverage.ind]=paste(c("coverage"),1:length(sample.ids),sep="" )
  names(newdat)[numCs.ind]   =paste(c("numCs"),1:length(sample.ids),sep="" )
  names(newdat)[numTs.ind]   =paste(c("numTs"),1:length(sample.ids),sep="" )
  
  if(!save.db) {
    
    new("methylBase",as.data.frame(newdat),sample.ids=sample.ids,
         assembly=methylObj@assembly,context=methylObj@context,
         treatment=treatment,coverage.index=coverage.ind,
         numCs.index=numCs.ind,numTs.index=numTs.ind,
         destranded=methylObj@destranded, resolution=methylObj@resolution )
  
  } else {
    
    # catch additional args 
    args <- list(...)
    
    if( !( "dbdir" %in% names(args)) ){
      dbdir <- .check.dbdir(getwd())
    } else { dbdir <- .check.dbdir(args$dbdir) }
#                         if(!( "dbtype" %in% names(args) ) ){
#                           dbtype <- "tabix"
#                         } else { dbtype <- args$dbtype }
    if(!( "suffix" %in% names(args) ) ){
      suffix <- NULL
    } else { 
      suffix <- paste0("_",args$suffix)
    }
    
    # create methylBaseDB
    makeMethylBaseDB(df=as.data.frame(newdat),dbpath=dbdir,
                     dbtype="tabix",sample.ids=sample.ids,
                     assembly=methylObj@assembly,context=methylObj@context,
                     treatment=treatment,coverage.index=coverage.ind,
                     numCs.index=numCs.ind,numTs.index=numTs.ind,
                     destranded=methylObj@destranded, 
                     resolution=methylObj@resolution,suffix=suffix )
  }
})


#' @rdname reorganize-methods
#' @aliases reorganize,methylRawList-method
setMethod("reorganize", signature(methylObj="methylRawList"),
                    function(methylObj,sample.ids,treatment,save.db=FALSE,...){
                      
  #sample.ids length and treatment length should be equal
  if(length(sample.ids) != length(treatment) ){
    stop("length of sample.ids should be equal to treatment")
  }
  
  # get ids from the list of methylRaw 
  orig.ids=sapply(methylObj,function(x) x@sample.id) 
  if( ! all(sample.ids %in% orig.ids) ){
    stop("provided sample.ids is not a subset of the sample ids of the object")
  }
  
  # get the column order in the original matrix
  col.ord=order(match(orig.ids,sample.ids))[1:length(sample.ids)] 
  if(!save.db) {
    
    outList=list()    
    for(i in 1:length(sample.ids)){
      outList[[i]]=methylObj[[ col.ord[i]  ]]
    }

    new("methylRawList",outList,treatment=treatment)
  
  } else {
    
    # catch additional args 
    args <- list(...)
    
    if( !( "dbdir" %in% names(args)) ){
      dbdir <- .check.dbdir(getwd())
    } else { dbdir <- .check.dbdir(args$dbdir) }
#                         if(!( "dbtype" %in% names(args) ) ){
#                           dbtype <- "tabix"
#                         } else { dbtype <- args$dbtype }
    if(!( "suffix" %in% names(args) ) ){
      suffix <- NULL
    } else { 
      suffix <- paste0("_",args$suffix)
    }


    outList=list()    
    for(i in 1:length(sample.ids)){
      # create methylRawDB
      obj <- makeMethylRawDB(df=getData(methylObj[[ col.ord[i]  ]]),
                             dbpath=dbdir,dbtype="tabix",
                  sample.id=paste0(methylObj[[ col.ord[i]  ]]@sample.id,suffix),
                             assembly=methylObj[[ col.ord[i]  ]]@assembly,
                             context=methylObj[[ col.ord[i]  ]]@context,
                             resolution=methylObj[[ col.ord[i]  ]]@resolution)
      obj@sample.id <- methylObj[[ col.ord[i]  ]]@sample.id
      outList[[i]]=obj
    }
    
    new("methylRawListDB",outList,treatment=treatment)
    
    
  }
        
})
          
          
          

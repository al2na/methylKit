
valid.methylRawDB <- function(object) {
  
  
  #data=getData(object,nrow=5)
  check1=( (object@resolution == "base") | (object@resolution == "region") )
  check2=file.exists(object@dbpath)
  if(check1 & check2){
    return(TRUE)
  }
  else if (! check1 ){
    message("resolution slot has to be either 'base' or 'region':",
            "other values not allowed")
    FALSE
  }
  else if(! check2){
    message("The DB file can not be found, check the value of 'dbpath'")
    FALSE
  }
  
}


# methylRawDB
# dbpath -> filepath to flat file database, in this case tabix file
# num.records -> number of lines in tabix file
# dbtype currently only "tabix" allowed
setClass("methylRawDB", slots=list(dbpath="character",num.records="numeric",
  sample.id = "character", assembly = "character",context="character",
  resolution="character",dbtype="character"),validity=valid.methylRawDB)





# PRIVATE function:
# makes a methylRawDB object from a df
# it is called from modRead function or whenever this functionality is needed
makeMethylRawDB<-function(df,dbpath,dbtype,
                          sample.id, assembly ,context,
                          resolution){
  filepath=paste0(dbpath,"/",sample.id,".txt")
  df <- df[with(df,order(chr,start,end)),]
  df2tabix(df,filepath)
  num.records=Rsamtools::countTabix(paste0(filepath,".bgz"))[[1]] ## 
  
  new("methylRawDB",dbpath=paste0(filepath,".bgz"),num.records=num.records,
  sample.id = sample.id, assembly = assembly,context=context,
  resolution=resolution,dbtype=dbtype)
}

setMethod("show", "methylRawDB", function(object) {
  
  cat("methylRawDB object with",object@num.records,"rows\n--------------\n")
  print(headTabix(object@dbpath,nrow = 6,return.type = "data.frame"))
  cat("--------------\n")
  cat("sample.id:",object@sample.id,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("resolution:", object@resolution,"\n")
  cat("dbtype:", object@dbtype,"\n")
  cat("dbpath:",object@dbpath,"\n\n")
  
})


selectByOverlap<-function(object, ranges,inmemory=TRUE){
  
}

# methylRawListDB

valid.methylRawListDB <- function(object) {
  
  
  # if all elements are methyl
  if(!all(sapply(object,class)=="methylRawDB")){
    FALSE
  }
  else if ( length(object) != length(object@treatment) ){
    FALSE
  }
  else{
    TRUE
  }
  
}

setClass("methylRawListDB", slots=list(treatment = "vector"),contains = "list",
         validity=valid.methylRawListDB)

setMethod("show", "methylRawListDB", function(object) {
  
  cat("methylRawListDB object with",length(object),"methylRawDB objects\n\n")
  
  lapply(object,show)
  cat("treament:", object@treatment,"\n")
  
})






# methylBaseDB

#' @name methylBaseDB-class
#' @aliases methylBaseDB
#' @docType class
#' @rdname methylBaseDB-class
#' @export
setClass("methylBaseDB",contains="data.frame",representation(
  sample.ids = "character", assembly = "character",context = "character",treatment="numeric",coverage.index="numeric",
  numCs.index="numeric",numTs.index="numeric",destranded="logical",resolution = "character"))


# accessors ---------------------------------------------------------------


#' @rdname getAssembly-methods
#' @aliases getAssembly,methylRawDB-method
setMethod("getAssembly", signature="methylRawDB", definition=function(x) {
  return(x@assembly)
})

#' @rdname getAssembly-methods
#' @aliases getAssembly,methylRawListDB-method
setMethod("getAssembly", signature="methylRawListDB", definition=function(x) {
  return(lapply(x,getAssembly))
})

#' @rdname getContext-methods
#' @aliases getContext,methylRawDB-method
setMethod("getContext", signature="methylRawDB", definition=function(x) {
  return(x@context)
})

#' @rdname getContext-methods
#' @aliases getContext,methylRawListDB-method
setMethod("getContext", signature="methylRawListDB", definition=function(x) {
  return(lapply(x,getContext))
})

#' @rdname getData-methods
#' @aliases getData,methylRawDB-method
setMethod("getData", signature="methylRawDB", definition=function(x) {
  return(headTabix(tbxFile = x@dbpath, nrow = x@num.records, return.type = "data.frame"))
})

#' @rdname getData-methods
#' @aliases getData,methylRawListDB-method
setMethod("getData", signature="methylRawListDB", definition=function(x) {
  return(lapply(x,getData))
})




# subset classes ----------------------------------------------------------

#' @aliases select,methylRawDB-method
#' @rdname select-methods
setMethod("select", "methylRawDB",
          function(x, i)
          {
            
            new("methylRaw",getData(x)[i,],sample.id=x@sample.id,
                assembly=x@assembly,
                context=x@context,
                resolution=x@resolution)
          }
          
          
)

#' @examples
#' 
#' # This will get chromomsomes, will return a factor
#' # That means the resulting object will ceases to be a methylKit object
#' chrs=methylRawDB.obj[][[2]]
#' 
#' @aliases [,methylRawDB-method
#' @rdname extract-methods
setMethod("[", signature(x="methylRawDB", i = "ANY", j="ANY"),  
          function(x,i,j){
            #cat(missing(i),"\n",missing(j),"\n",missing(drop))
            if(!missing(j)){
              stop(paste("subsetting on columns is not allowed for",class(x),
                         "object\nif you want to do extraction on the data part", 
                         "of the object use getData() first"),
                   call. = FALSE)
            }
            select(x,i)
          }
)




# get/set values ----------------------------------------------------------



# get/set treatment vector

# get/set sample names
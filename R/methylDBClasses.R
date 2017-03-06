
#---------------------------------------------------------------------------------------
# regular R functions to be used in S4 functions


# set column names for methylRawDB and methylBaseDB data aquired from 
# flat file database
# @param df data.frame containing methylRaw or methylBase data
# @param methylDBclass 
.setMethylDBNames <- function(df,
                              methylDBclass=c("methylDB","methylBaseDB",
                                              "methylDiffDB")){
  
  if(nrow(df) == 0) return(df)
  
  if(missing(methylDBclass)){
        
    if( length(df) == 7 & unique(sapply(df,class)[5:7])=="integer"){
      setnames(x = df,old = names(df), 
                           new = c("chr","start","end","strand",
                                   "coverage","numCs","numTs"))
      
    } else if( length(df) == 7 & unique(sapply(df,class)[5:7])=="numeric"){
      setnames(x = df,old = names(df), 
                           new = c("chr","start","end","strand",
                                   "pvalue","qvalue","meth.diff")) 

    } else if( length(df) > 7){
      setnames(x = df,old = names(df)[1:4], 
                           new = c("chr","start","end","strand"))
      # get indices of coverage,numCs and numTs in the data frame 
      numsamples = (length(df)-4)/3
      coverage.ind=seq(5,by=3,length.out=numsamples)
      numCs.ind   =coverage.ind+1
      numTs.ind   =coverage.ind+2
      
      # change column names
      setnames(df,names(df)[coverage.ind], 
                           paste(c("coverage"),1:numsamples,sep="" ))
      setnames(df,names(df)[numCs.ind], 
                           paste(c("numCs"),1:numsamples,sep="" ))
      setnames(df,names(df)[numTs.ind], 
                           paste(c("numTs"),1:numsamples,sep="" ))
      
    } 
    
    #return(df)
    
  } else {
    
    if( methylDBclass == "methylRawDB" ){
      setnames(x = df,old = names(df), 
                           new = c("chr","start","end","strand",
                                   "coverage","numCs","numTs"))
    
    } else if ( methylDBclass == "methylBaseDB"){
      setnames(x = df,old = names(df)[1:4], 
                           new = c("chr","start","end","strand"))
      # get indices of coverage,numCs and numTs in the data frame 
      numsamples = (length(df)-4)/3
      coverage.ind=seq(5,by=3,length.out=numsamples)
      numCs.ind   =coverage.ind+1
      numTs.ind   =coverage.ind+2
      
      # change column names
      setnames(df,names(df)[coverage.ind], 
                           paste(c("coverage"),1:numsamples,sep="" ))
      setnames(df,names(df)[numCs.ind], 
                           paste(c("numCs"),1:numsamples,sep="" ))
      setnames(df,names(df)[numTs.ind], 
                           paste(c("numTs"),1:numsamples,sep="" ))
      
    } else if( methylDBclass == "methylDiffDB" ){
      setnames(x = df,old = names(df), 
                           new = c("chr","start","end","strand",
                                   "pvalue","qvalue","meth.diff"))
    
    #return(df)
    }
  }
}


# end of regular functions to be used in S4 functions
#---------------------------------------------------------------------------------------


# methylRawDB -------------------------------------------------------------



valid.methylRawDB <- function(object) {
  
  
  #data=getData(object,nrow=5)
  check1=( (object@resolution == "base") | (object@resolution == "region") )
  check2=file.exists(object@dbpath)
  check3=(length(object@sample.id) == 1)
  if (check2) check4 = (nrow(headTabix(object@dbpath)) != 0)
  if (check2) check5 = ( ncol(headTabix(object@dbpath)) == 7 )
  
  if (! check1 ){
    message("resolution slot has to be either 'base' or 'region':",
            "other values not allowed")
    FALSE
  }
  else if(! check2){
    cat("The DB file can not be found, check the value of 'dbpath'")
    FALSE
  }
  else if (! check3 ){
    message("object has more than one sample id:",object@sample.id,"\n",
        "only one allowed\n")
    FALSE
  } 
  else if(! check4) {
    
    ## NOTE: this check is done to hinder further issues with empty tabix files
    message("The tabix file for sample ",object@sample.id,
            " does not contain any data.\n",
            "Consider deleting file and associated index:\n",object@dbpath)
    FALSE
  }
  else if (! check5 ){
    cat("object does not have 7 columns, but has", 
        ncol(headTabix(object@dbpath)), "columns.",
        "This cannot be a methylRawDB Tabix file.\n")
    FALSE
  } 
  else {
    TRUE
  }
  
}


#' An S4 class for storing raw methylation data as flat file database.
#'
#' This object stores the raw mehylation data that is read in through read 
#' function as flat file database.The raw methylation data is basically
#' percent methylation values and read coverage values per genomic base/region.
#'
#' @section Slots:\describe{
#'  \item{\code{dbpath}:}{path to flat file database }
#'  \item{\code{num.records}:}{number of records (lines) in the object}
#'  \item{\code{sample.id}:}{string for an identifier of the sample}
#'  \item{\code{assembly}:}{string for genome assembly, ex: hg18,hg19,mm9}
#'  \item{\code{context}:}{ methylation context string, ex: CpG,CpH,CHH, etc.}
#'  \item{\code{resolution}:}{ resolution of methylation information, 'base' or 
#'  'region'}
#'  \item{\code{dbtype}:}{string for type of the flat file database, ex: tabix}
#'                 }
#' @section Details:
#' \code{methylRawDB} is created via \code{\link{read}} function and has the 
#' same functionality as \code{\link{methylRaw}} class, 
#' but the data is saved in a flat database file and therefore allocates less 
#' space in memory.
#' 
#' @section Subsetting:
#'  In the following code snippets, \code{x} is a \code{methylRawDB}.
#'  Subsetting by \code{x[i,]} will produce a new \code{methylRaw} object if 
#'  subsetting is done on
#'  rows. Column subsetting is not directly allowed to prevent errors in the 
#'  downstream analysis. see ?methylKit[ .
#'  \code{x[]} will return the \code{methylRawDB} object as new \code{methylRaw} object
#' 
#' @section Accessors: 
#' The following functions provides access to data slots of methylDiffDB:
#' - \code{\link{getData}}: get the data slot from the methylKit objects,
#' - \code{\link{getAssembly}}: get assembly of the genome,
#' - \code{\link{getContext}}: get the context of methylation
#' 
#' @section Coercion:
#'   \code{methylRawDB} object can be coerced to:
#'   \code{\link[GenomicRanges:GRanges-class]{GRanges}} object via \code{\link{as}} function.
#'   \code{\link{methylRaw}} object via \code{\link{as}} function.
#' 
#' @examples
#' 
#' # example of a raw methylation data contained as a text file
#' read.table(system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#' header=TRUE,nrows=5)
#' 
#' 
#' methylRawDB.obj <- methRead(
#' system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#'                         sample.id = "ctrl1", assembly = "hg18",
#'                         dbtype = "tabix", dbdir = "methylDB")
#' 
#' # example of a methylRawDB object
#' methylRawDB.obj
#' str(methylRawDB.obj)
#' 
#' library(GenomicRanges)
#' 
#' #coercing methylRawDB object to GRanges object
#' my.gr=as(methylRawDB.obj,"GRanges")
#' 
#' #coercing methylRawDB object to methylRaw object
#' myRaw=as(methylRawDB.obj,"methylRaw")
#' 
#' # remove Database again
#' rm(methylRawDB.obj)
#' unlink("methylDB",recursive=TRUE)
#' 
#' @name methylRawDB-class
#' @aliases methylRawDB
#' @docType class
#' @rdname methylRawDB-class
#' @export
setClass("methylRawDB", slots=list(dbpath="character",num.records="numeric",
  sample.id = "character", assembly = "character",context="character",
  resolution="character",dbtype="character"),validity=valid.methylRawDB)


# PRIVATE function:
# makes a methylRawDB object from a df
# it is called from read function or whenever this functionality is needed
makeMethylRawDB<-function(df,dbpath,dbtype,
                          sample.id, assembly ,context,
                          resolution,filename=NULL){
  if(dbtype != "tabix"){
    stop("unknown 'dbtype' argument provided, ",
         "currently only 'tabix' is accepted.")
  }
  
  # first check for filename length
  if(is.null(filename)) {
    filename <- paste0(dbpath,"/",sample.id,".txt")
    if( nchar( basename(filename) ) >= 245 ) {
      stop(paste("Generic Filename too long,\n",
                 "please manually provide filename.")) 
           }
  }
  
  # then we write the creation date ..
  write(paste0("#Date:",format(Sys.time(),'%Y-%m-%d %H:%M:%S')),
        file = filename)
  # add the class of the object
  write(paste0("#Class:","methylRawDB"),
        file = filename, append = TRUE)
  # and the slots as comments 
  
  num.records = nrow(df)

  # prepare slots as comments
  tabixHead <- paste0("#NR:",num.records,"\n",
                      "#SI:",sample.id,"\n",
                      "#AS:",assembly,"\n",
                      "#CT:",context,"\n",
                      "#RS:",resolution,"\n",
                      "#DT:",dbtype)
  
  # write slots to file
  write(tabixHead,
        file = filename ,append = TRUE)
  
  # sort the data
  df <- df[with(df,order(chr,start,end)),]
  # then we write the data
  write.table(x = df,
              file = filename, quote = FALSE,
              append = TRUE,col.names = FALSE,
              row.names = FALSE,sep = "\t")
  
  # and make tabix out of file
  makeMethTabix(filename,rm.file = FALSE)
  
  #filepath=paste0(dbpath,"/",sample.id,".txt")
  #df <- df[with(df,order(chr,start,end)),]
  #df2tabix(df,filepath) 
  #num.records=Rsamtools::countTabix(paste0(filepath,".bgz"))[[1]] ## 

  #new("methylRawDB",dbpath=paste0(filepath,".bgz"),num.records=num.records,
  new("methylRawDB",dbpath=paste0(filename,".bgz"),num.records=num.records,
  sample.id = sample.id, assembly = assembly,context=context,
  resolution=resolution,dbtype=dbtype)

}

# PRIVATE function:
# creates a methylRawDB object from a flat file database
# it is called from read function or whenever this functionality is needed
readMethylRawDB<-function(dbpath,dbtype,
                          sample.id, assembly ,context,
                          resolution,skip=0){

  if(!file.exists(paste0(dbpath,".tbi"))) 
  {  
    Rsamtools::indexTabix(dbpath,seq=1, start=2, end=3,
                          skip=skip, comment="#", zeroBased=FALSE)
  }
  num.records=Rsamtools::countTabix(dbpath)[[1]] ## 
  
  obj <- new("methylRawDB",dbpath=normalizePath(dbpath),num.records=num.records,
      sample.id = sample.id, assembly = assembly,context=context,
      resolution=resolution,dbtype=dbtype)
  
  if(valid.methylRawDB(obj)) {return(obj)}   
}

# PRIVATE function:
# creates a methylRawDB object from a flat file database and reads header content
# it is called from read function or whenever this functionality is needed
loadMethylRawDB<-function(dbpath,skip=0){
  
  if(!file.exists(paste0(dbpath,".tbi"))) 
  {  
    Rsamtools::indexTabix(dbpath,seq=1, start=2, end=3,
                          skip=skip, comment="#", zeroBased=FALSE)
  }
  
  h <- readTabixHeader(dbpath)
  
  if(! any(h$class == c("methylRaw", "methylRawDB") ) ) {
    stop(
      paste("Tabix file does not originate from methylRaw or methylRawDB.\n",
               "Please provide the correct class of data.")
      )
    }
  
  obj <- new("methylRawDB", dbpath=normalizePath(dbpath),
             num.records=h$num.records, sample.id = h$sample.ids, 
             assembly = h$assembly,context=h$context, 
             resolution=h$resolution,dbtype=h$dbtype)
  
  if(valid.methylRawDB(obj)) {return(obj)}   
}

# methylRawListDB

valid.methylRawListDB <- function(object) {
  
  
  # if all elements are methyl
  if(!all(sapply(object,class)=="methylRawDB")){
    message("It seems one the methylRawDB objects is invalid.")
    FALSE
  }
  else if ( length(object) != length(object@treatment) ){
    cat("The number of samples is different from the number of treatments, ","
            check the length of 'treatment'")
    FALSE
  }
  else{
    TRUE
  }
  
}


#' An S4 class for holding a list of methylRawDB objects.
#'
#' This class stores the list of  \code{\link{methylRawDB}} objects.
#' Functions such as \code{lapply} can be used on this list. It extends
#'  \code{\link[base]{list}} class. This object is primarily produced
#' by \code{\link{methRead}} function.
#'
#' @section Slots:\describe{
#'                  \item{\code{treatment}}{numeric vector denoting control
#'                   and test samples}
#'                  \item{\code{.Data}}{a list of \code{\link{methylRawDB}} objects  } 
#'                }
#'                
#' @section Constructor:\describe{
#'                  \item{\code{methylRawListDB(...)}}{combine multiple methylRawDB
#'                  objects supplied in ... into a methylRawListDB object.}
#'                }
#'                
#' @examples
#' file.list=list( system.file("extdata", "test1.myCpG.txt", 
#' package = "methylKit"),
#'                 system.file("extdata", "test2.myCpG.txt", 
#'                 package = "methylKit"),
#'                 system.file("extdata", "control1.myCpG.txt", 
#'                 package = "methylKit"),
#'                 system.file("extdata", "control2.myCpG.txt", 
#'                 package = "methylKit") )
#' 
#' methylRawListDB.obj <- methRead(file.list,
#' sample.id = list("test1","test2","ctrl1","ctrl2"),
#'                             assembly = "hg18",treatment = c(1,1,0,0),
#'                             dbtype = "tabix",dbdir = "methylDB")
#' 
#' #applying functions designed for methylRawDB on methylRawListDB object
#' lapply(methylRawListDB.obj,"getAssembly")
#' 
#' 
#' # remove Database again
#' rm(methylRawListDB.obj)
#' unlink("methylDB",recursive=TRUE)
#'
#' @name methylRawListDB-class
#' @aliases methylRawListDB
#' @docType class
#' @rdname methylRawListDB-class
#' @export
setClass("methylRawListDB", slots=list(treatment = "vector"),contains = "list",
         validity=valid.methylRawListDB)


### constructor function
methylRawListDB <- function(...,treatment) {
  
  listData <- list(...)
  ## check if input is really of type methylRaw 
  if (length(listData) == 0L) {
    stop("no methylRawDB object given.")
  } else {
    if (length(listData) == 1L && is.list(listData[[1L]]))
      listData <- listData[[1L]]
    if (!all(sapply(listData, is, "methylRawDB")))
      stop("all elements in '...' must be methylRawDB objects")
    
    ## just merge if valid
    mrl <- new("methylRawListDB", listData, treatment = treatment)
    if(valid.methylRawListDB(mrl)) return(mrl)
  }
}

# methylBaseDB ------------------------------------------------------------


valid.methylBaseDB <- function(object) {
  
  
  check1=( (object@resolution == "base") | (object@resolution == "region") )
  check2=file.exists(object@dbpath)
  check3=( length(object@sample.ids) == length(object@treatment) )
  numsamples = (length(df)-4)/3
  if (check2) {
    df <- headTabix(object@dbpath)
    numsamples = (length(df)-4)/3
    check4 = ( numsamples == length(object@sample.ids) )
  }
  
  if (! check1 ){
    cat("resolution slot has to be either 'base' or 'region':",
            "other values not allowed")
    FALSE
  }
  else if(! check2){
    cat("The DB file can not be found, check the value of 'dbpath'")
    FALSE
  }
  else if(! check3){
    cat("The number of samples is different from the number of treatments,",
            " check the length of 'sample.ids' & 'treatment'")
    FALSE
  }
  else if(! check4){
    cat("The number of samples is different from the number of sample names,",
            " check the length of 'sample.ids' & 'treatment'")
    FALSE
  }
  else {
    TRUE
  }
  
}


#' An S4 class for storing methylation events sampled in multiple experiments 
#' as flat file database
#'
#' This class is designed to contain methylation information such as coverage, 
#' number of methylated bases, etc...
#' The class creates an object that holds methylation information and genomic 
#' location as flat file database.
#' The object belonging to this class is produced by \code{\link{unite}} 
#' function.
#'          
#' @section Slots:\describe{
#'                  \item{\code{dbpath}:}{path to flat file database(s) }
#'                  \item{\code{num.records}:}{number of records (lines) 
#'                  in the object}
#'                  \item{\code{sample.ids}:}{character vector for ids of 
#'                  samples in the object}
#'                  \item{\code{assembly}:}{name of the genome assembly}
#'                  \item{\code{context}:}{context of methylation. 
#'                  Ex: CpG,CpH,CHH, etc}
#'                  \item{\code{treatment}:}{treatment vector denoting which 
#'                  samples are test and control}
#'                  \item{\code{coverage.index}:}{vector denoting which 
#'                  columns in the data correspond to 
#'                  coverage values}
#'                  \item{\code{numCs.index}:}{vector denoting which columns 
#'                  in the data correspond to 
#'                  number of methylatedCs values}
#'                  \item{\code{numTs.index}:}{vector denoting which columns
#'                   in the data correspond to 
#'                  number of unmethylated Cs values}
#'                  \item{\code{destranded}:}{ logical value. 
#'                  If \code{TRUE} object is destranded, 
#'                  if \code{FALSE} it is not.}
#'                  \item{\code{resolution}:}{ resolution of methylation 
#'                  information, 
#'                  allowed values: 'base' or 'region'}
#'                  \item{\code{dbtype}:}{string for type of the flat file 
#'                  database, ex: tabix}
#' }
#' 
#' @section Details:
#' \code{methylBaseDB} class has the same functionality as
#'  \code{\link{methylBase}} class, 
#' but the data is saved in a flat database file and therefore allocates
#'  less space in memory.
#' 
#' 
#' @section Subsetting:
#'  In the following code snippets, \code{x} is a \code{methylBaseDB}.
#'  Subsetting by \code{x[i,]} will produce a new \code{methylBase} object 
#'  if subsetting is done on
#'  rows. Column subsetting is not directly allowed to prevent errors in the 
#'  downstream analysis. see ?methylKit[ .
#' 
#' @section Accessors: 
#' The following functions provides access to data slots of methylDiffDB:
#' - \code{\link{getData}}: get the data slot from the methylKit objects,
#' - \code{\link{getAssembly}}: get assembly of the genome,
#' - \code{\link{getContext}}: get the context of methylation
#' 
#' 
#' @section Coercion:
#'   \code{methylBaseDB} object can be coerced to:
#'   \code{\link[GenomicRanges:GRanges-class]{GRanges}} object via \code{\link{as}} function.
#'   \code{\link{methylBase}} object via \code{\link{as}} function. 
#' 
#' @examples
#' data(methylKit)
#' methylBaseDB.obj <- unite(methylRawList.obj,save.db=TRUE,dbdir="methylDB")
#' library(GenomicRanges)
#' my.gr=as(methylBaseDB.obj,"GRanges")
#' 
#' 
#' # remove Database again
#' rm(methylBaseDB.obj)
#' unlink("methylDB",recursive=TRUE)
#' 
#' 
#' @name methylBaseDB-class
#' @aliases methylBaseDB
#' @docType class
#' @rdname methylBaseDB-class
#' @export
setClass("methylBaseDB",slots=list(dbpath = "character", 
                                   num.records = "numeric", 
                                   sample.ids = "character",
                                   assembly = "character", 
                                   context = "character", resolution = "character",
                                   dbtype = "character", 
                                   treatment = "numeric", coverage.index = "numeric",
                                   numCs.index = "numeric", 
                                   numTs.index = "numeric", 
                                   destranded = "logical"),
         validity = valid.methylBaseDB)

# PRIVATE function:
# makes a methylBaseDB object from a df
# it is called from read function or whenever this functionality is needed
makeMethylBaseDB<-function(df,dbpath,dbtype,
                           sample.ids, assembly ,context,
                           resolution,treatment,coverage.index,
                           numCs.index,numTs.index,destranded,
                           suffix=NULL
                           )
  {
  # new tabix file is named by "metyhlBase"+tmpstring, works for now
  # if additional suffix is passed, tmpstring is skipped 
  filepath <- paste0(ifelse(is.null(suffix),
                            yes = tempfile(pattern = "methylBase_",tmpdir = dbpath),
                            no = paste0(dbpath,"/methylBase",suffix)),
                     ".txt")
  
  # new tabix file is named by concatenation of sample.ids, works for now
  # filepath=paste0(dbpath,"/",paste0(sample.ids,collapse = "_"),suffix,".txt") 
  #print(filepath)
  df <- df[with(df,order(chr,start,end)),]
  df2tabix(df,filepath)
  num.records=Rsamtools::countTabix(paste0(filepath,".bgz"))[[1]] ## 

  new("methylBaseDB",dbpath=paste0(filepath,".bgz"),num.records=num.records,
      sample.ids = sample.ids, assembly = assembly,context=context,
      resolution=resolution,dbtype=dbtype,treatment=treatment,
      coverage.index=coverage.index,numCs.index=numCs.index,numTs.index=numTs.index,
      destranded=destranded)
}

# PRIVATE function:
# reads a methylBaseDB object from flat file database
# it is called from read function or whenever this functionality is needed
readMethylBaseDB<-function(dbpath,dbtype,
                           sample.ids, assembly ,context,
                           resolution,treatment,destranded,skip=0){
  
  if(!file.exists(paste0(dbpath,".tbi"))) 
  {  
    Rsamtools::indexTabix(dbpath,seq=1, start=2, end=3,
                          skip=skip, comment="#", zeroBased=FALSE)
  }
  num.records=Rsamtools::countTabix(dbpath)[[1]] ## 
  
  # determine postion of coverage/numCs/numTs indices
  df <- headTabix(dbpath)
  numsamples = (length(df)-4)/3
  coverage.ind=seq(5,by=3,length.out=numsamples)
  numCs.ind   =coverage.ind+1
  numTs.ind   =coverage.ind+2
  
  obj <- new("methylBaseDB",dbpath=normalizePath(dbpath),num.records=num.records,
              sample.ids = sample.ids, assembly = assembly,context=context,
              resolution=resolution,dbtype=dbtype,treatment=treatment,
              coverage.index=coverage.ind,numCs.index=numCs.ind,numTs.index=numTs.ind,
              destranded=destranded)
 
  if(valid.methylBaseDB(obj)) {return(obj)}   
    
}

# PRIVATE function:
# creates a methylBaseDB object from a flat file database and reads header content
# it is called from read function or whenever this functionality is needed
loadMethylBaseDB<-function(dbpath,skip=0){
  
  if(!file.exists(paste0(dbpath,".tbi"))) 
  {  
    Rsamtools::indexTabix(dbpath,seq=1, start=2, end=3,
                          skip=skip, comment="#", zeroBased=FALSE)
  }
  
  h <- readTabixHeader(dbpath)
  
  if(! any(h$class == c("methylBase", "methylBaseDB") ) ) {
    stop("Tabix file does not originate from methylBase or methylBaseDB.\nPlease provide the correct class of data.")
  }
  
  obj <- new("methylBaseDB",dbpath=normalizePath(dbpath),num.records=h$num.records,
             sample.ids = h$sample.ids, assembly = h$assembly,context=h$context,
             resolution=h$resolution,dbtype=h$dbtype,treatment=h$treatment,
             coverage.index=h$coverage.ind,numCs.index=h$numCs.ind,numTs.index=h$numTs.ind,
             destranded=h$destranded)
  
  if(valid.methylBaseDB(obj)) {return(obj)}    
}

# methylDiffDB -------------------------------------------------------


valid.methylDiffDB <- function(object) {
  
  
  check1=( (object@resolution == "base") | (object@resolution == "region") )
  check2=file.exists(object@dbpath)
  check3=( length(object@sample.ids) == length(object@treatment) )
  
  if(check1 & check2 & check3 ){
    return(TRUE)
  }
  else if (! check1 ){
    cat("resolution slot has to be either 'base' or 'region':",
            "other values not allowed")
    FALSE
  }
  else if(! check2){
    cat("The DB file can not be found, check the value of 'dbpath'")
    FALSE
  }
  else if(! check3){
    cat("The number of samples is different from the number of treatments, ",
            "check the length of 'treatment'")
    FALSE
  }
  
}

#' An S4 class that holds differential methylation information as flat file database
#'
#' This class is designed to hold statistics and locations for differentially 
#' methylated regions/bases as flat file database.
#' \code{\link{calculateDiffMeth}} function returns an object 
#' with \code{methylDiffDB} class.
#'          
#' @section Slots:\describe{
#'    \item{\code{dbpath}:}{path to flat file database(s) }
#'    \item{\code{num.records}:}{number of records (lines) in the object}
#'    \item{\code{sample.ids}}{ids/names of samples in a vector}
#'    \item{\code{assembly}}{a name of genome assembly, such as :hg18,mm9, etc}
#'    \item{\code{context}}{numeric vector identifying which samples are which
#'    group }
#'    \item{\code{treatment}}{numeric vector identifying which samples are which
#'     group }
#'    \item{\code{destranded}}{logical denoting if methylation inormation is
#'     destranded or not}
#'    \item{\code{resolution}}{string either 'base' or 'region' defining the 
#'    resolution of methylation information}
#'    \item{\code{dbtype}:}{string for type of the flat file database, ex: tabix}
#'
#' }
#' 
#' @section Details:
#' \code{methylDiffDB} class has the same functionality as 
#' \code{\link{methylDiff}} class, 
#' but the data is saved in a flat database file and therefore 
#' allocates less space in memory.
#' 
#' 
#' @section Subsetting:
#'  In the following code snippets, \code{x} is a \code{methylDiffDB}.
#'  Subsetting by \code{x[i,]} will produce a new object if subsetting is done 
#'  on rows. Column subsetting is not directly allowed to prevent errors in the 
#'  downstream analysis. see ?methylKit[ .
#' 
#' @section Coercion:
#'   \code{methylDiffDB} object can be coerced to:
#'   \code{\link[GenomicRanges:GRanges-class]{GRanges}} object via \code{\link{as}} function.
#'   \code{\link{methylDiff}} object via \code{\link{as}} function.
#'    
#' @section Accessors: 
#' The following functions provides access to data slots of methylDiffDB:
#' - \code{\link{getData}}: get the data slot from the methylKit objects,
#' - \code{\link{getAssembly}}: get assembly of the genome,
#' - \code{\link{getContext}}: get the context of methylation
#' 
#' @examples
#' data(methylKit)
#' 
#' methylDiffDB.obj <- calculateDiffMeth(methylBase.obj,save.db=TRUE,dbdir="methylDB")
#' 
#' library(GenomicRanges)
#' my.gr=as(methylDiffDB.obj,"GRanges")
#' 
#' # remove Database again
#' rm(methylDiffDB.obj)
#' unlink("methylDB",recursive=TRUE)
#' 
#' @name methylDiffDB-class
#' @aliases methylDiffDB
#' @rdname methylDiffDB-class
#' @export
#' @docType class
setClass("methylDiffDB",slots = list(dbpath= "character",num.records= "numeric",sample.ids = "character", 
                                     assembly = "character",context = "character",dbtype = "character", 
                                     treatment="numeric",destranded="logical",resolution="character"),validity = valid.methylDiffDB)

# PRIVATE function:
# makes a methylDiffDB object from a df
# it is called from read function or whenever this functionality is needed
makeMethylDiffDB<-function(df,dbpath,dbtype,
                           sample.ids, assembly ,context,
                           resolution,treatment,destranded,
                           suffix=NULL){
  
  # new tabix file is named by "metyhlBase"+tmpstring, works for now
  # if additional suffix is passed, tmpstring is skipped 
  filepath <- paste0(ifelse(is.null(suffix),
                            yes = tempfile(pattern = "methylDiff_",tmpdir = dbpath),
                            no = paste0(dbpath,"/methylDiff",suffix)),
                     ".txt")
  
  # # new tabix file is named by concatenation of sample.ids, works for now
  # filepath=paste0(dbpath,"/",paste0(sample.ids,collapse = "_"),suffix,".txt")
  df <- df[with(df,order(chr,start,end)),]
  df2tabix(df,filepath)
  num.records=Rsamtools::countTabix(paste0(filepath,".bgz"))[[1]] ## 

  new("methylDiffDB",dbpath=paste0(filepath,".bgz"),num.records=num.records,
      sample.ids = sample.ids, assembly = assembly,context=context,
      resolution=resolution,dbtype=dbtype,treatment=treatment,
      destranded=destranded)
}

# PRIVATE function:
# reads a methylDiffDB object from flat file database
# it is called from read function or whenever this functionality is needed
readMethylDiffDB<-function(dbpath,dbtype,
                           sample.ids, assembly ,context,
                           resolution,treatment,destranded,skip=0){
  
  if(!file.exists(paste0(dbpath,".tbi"))) 
  {  
    Rsamtools::indexTabix(dbpath,seq=1, start=2, end=3,
                          skip=skip, comment="#", zeroBased=FALSE)
  }
  num.records=Rsamtools::countTabix(dbpath)[[1]] ## 
  
  obj <- new("methylDiffDB",dbpath=normalizePath(dbpath),num.records=num.records,
      sample.ids = sample.ids, assembly = assembly,context=context,
      resolution=resolution,dbtype=dbtype,treatment=treatment,
      destranded=destranded)
  
  if(valid.methylDiffDB(obj)) {return(obj)} 
}

# PRIVATE function:
# creates a methylDiffDB object from a flat file database and reads header content
# it is called from read function or whenever this functionality is needed
loadMethylDiffDB<-function(dbpath,skip=0){
  
  if(!file.exists(paste0(dbpath,".tbi"))) 
  {  
    Rsamtools::indexTabix(dbpath,seq=1, start=2, end=3,
                          skip=skip, comment="#", zeroBased=FALSE)
  }
  
  h <- readTabixHeader(dbpath)
  
  if(! any(h$class == c("methylDiff", "methylDiffDB") ) ) {
    stop("Tabix file does not originate from methylDiff or methylDiffDB.\nPlease provide the correct class of data.")
  }
  
  obj <-   new("methylDiffDB",dbpath=normalizePath(dbpath),num.records=h$num.records,
               sample.ids = h$sample.ids, assembly = h$assembly,context=h$context,
               resolution=h$resolution,dbtype=h$dbtype,treatment=h$treatment,
               destranded=h$destranded)
  
  if(valid.methylDiffDB(obj)) {return(obj)}    
}

# coercion functions ------------------------------------------------------

setAs("methylRawDB", "GRanges", function(from)
{
  gr <- headTabix(tbxFile = from@dbpath, nrow = from@num.records, return.type = "GRanges")
  names(GenomicRanges::mcols(gr)) <- c("coverage","numCs","numTs") 
  return(gr)
})

setAs("methylBaseDB", "GRanges", function(from)
{
  from=getData(from)
  GRanges(seqnames=as.character(from$chr),ranges=IRanges(start=from$start, end=from$end),
          strand=from$strand, 
          data.frame(from[,5:ncol(from)])
  )
#   gr <- headTabix(tbxFile = from@dbpath, nrow = from@num.records, return.type = "GRanges")
#   names(GenomicRanges::mcols(gr)) <- c("coverage","numCs","numTs") 
#   return(gr)
})

setAs("methylDiffDB", "GRanges", function(from)
{
  gr <- headTabix(tbxFile = from@dbpath, nrow = from@num.records, 
                  return.type = "GRanges")
  names(GenomicRanges::mcols(gr)) <- c("pvalue","qvalue","meth.diff") 
  return(gr)
})

## coerce methylDB to methyl-obj
setAs("methylRawDB","methylRaw", function(from)
{
  return(from[])
})

setAs("methylRawListDB","methylRawList", function(from)
{
  outList = lapply(from,as,"methylRaw")
  new("methylRawList", outList,treatment=from@treatment)
})

setAs("methylBaseDB","methylBase", function(from)
{
  return(from[])
})

setAs("methylDiffDB","methylDiff", function(from)
{
  return(from[])
})



## coerce methyl-obj to methylDB

#' coerce methylKit objects from memory to flat file database objects
#' 
#' The function converts in-memory methylKit objects to methylDB objects
#' 
#' @param obj an \code{\link{methylBase}},
#' \code{\link{methylRaw}},
#' \code{\link{methylRawList}} or an \code{\link{methylDiff}} object
#' @param dbdir directory where flat file database(s) should be stored, 
#'              defaults to getwd(), working directory.
#' @examples 
#' \dontrun{
#' data(methylKit)
#' 
#' 
#' makeMethylDB(methylBase.obj,"my/path")
#' }
#' 
#' @return an \code{\link{methylBaseDB}},\code{\link{methylRawDB}},
#' \code{\link{methylRawListDB}} or an \code{\link{methylDiffDB}} object
#' @export
#' @docType methods
#' @rdname makeMethylDB-methods
setGeneric("makeMethylDB", def=function(obj,dbdir=getwd()) 
  standardGeneric("makeMethylDB"))

#' @rdname makeMethylDB-methods
#' @aliases makeMethylDB,methylBase-methods
setMethod("makeMethylDB", signature="methylBase", definition=function(obj,dbdir) 
  {
  dbdir <- .check.dbdir(dbdir)
  objdb <- makeMethylBaseDB(df=getData(obj),dbpath=dbdir,dbtype="tabix",
                   sample.ids=obj@sample.ids,
                   assembly=obj@assembly,context=obj@context,
                   treatment=obj@treatment,coverage.index=obj@coverage.index,
                   numCs.index=obj@numCs.index,numTs.index=obj@numTs.index,
                   destranded=obj@destranded,
                   resolution=obj@resolution)
  message(paste0("flatfile located at: ",getDBPath(objdb)))
  
  return(objdb)
})

#' @rdname makeMethylDB-methods
#' @aliases makeMethylDB,methylRaw-methods
setMethod("makeMethylDB", signature="methylRaw", definition=function(obj,dbdir) {
  dbdir <- .check.dbdir(dbdir)
  objdb <- makeMethylRawDB(df=getData(obj), dbpath=dbdir, dbtype="tabix", 
                  sample.id=obj@sample.id,
                  assembly=obj@assembly, context=obj@context, 
                  resolution=obj@resolution)
  message(paste0("flatfile located at: ",getDBPath(objdb)))
  return(objdb)
})

#' @rdname makeMethylDB-methods
#' @aliases makeMethylDB,methylDiff-methods
setMethod("makeMethylDB", signature="methylDiff", 
          definition=function(obj,dbdir) {
  dbdir <- .check.dbdir(dbdir)
  suffix <- "_diffMeth"
  objdb <- makeMethylDiffDB(df=getData(obj), dbpath=dbdir, dbtype="tabix", 
                   sample.ids=obj@sample.ids,
                   assembly=obj@assembly,context=obj@context,
                   destranded=obj@destranded,treatment=obj@treatment,
                   resolution=obj@resolution,suffix=suffix)
  message(paste0("flatfile located at: ",getDBPath(objdb)))
  return(objdb)
})

#' @rdname makeMethylDB-methods
#' @aliases makeMethylDB,methylRawList-methods
setMethod("makeMethylDB", signature="methylRawList", 
          definition=function(obj,dbdir) {
  dbdir <- .check.dbdir(dbdir)
  outList <- lapply(obj,makeMethylDB,dbdir)
  objdb <- new("methylRawListDB",outList,treatment=obj@treatment)
  return(objdb)
})

# accessors ---------------------------------------------------------------


#' @rdname getAssembly-methods
#' @aliases getAssembly,methylRawDB-method
setMethod("getAssembly", signature="methylRawDB", definition=function(x) {
  return(x@assembly)
})

#' @rdname getAssembly-methods
#' @aliases getAssembly,methylBaseDB-method
setMethod("getAssembly", signature="methylBaseDB", definition=function(x) {
  return(x@assembly)
}) 

#' @rdname getAssembly-methods
#' @aliases getAssembly,methylDiffDB-method
setMethod("getAssembly", signature="methylDiffDB", definition=function(x) {
  return(x@assembly)
})

#' @rdname getContext-methods
#' @aliases getContext,methylRawDB-method
setMethod("getContext", signature="methylRawDB", definition=function(x) {
  return(x@context)
})

#' @rdname getContext-methods
#' @aliases getContext,methylBaseDB-method
setMethod("getContext", signature="methylBaseDB", definition=function(x) {
  return(x@context)
})

#' @rdname getContext-methods
#' @aliases getContext,methylDiffDB-method
setMethod("getContext", signature="methylDiffDB", definition=function(x) {
  return(x@context)
})

#' @rdname getData-methods
#' @aliases getData,methylRawDB-method
setMethod("getData", signature="methylRawDB", definition=function(x) {
  df <- headTabix(tbxFile = x@dbpath, nrow = x@num.records, 
                  return.type = "data.frame")
  .setMethylDBNames(df,"methylRawDB")
  
  return(df)
})

#' @rdname getData-methods
#' @aliases getData,methylBaseDB-method
setMethod("getData", signature="methylBaseDB", definition=function(x) {
  df <- headTabix(tbxFile = x@dbpath, nrow = x@num.records, 
                  return.type = "data.frame")
 .setMethylDBNames(df,"methylBaseDB")
  
  return(df)
})

#' @rdname getData-methods
#' @aliases getData,methylDiffDB-method
setMethod(f="getData", signature="methylDiffDB", definition=function(x) {
  df <- headTabix(tbxFile = x@dbpath, nrow = x@num.records, 
                  return.type = "data.frame")
  .setMethylDBNames(df,"methylDiffDB")
  
  return(df)
})

# show functions ----------------------------------------------------------

#' @rdname show-methods
#' @aliases show,methylRawDB
setMethod("show", "methylRawDB", function(object) {
  
  cat("methylRawDB object with",object@num.records,"rows\n--------------\n")
  print(.setMethylDBNames(headTabix(object@dbpath,nrow = 6,
                                    return.type = "data.frame"),
                          methylDBclass = "methylRawDB"))
  cat("--------------\n")
  cat("sample.id:",object@sample.id,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("resolution:", object@resolution,"\n")
  cat("dbtype:", object@dbtype,"\n")
  #cat("dbpath:",object@dbpath,"\n")
  cat("\n")
})

#' @rdname show-methods
#' @aliases show,methylRawListDB
setMethod("show", "methylRawListDB", function(object) {
  
  cat("methylRawListDB object with",length(object),"methylRawDB objects\n\n")
  
  lapply(object,show)
  cat("treament:", object@treatment,"\n")
  
})

#' @rdname show-methods
#' @aliases show,methylBaseDB
setMethod("show", "methylBaseDB", function(object) {
  
  cat("methylBaseDB object with",object@num.records,"rows\n--------------\n")
  print(.setMethylDBNames(headTabix(object@dbpath,
                                    nrow = 6,return.type = "data.frame"),
                          methylDBclass = "methylBaseDB"))
  cat("--------------\n")
  cat("sample.ids:",object@sample.ids,"\n")
  cat("destranded",object@destranded,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("treament:", object@treatment,"\n")
  cat("resolution:", object@resolution,"\n")
  cat("dbtype:", object@dbtype,"\n")
  #cat("dbpath:",object@dbpath,"\n")
  cat("\n")
  
})

#' @rdname show-methods
#' @aliases show,methylDiffDB
setMethod("show", "methylDiffDB", function(object) {
  
  cat("methylDiffDB object with",object@num.records,"rows\n--------------\n")
  print(.setMethylDBNames(headTabix(object@dbpath,nrow = 6,
                                    return.type = "data.frame"),
                          methylDBclass = "methylDiffDB"))
  cat("--------------\n")
  cat("sample.ids:",object@sample.ids,"\n")
  cat("destranded",object@destranded,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("treament:", object@treatment,"\n")
  cat("resolution:", object@resolution,"\n")
  cat("dbtype:", object@dbtype,"\n")
  #cat("dbpath:",object@dbpath,"\n")
  cat("\n")
})

# subset classes ----------------------------------------------------------

#' @aliases select,methylRawDB-method
#' @rdname select-methods
setMethod("select", "methylRawDB",
          function(x, i)
          {
            if(missing(i)) { i <- 1:x@num.records }
            df <- headTabix(x@dbpath,nrow = max(i),return.type = "data.frame")
            .setMethylDBNames(df,"methylRawDB")
            
            if( max(i) > x@num.records  )
              stop("subscript contains out-of-bounds indices")
            
            new("methylRaw",df[i,],
                sample.id=x@sample.id,
                assembly=x@assembly,
                context=x@context,
                resolution=x@resolution)
          }
)


#' @aliases select,methylBaseDB-method
#' @rdname select-methods
setMethod("select", "methylBaseDB",
          function(x, i)
          {
            if(missing(i)) { i <- 1:x@num.records }
            df <- headTabix(x@dbpath,nrow = max(i),return.type = "data.frame")
            .setMethylDBNames(df,"methylBaseDB")
            
            if( max(i) > x@num.records  )
              stop("subscript contains out-of-bounds indices")
            
            new("methylBase",df[i,],
                sample.ids = x@sample.ids, 
                assembly = x@assembly,
                context = x@context,
                treatment=x@treatment,
                coverage.index=x@coverage.index,
                numCs.index=x@numCs.index,
                numTs.index=x@numTs.index,
                destranded=x@destranded,
                resolution =x@resolution)
          }
)

#' @aliases select,methylDiffDB-method
#' @rdname select-methods
setMethod("select", "methylDiffDB",
          function(x, i)
          {
            if(missing(i)) { i <- 1:x@num.records }
            df <- headTabix(x@dbpath,nrow = max(i),return.type = "data.frame")
            .setMethylDBNames(df,"methylDiffDB")
            
            if( max(i) > x@num.records  )
              stop("subscript contains out-of-bounds indices")
            
            new("methylDiff",df[i,],
                sample.ids = x@sample.ids,
                assembly = x@assembly,
                context = x@context,
                treatment=x@treatment,
                destranded=x@destranded,
                resolution=x@resolution)
          }
)


#' @aliases [,methylRawDB,ANY,ANY,ANY-method
#' @rdname extract-methods
#' @aliases extract,methylRawDB,ANY-method
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

#' @aliases [,methylBaseDB,ANY,ANY,ANY-method
#' @aliases extract,methylBaseDB,ANY-method
#' @rdname extract-methods
setMethod("[",signature(x="methylBaseDB", i = "ANY", j="ANY"), 
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

#' @rdname extract-methods
#' @aliases [,methylDiffDB,ANY,ANY,ANY-method
#' @aliases extract,methylDiffDB,ANY-method
setMethod("[",signature(x="methylDiffDB", i="ANY", j="ANY"),
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


# select by range ---------------------------------------------------------



#' @aliases selectByOverlap,methylRawDB-method
#' @rdname selectByOverlap-methods
setMethod("selectByOverlap", c("methylRawDB","GRanges"),
          function(object, ranges){
            
    if(missing(ranges) | class(ranges)!="GRanges") {
      stop("No ranges specified or given ranges object not of class GRanges, ", 
           "please check your input!")
      }
            
            df <-  getTabixByOverlap(tbxFile = object@dbpath,granges = ranges, 
                                     return.type = "data.frame") 
            
            

            obj <- new("methylRaw",.setMethylDBNames(df,"methylRawDB"),
                       sample.id=object@sample.id,
                       assembly=object@assembly,
                       context=object@context,
                       resolution=object@resolution)
            obj
          }
)


#' @aliases selectByOverlap,methylRawListDB-method
#' @rdname selectByOverlap-methods
setMethod("selectByOverlap", c("methylRawListDB","GRanges"),
          function(object, ranges){
            
    if(missing(ranges) | class(ranges)!="GRanges") {
      stop("No ranges specified or given ranges object not of class GRanges, ", 
           "please check your input!")
    }
    
    new.list <- lapply(object,selectByOverlap,ranges)
    
    new("methylRawList",new.list,treatment=object@treatment)
      
  }
)

#' @aliases selectByOverlap,methylBaseDB-method
#' @rdname selectByOverlap-methods
setMethod("selectByOverlap", c("methylBaseDB","GRanges"),
          function(object, ranges){
  
if(missing(ranges) | class(ranges)!="GRanges") {
  stop("No ranges specified or given ranges object not of class GRanges, ",
       "please check your input!")
}

df <-  getTabixByOverlap(tbxFile = object@dbpath,granges = ranges, 
                         return.type = "data.frame") 

new("methylBase",.setMethylDBNames(df,"methylBaseDB"),
    sample.ids = object@sample.ids, 
    assembly = object@assembly,
    context = object@context,
    treatment=object@treatment,
    coverage.index=object@coverage.index,
    numCs.index=object@numCs.index,
    numTs.index=object@numTs.index,
    destranded=object@destranded,
    resolution =object@resolution)

}
)

#' @aliases selectByOverlap,methylDiffDB-method
#' @rdname selectByOverlap-methods
setMethod("selectByOverlap", c("methylDiffDB","GRanges"),
          function(object, ranges){
            
  if(missing(ranges) | class(ranges)!="GRanges") {
    stop("No ranges specified or given ranges object not of class GRanges, ",
         "please check your input!")
  }
  
  df <-  getTabixByOverlap(tbxFile = object@dbpath,
                           granges = ranges, return.type = "data.frame") 
  
  new("methylDiff",.setMethylDBNames(df,"methylDiffDB"),
      sample.ids = object@sample.ids,
      assembly = object@assembly,
      context = object@context,
      treatment=object@treatment,
      destranded=object@destranded,
      resolution=object@resolution)
  
}
)


# get/set values ----------------------------------------------------------


#  Get or Set treatment vector of the methylRawListDB or methylBaseDB object
#  
#  The function returns or replaces the treatment vector stored in any of the
#  \code{\link{methylBaseDB}},\code{\link{methylRawListDB}} objects
#  
#  @param x an \code{\link{methylBaseDB}} or \code{\link{methylRawListDB}}
#    object
#  @param value a valid replacement for the treatment vector of the object
#  @usage getTreatment(x) getTreatment(x) <- value
#  @examples
#  
#  data(methylKit)
#  
#  # The treatment vector can be printed ..
#  getTreatment(methylBaseDB.obj)
#   
#  # .. or replaced with a new one  
#  newObj <- methylBaseDB.obj
#  getTreatment(newObj) <- c(1,2,3,4)
#  getTreatment(newObj)


#' @rdname getTreatment-methods
#' @aliases getTreatment,methylRawListDB-method
setMethod("getTreatment", signature = "methylRawListDB", function(x) {
    return(x@treatment)
})

#' @rdname getTreatment-methods
#' @aliases getTreatment,methylRawListDB-method
setReplaceMethod("getTreatment", signature = "methylRawListDB", function(x, value) {
  
  if(! ( length(x@treatment) == length(value) ) ){
    stop("The new treatment vector is not valid, check the length of input")
  } else {
    x@treatment <- value
    return(x)
  }
})


#' @rdname getTreatment-methods
#' @aliases getTreatment,methylBaseDB-method
setMethod("getTreatment", signature = "methylBaseDB", function(x) {
  return(x@treatment)
})

#' @rdname getTreatment-methods
#' @aliases getTreatment,methylBaseDB-method
setReplaceMethod("getTreatment", signature = "methylBaseDB", function(x, value) {
  
  if(! ( length(x@treatment) == length(value) ) ){
    stop("The new treatment vector is not valid, check the length of input")
  } else {
    x@treatment <- value
    return(x)
  }
  
})


#' @rdname getTreatment-methods
#' @aliases getTreatment,methylDiffDB-method
setMethod("getTreatment", signature = "methylDiffDB", function(x) {
  return(x@treatment)
})

#' @rdname getTreatment-methods
#' @aliases getTreatment,getTreatment,methylDiffDB-method
setReplaceMethod("getTreatment", signature = "methylDiffDB", function(x, value) {
  
  if(! ( length(x@treatment) == length(value) ) ){
    stop("The new treatment vector is not valid, check the length of input")
  } else {
    x@treatment <- value
    return(x)
  }
  
})






#  Get or Set sample ids of the methylDB objects
#  
#  The function returns or replaces the sample-ids stored in any of the 
# \code{\link{methylRawDB}}, \code{\link{methylBaseDB}},
# \code{\link{methylRawListDB}} or \code{\link{methylDiffDB}} objects 
#  
#  @param x an \code{\link{methylBaseDB}},\code{\link{methylRawListDB}} or
#           \code{\link{methylDiffDB}} object
#  @param value a valid replacement for the sample-ids of the object 
#  @usage 
#  getSampleID(x)
#  getSampleID(x) <- value
#  @examples
#  
#  data(methylKit)
#  
#  #The Sample-Ids can be printed ..
#  getSampleID(methylBaseDB.obj)
#  
#  # .. or replaced. 
#  newObj <- methylBaseDB.obj
#  getSampleID(newObj) <- c("sample1","sample2","sample3","sample4")
#  getSampleID(newObj)
#  

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylRawListDB-method
setMethod("getSampleID", signature = "methylRawListDB", function(x) {
  names <- vapply(x,function(z) z@sample.id,FUN.VALUE = "character")
  return(names)
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylRawListDB-method
setReplaceMethod("getSampleID", signature = "methylRawListDB", function(x, value) {
  
  if(! ( length(getSampleID(x)) == length(value) ) ){
    stop("The vector of new sample ids is not valid, check the length of input")
  } else {
    treatment <- x@treatment
    x <- mapply(`getSampleID<-`, x, value)
    x <- new("methylRawListDB",x,treatment=treatment)
    return(x)
  }
  
})


#' @rdname getSampleID-methods
#' @aliases getSampleID,methylBaseDB-method
setMethod("getSampleID", signature = "methylBaseDB", function(x) {
  return(x@sample.ids)
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylBaseDB-method
setReplaceMethod("getSampleID", signature = "methylBaseDB", function(x, value) {
  
  if(! ( length(x@sample.ids) == length(value) ) ){
    stop("The vector of new sample ids is not valid, check the length of input")
  } else {
    x@sample.ids <- value
    return(x)
  }
  
})


#' @rdname getSampleID-methods
#' @aliases getSampleID,methylRawDB-method
setMethod("getSampleID", signature = "methylRawDB", function(x) {
  return(x@sample.id)
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylRawDB-method
setReplaceMethod("getSampleID", signature = "methylRawDB", function(x, value) {
  
  if(! ( length(value) == 1 ) ){
    stop("The vector of new sample ids is not valid, check the length of input")
  } else {
    x@sample.id <- value
    return(x)
  }
  
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylDiffDB-method
setMethod("getSampleID", signature = "methylDiffDB", function(x) {
  return(x@sample.id)
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylDiffDB-method
setReplaceMethod("getSampleID", signature = "methylDiffDB", function(x, value) {
  
  if(! ( length(value) == length(x@sample.ids) ) ){
    stop("The vector of new sample ids is not valid, check the length of input")
  } else {
    x@sample.ids <- value
    return(x)
  }
  
})


#' Get path to database of the methylDB objects
#' 
#' The function returns the path to the flat file database that stores the data 
#' of the \code{\link{methylRawDB}}, 
#' \code{\link{methylRawListDB}}, \code{\link{methylBaseDB}} or 
#' \code{\link{methylDiffDB}} objects. 
#'  
#' 
#' @param x an \code{\link{methylBaseDB}},\code{\link{methylRawDB}},
#' \code{\link{methylRawListDB}} or \code{\link{methylDiffDB}} object
#' @examples
#' 
#' data(methylKit)
#' 
#' methylBaseDB.obj <- unite(methylRawList.obj,save.db=TRUE,dbdir="methylDB")
#' 
#' 
#' #The path to the database is returned
#' getDBPath(methylBaseDB.obj)
#' 
#' 
#' # remove Database again
#' rm(methylBaseDB.obj)
#' unlink("methylDB",recursive=TRUE)
#' 
#' @export
#' @docType methods
#' @rdname getDBPath-methods
setGeneric("getDBPath", def=function(x) standardGeneric("getDBPath"))

#' @rdname getDBPath-methods
#' @aliases getDBPath,methylRawListDB-method
setMethod("getDBPath", signature = "methylRawListDB", function(x) {
  names <- vapply(x,function(z) z@dbpath,FUN.VALUE = "character")
  return(names)
})

#' @rdname getDBPath-methods
#' @aliases getDBPath,methylBaseDB-method
setMethod("getDBPath", signature = "methylBaseDB", function(x) {
  return(x@dbpath)
})

#' @rdname getDBPath-methods
#' @aliases getDBPath,methylRawDB-method
setMethod("getDBPath", signature = "methylRawDB", function(x) {
  return(x@dbpath)
})

#' @rdname getDBPath-methods
#' @aliases getDBPath,methylDiffDB-method
setMethod("getDBPath", signature = "methylDiffDB", function(x) {
  return(x@dbpath)
})


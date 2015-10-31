
#---------------------------------------------------------------------------------------
# regular R functions to be used in S4 functions


#' set column names for methylRawDB and methylBaseDB data aquired from flat file database
#' @param df data.frame containing methylRaw or methylBase data
#' @param methylDBclass 
.setMethylDBNames <- function(df,methylDBclass=c("methylDB","methylBaseDB","methylDiffDB")){
  
  if(missing(methylDBclass)){

    if( length(df) == 7 & unique(sapply(df[,5:7],class))=="integer"){
      
      data.table::setnames(x = df,old = names(df), new = c("chr","start","end","strand","coverage","numCs","numTs"))
      
    } else if( length(df) == 7 & unique(sapply(df[,5:7],class))=="numeric"){
      
      data.table::setnames(x = df,old = names(df), new = c("chr","start","end","strand","pvalue","qvalue","meth.diff")) 

    } else if( length(df) > 7){
      
      data.table::setnames(x = df,old = names(df)[1:4], new = c("chr","start","end","strand"))
      # get indices of coverage,numCs and numTs in the data frame 
      numsamples = (length(df)-4)/3
      coverage.ind=seq(5,by=3,length.out=numsamples)
      numCs.ind   =coverage.ind+1
      numTs.ind   =coverage.ind+2
      
      # change column names
      data.table::setnames(df,names(df)[coverage.ind], paste(c("coverage"),1:numsamples,sep="" ))
      data.table::setnames(df,names(df)[numCs.ind], paste(c("numCs"),1:numsamples,sep="" ))
      data.table::setnames(df,names(df)[numTs.ind], paste(c("numTs"),1:numsamples,sep="" ))
      
    } 
    
    #return(df)
    
  } else {
    
    if( methylDBclass == "methylRawDB" ){
    
      data.table::setnames(x = df,old = names(df), new = c("chr","start","end","strand","coverage","numCs","numTs"))
    
    } else if ( methylDBclass == "methylBaseDB"){
    
      data.table::setnames(x = df,old = names(df)[1:4], new = c("chr","start","end","strand"))
      # get indices of coverage,numCs and numTs in the data frame 
      numsamples = (length(df)-4)/3
      coverage.ind=seq(5,by=3,length.out=numsamples)
      numCs.ind   =coverage.ind+1
      numTs.ind   =coverage.ind+2
      
      # change column names
      data.table::setnames(df,names(df)[coverage.ind], paste(c("coverage"),1:numsamples,sep="" ))
      data.table::setnames(df,names(df)[numCs.ind], paste(c("numCs"),1:numsamples,sep="" ))
      data.table::setnames(df,names(df)[numTs.ind], paste(c("numTs"),1:numsamples,sep="" ))
      
    } else if( methylDBclass == "methylDiffDB" ){
      
      data.table::setnames(x = df,old = names(df), new = c("chr","start","end","strand","pvalue","qvalue","meth.diff"))
    
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


#' An S4 class for storing raw methylation data as flat file database.
#'
#' This object stores the raw mehylation data that is read in through read 
#' function as flat file database.The raw methylation data is basically
#' percent methylation values and read coverage values per genomic base/region.
#'
#' @section Slots:\describe{
#'  \item{\code{dbpath}:}{path to flat file database(s) }
#'  \item{\code{num.records}:}{number of records (lines) in the object}
#'  \item{\code{sample.id}:}{string for an identifier of the sample}
#'  \item{\code{assembly}:}{string for genome assembly, ex: hg18,hg19,mm9}
#'  \item{\code{context}:}{ methylation context string, ex: CpG,CpH,CHH, etc.}
#'  \item{\code{resolution}:}{ resolution of methylation information, 'base' or 
#'  'region'}
#'  \item{\code{dbtype}:}{string for type of the flat file database, ex: tabix}
#'                 }
#' @section Details:
#' \code{methylRawDB} is created via \code{\link{read}} function and has the same functionality as \code{\link{methylRaw}} class, 
#' but the data is saved in a flat database file and therefore allocates less space in memory.
#' 
#' @section Subsetting:
#'  In the following code snippets, \code{x} is a \code{methylRawDB}.
#'  Subsetting by \code{x[i,]} will produce a new \code{methylRaw} object if subsetting is done on
#'  rows. Column subsetting is not directly allowed to prevent errors in the 
#'  downstream analysis. see ?methylKit[ .
#' 
#' @section Accessors:
#' The following functions provides access to data slots of methylRawDB:
#' \code{\link[methylKit]{getData}},\code{\link[methylKit]{getAssembly}},
#' \code{\link[methylKit]{getContext}}
#' 
#' @section Coercion:
#'   \code{methylRawDB} object can be coerced to 
#'   \code{\link[GenomicRanges]{GRanges}} object via \code{\link{as}} function.
#' 
#' @examples
#' 
#' # example of a raw methylation data contained as a text file
#' read.table(system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#' header=TRUE,nrows=5)
#' 
#' methylRawListDB.obj <- read(file.list,sample.id = list("test1","test2","ctrl1","ctrl2"),
#'                             assembly = "hg18",treatment = c(1,1,0,0),
#'                             dbtype = "tabix",dbdir = "methylDB")
#' 
#' # example of a methylRawDB object
#' methylRawListDB.obj[[1]]
#' str(methylRawListDB.obj[[1]])
#' 
#' library(GenomicRanges)
#' 
#' #coercing methylRawDB object to GRanges object
#' my.gr=as(methylRawListDB.obj[[1]],"GRanges")
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
                          resolution){
  filepath=paste0(dbpath,"/",sample.id,".txt")
  df <- df[with(df,order(chr,start,end)),]
  df2tabix(df,filepath)
  num.records=Rsamtools::countTabix(paste0(filepath,".bgz"))[[1]] ## 
  
  new("methylRawDB",dbpath=paste0(filepath,".bgz"),num.records=num.records,
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
  
  new("methylRawDB",dbpath=normalizePath(dbpath),num.records=num.records,
      sample.id = sample.id, assembly = assembly,context=context,
      resolution=resolution,dbtype=dbtype)
}


# PRIVATE function:
# selects records in \code{object} of class \code{methylRawDB} 
# that lie inside the regions given by \code{ranges} of class \code{GRanges}. 
selectByOverlap<-function(object, ranges){
  
  return( getTabixByOverlap(tbxFile = object@dbpath,granges = ranges, return.type = "data.frame") )
  
}
  

# methylRawListDB

valid.methylRawListDB <- function(object) {
  
  
  # if all elements are methyl
  if(!all(sapply(object,class)=="methylRawDB")){
    FALSE
  }
  else if ( length(object) != length(object@treatment) ){
    message("The number of samples is different from the number of treatments, check the length of 'treatment'")
    FALSE
  }
  else{
    TRUE
  }
  
}


#' An S4 class for holding a list of methylRawDB objects.
#'
#' This class stores the list of  \code{\link[methylKit]{methylRawDB}} objects.
#' Functions such as \code{lapply} can be used on this list. It extends
#'  \code{\link[base]{list}} class. This object is primarily produced
#' by \code{\link[methylKit]{read}} function.
#'
#' @section Slots:\describe{
#'                  \item{\code{treatment}}{numeric vector denoting control and test samples}
#'                  \item{\code{.Data}}{a list of \code{\link{methylRawDB}} objects  } 
#'                }
#'                
#' @examples
#' file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
#'                 system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
#'                 system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#'                 system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
#' 
#' methylRawListDB.obj <- read(file.list,sample.id = list("test1","test2","ctrl1","ctrl2"),
#'                             assembly = "hg18",treatment = c(1,1,0,0),
#'                             dbtype = "tabix",dbdir = "methylDB")
#' 
#' #applying functions designed for methylRawDB on methylRawListDB object
#' lapply(methylRawListDB.obj,"getAssembly")
#'
#' @name methylRawListDB-class
#' @aliases methylRawListDB
#' @docType class
#' @rdname methylRawListDB-class
#' @export
setClass("methylRawListDB", slots=list(treatment = "vector"),contains = "list",
         validity=valid.methylRawListDB)


#' @aliases filterByCoverage,methylRawDB-method
#' @rdname filterByCoverage-methods
setMethod("filterByCoverage", signature(methylObj="methylRawDB"),
          function(methylObj,lo.count,lo.perc,hi.count,hi.perc){
            if( is.null(lo.count) & is.null(lo.perc) & is.null(hi.count) & is.null(hi.perc) ){return(methylObj)}
            
            filter <- function(data,lo.count,lo.perc,hi.count,hi.perc) {
              
              .setMethylDBNames(data)
            
              #figure out which cut-offs to use, maybe there is more elagent ways, quick&dirty works for now
              if(is.numeric(lo.count) ){lo.count=lo.count}
              if(is.numeric(lo.perc)){lo.count=quantile(data$coverage,lo.perc/100)}
              if(is.numeric(hi.count)){hi.count=hi.count}
              if(is.numeric(hi.perc)){hi.count=quantile(data$coverage,hi.perc/100)}
              
              if(is.numeric(lo.count)){data=data[data$coverage>=lo.count,]}
              if(is.numeric(hi.count)){data=data[data$coverage<hi.count,]}
              
              return(data)
              
            }
            
            dir <- dirname(methylObj@dbpath)
            filename <- paste(basename(tools::file_path_sans_ext(methylObj@dbpath)),"filtered",sep="_")
            
            newdbpath <- applyTbxByChunk(methylObj@dbpath, dir=dir,filename = filename, 
                                         return.type = "tabix", FUN = filter, lo.count=lo.count, lo.perc=lo.perc, 
                                         hi.count=hi.count, hi.perc=hi.perc)
            
            readMethylRawDB(dbpath = newdbpath,dbtype = "tabix",
                            sample.id = methylObj@sample.id,
                            assembly = methylObj@assembly, context = methylObj@context,
                            resolution = methylObj@resolution)
            
          })

#' @aliases filterByCoverage,methylRawListDB-method
#' @rdname filterByCoverage-methods
setMethod("filterByCoverage", signature(methylObj="methylRawListDB"),
          function(methylObj,lo.count,lo.perc,hi.count,hi.perc){

              new.list=lapply(methylObj,filterByCoverage,lo.count,lo.perc,hi.count,hi.perc)
              new("methylRawListDB", new.list,treatment=methylObj@treatment)
            
          })

#' @rdname getCoverageStats-methods
#' @aliases getCoverageStats,methylRawDB-method
setMethod("getCoverageStats", "methylRawDB",
          function(object,plot,both.strands,labels,...){
            
            tmp = applyTbxByChunk(object@dbpath,return.type = "data.table", 
                                  FUN = function(x) { .setMethylDBNames(x); return(x[,.(strand,coverage)])} )
            
            if(!plot){
              qts=seq(0,0.9,0.1) # get quantiles
              qts=c(qts,0.95,0.99,0.995,0.999,1)                          
              
              if(both.strands){
                
                
                plus.cov=tmp[strand=="+",coverage]
                mnus.cov=tmp[strand=="-",coverage]
                
                cat("read coverage statistics per base\n\n")
                cat("FORWARD STRAND:\n")
                cat("summary:\n")
                print( summary( plus.cov ) )
                cat("percentiles:\n")
                print(quantile( plus.cov,p=qts ))
                cat("\n\n")
                cat("REVERSE STRAND:\n")
                cat("summary:\n")
                print( summary( mnus.cov ) )
                cat("percentiles:\n")
                print(quantile( mnus.cov,p=qts ))
                cat("\n")                          
              }else{
                
                all.cov=tmp[,coverage]
                
                cat("read coverage statistics per base\n")
                cat("summary:\n")
                print( summary( all.cov ) )
                cat("percentiles:\n")
                print(quantile( all.cov,p=qts ))
                cat("\n")
              }
              
            }else{
              if(both.strands){   
                plus.cov=tmp[strand=="+",coverage]
                mnus.cov=tmp[strand=="-",coverage]
                
                par(mfrow=c(1,2))
                if(labels){
                  a=hist(log10(plus.cov),plot=F)
                  my.labs=as.character(round(100*a$counts/length(plus.cov),1))
                }else{my.labs=F}
                
                hist(log10(plus.cov),col="chartreuse4",
                     xlab=paste("log10 of read coverage per",object@resolution),
                     main=paste("Histogram of", object@context, "coverage: Forward strand"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
                if(labels){
                  a=hist(log10(mnus.cov),plot=F)
                  my.labs=as.character(round(100*a$counts/length(mnus.cov),1))
                }else{my.labs=F}
                a=hist(log10(mnus.cov),plot=F)
                hist(log10(mnus.cov),col="chartreuse4",
                     xlab=paste("log10 of read coverage per",object@resolution),
                     main=paste("Histogram of", object@context, "coverage: Reverse strand"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
              }else{
                all.cov=tmp[,coverage]
                if(labels){
                  a=hist(log10(all.cov),plot=F)
                  my.labs=as.character(round(100*a$counts/length(all.cov),1))
                }else{my.labs=F}                          
                
                hist(log10(all.cov),col="chartreuse4",
                     xlab=paste("log10 of read coverage per",object@resolution),
                     main=paste("Histogram of", object@context, "coverage"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
              }
              
              
              
            }
            
            
            
          })

#' @rdname getMethylationStats-methods
#' @aliases getMethylationStats,methylRawDB-method
setMethod("getMethylationStats", "methylRawDB",
          function(object,plot,both.strands,labels,...){
            
            tmp = applyTbxByChunk(object@dbpath,return.type = "data.table", 
                                  FUN = function(x) { .setMethylDBNames(x); return(x[,.(strand,coverage,numCs)])} )
            
            
            plus.met=100* tmp[strand=="+",numCs/coverage]
            mnus.met=100* tmp[strand=="-",numCs/coverage]
            all.met =100* tmp[,numCs/coverage]
            
            if(!plot){
              qts=seq(0,0.9,0.1) # get quantiles
              qts=c(qts,0.95,0.99,0.995,0.999,1)                          
              
              if(both.strands){       
                
                cat("methylation statistics per base\n\n")
                cat("FORWARD STRAND:\n")
                cat("summary:\n")
                print( summary( plus.met ) )
                cat("percentiles:\n")
                print(quantile( plus.met,p=qts ))
                cat("\n\n")
                cat("REVERSE STRAND:\n")
                cat("summary:\n")
                print( summary( mnus.met ) )
                cat("percentiles:\n")
                print(quantile( mnus.met,p=qts ))
                cat("\n")                          
              }else{
                
                
                cat("methylation statistics per base\n")
                cat("summary:\n")
                print( summary( all.met ) )
                cat("percentiles:\n")
                print(quantile( all.met,p=qts ))
                cat("\n")
              }
              
            }else{
              if(both.strands){   
                
                par(mfrow=c(1,2))
                if(labels){
                  a=hist((plus.met),plot=F,...)
                  my.labs=as.character(round(100*a$counts/length(plus.met),1))
                }else{my.labs=FALSE}
                hist((plus.met),col="cornflowerblue",
                     xlab=paste("% methylation per",object@resolution),
                     main=paste("Histogram of %", object@context,"methylation: Forward strand"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
                if(labels){                          
                  a=hist((mnus.met),plot=F,...)
                  my.labs=as.character(round(100*a$counts/length(mnus.met),1))
                }
                else{my.labs=FALSE}
                
                hist((mnus.met),col="cornflowerblue",
                     xlab=paste("% methylation per",object@resolution),
                     main=paste("Histogram of %", object@context,"methylation: Reverse strand"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
              }else{
                if(labels){                          
                  
                  a=hist((all.met),plot=F,...)
                  my.labs=as.character(round(100*a$counts/length(all.met),1))
                }else{my.labs=FALSE}
                hist((all.met),col="cornflowerblue",
                     xlab=paste("% methylation per",object@resolution),
                     main=paste("Histogram of %", object@context,"methylation"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
              }
              
              
            }
            
            
          })


#' @rdname adjust.methylC
#' @aliases adjust.methylC,methylRawDB,methylRawDB-method
setMethod("adjust.methylC", c("methylRawDB","methylRawDB"),function(mc,hmc){
  
  lst=new("methylRawListDB",list(mc,hmc),treatment=c(1,0))
  base=unite(lst)
  
  adjust <- function(data) {
    
    .setMethylDBNames(data)
    diff=(data$numCs1)-round(data$coverage1*(data$numCs2/data$coverage2))
    diff[diff<0]=0
    data$numCs1=diff
    data$numTs1=data$coverage1-data$numCs1
    return(data[1:7])
    
  }
  
  dir <- dirname(mc@dbpath) 
  filename <- paste(basename(tools::file_path_sans_ext(mc@dbpath)),"adjust",sep="_")
  
  newdbpath <- applyTbxByChunk(base@dbpath, dir=dir,filename = filename, 
                               return.type = "tabix", FUN = adjust)
  
  unlink(list.files(dirname(base@dbpath),pattern = basename(tools::file_path_sans_ext(base@dbpath)),full.names = TRUE))
  
  readMethylRawDB(dbpath = newdbpath,sample.id=mc@sample.id,
                  assembly=mc@assembly, context =mc@context,resolution=mc@resolution,
                  dbtype = mc@dbtype)
  
  
})


#' @rdname adjust.methylC
#' @aliases adjust.methylC,methylRawListDB,methylRawListDB-method
setMethod("adjust.methylC", c("methylRawListDB","methylRawListDB"),function(mc,hmc){
  
  # check lengths equal if not give error
  if(length(mc) != length(hmc)){stop("lengths of methylRawList objects should be same\n")}
  
  my.list=list()
  for(i in 1:length(mc)){
    my.list[[i]]=adjust.methylC(mc[[i]],hmc[[i]])
  }
  new("methylRawListDB",my.list,treatment=mc@treatment )
  
})


# methylBaseDB ------------------------------------------------------------


valid.methylBaseDB <- function(object) {
  
  
  check1=( (object@resolution == "base") | (object@resolution == "region") )
  check2=file.exists(object@dbpath)
  check3=( length(object@sample.ids) != length(object@treatment) )
  
  if(check1 & check2 ){
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
  else if(! check3){
    message("The number of samples is different from the number of treatments, check the length of 'treatment'")
    FALSE
  }
  
}




#' An S4 class for storing methylation events sampled in multiple experiments as flat file database
#'
#' This class is designed to contain methylation information such as coverage, number of methylated bases, etc...
#' The class creates an object that holds methylation information and genomic location as flat file database.
#' The object belonging to this class is produced by \code{\link{unite}} function.
#'          
#' @section Slots:\describe{
#'                  \item{\code{dbpath}:}{path to flat file database(s) }
#'                  
#'                  \item{\code{sample.ids}:}{character vector for ids of samples in the object}
#'
#'                  \item{\code{assembly}:}{name of the genome assembly}
#'
#'                  \item{\code{context}:}{context of methylation. Ex: CpG,CpH,CHH, etc}
#'
#'                  \item{\code{treatment}:}{treatment vector denoting which samples are test and control}
#'
#'                  \item{\code{coverage.index}:}{vector denoting which columns in the data correspons to coverage values}
#'
#'                  \item{\code{numCs.index}:}{vector denoting which columns in the data correspons to number of methylatedCs values}
#'                  \item{\code{numTs.index}:}{vector denoting which columns in the data correspons to number of unmethylated Cs values}
#'                  \item{\code{destranded}:}{ logical value. If \code{TRUE} object is destranded, if \code{FALSE} it is not.}
#'                  \item{\code{resolution}:}{ resolution of methylation information, allowed values: 'base' or 'region'}
#'                  \item{\code{dbtype}:}{string for type of the flat file database, ex: tabix}
#' }
#' 
#' @section Details:
#' \code{methylBaseDB} class has the same functionality as \code{\link{methylBase}} class, 
#' but the data is saved in a flat database file and therefore allocates less space in memory.
#' 
#' 
#' @section Subsetting:
#'  In the following code snippets, \code{x} is a \code{methylBaseDB}.
#'  Subsetting by \code{x[i,]} will produce a new \code{methylBase} object if subsetting is done on
#'  rows. Column subsetting is not directly allowed to prevent errors in the 
#'  downstream analysis. see ?methylKit[ .
#' 
#' @section Accessors:
#' The following functions provides access to data slots of \code{methylBaseDB}:
#' \code{\link[methylKit]{getData}},\code{\link[methylKit]{getAssembly}},
#' \code{\link[methylKit]{getContext}}
#' 
#' 
#' @section Coercion:
#'   \code{methylBaseDB} object can be coerced to \code{\link[GenomicRanges]{GRanges}} object via \code{\link{as}} function.
#' 
#' 
#' @examples
#' data(methylKit)
#' library(GenomicRanges)
#' my.gr=as(methylBaseDB.obj,"GRanges")
#' 
#' @name methylBaseDB-class
#' @aliases methylBaseDB
#' @docType class
#' @rdname methylBaseDB-class
#' @export
setClass("methylBaseDB",slots=list(dbpath = "character", num.records = "numeric", sample.ids = "character",
                                   assembly = "character", context = "character", resolution = "character",
                                   dbtype = "character", treatment = "numeric", coverage.index = "numeric",
                                   numCs.index = "numeric", numTs.index = "numeric", destranded = "logical"),
         validity = valid.methylBaseDB)

# PRIVATE function:
# makes a methylBaseDB object from a df
# it is called from read function or whenever this functionality is needed
makeMethylBaseDB<-function(df,dbpath,dbtype,
                           sample.ids, assembly ,context,
                           resolution,treatment,coverage.index,
                           numCs.index,numTs.index,destranded){
  
  # new tabix file is named by concatenation of sample.ids, works for now
  filepath=paste0(dbpath,"/",paste0(sample.ids,collapse = "_"),"_base.txt")
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
                           resolution,treatment,coverage.index,
                           numCs.index,numTs.index,destranded,skip=0){
  
  if(!file.exists(paste0(dbpath,".tbi"))) 
  {  
    Rsamtools::indexTabix(dbpath,seq=1, start=2, end=3,
                          skip=skip, comment="#", zeroBased=FALSE)
  }
  num.records=Rsamtools::countTabix(dbpath)[[1]] ## 
  
  new("methylBaseDB",dbpath=normalizePath(dbpath),num.records=num.records,
      sample.ids = sample.ids, assembly = assembly,context=context,
      resolution=resolution,dbtype=dbtype,treatment=treatment,
      coverage.index=coverage.index,numCs.index=numCs.index,numTs.index=numTs.index,
      destranded=destranded)
}


#' @rdname unite-methods
#' @aliases unite,methylRawListDB-method
setMethod("unite", "methylRawListDB",
          function(object,destrand,min.per.group){
            
            #check if assemblies,contexts and resolutions are same type NOT IMPLEMENTED   
            if( length(unique(vapply(object,function(x) x@context,FUN.VALUE="character"))) > 1)
            {
              stop("supplied methylRawList object have different methylation contexts:not all methylation events from the same bases")
            }
            if( length(unique(vapply(object,function(x) x@assembly,FUN.VALUE="character"))) > 1)
            {
              stop("supplied methylRawList object have different genome assemblies")
            }                     
            if( length(unique(vapply(object,function(x) x@resolution,FUN.VALUE="character"))) > 1)
            {
              stop("supplied methylRawList object have different methylation resolutions:some base-pair some regional")
            } 
            
            if( (!is.null(min.per.group)) &  ( ! is.integer( min.per.group ) )  ){stop("min.per.group should be an integer\ntry providing integers as 1L, 2L,3L etc.\n")}
            
            # destrand single objects contained in methylRawListDB
            if(destrand) { 
              
              destrandFun <- function(obj){
              
                if(obj@resolution == "base") {
                  dir <- dirname(obj@dbpath)
                  filename <- paste(basename(tools::file_path_sans_ext(obj@dbpath)),"destrand",sep="_")
                  # need to use .CpG.dinuc.unifyOld because output needs to be ordered
                  newdbpath <- applyTbxByChunk(obj@dbpath, dir=dir,filename = filename,return.type = "tabix", 
                                               FUN = function(x) { .CpG.dinuc.unifyOld(.setMethylDBNames(x) )})
                  
                  readMethylRawDB(dbpath = newdbpath,dbtype = "tabix",
                                  sample.id = obj@sample.id,
                                  assembly = obj@assembly, context = obj@context,
                                  resolution = obj@resolution)
                }else {obj}
                
              }
              new.list=lapply(object,destrandFun)
              object <- new("methylRawListDB", new.list,treatment=object@treatment)
              
              on.exit(unlink(list.files(dirname(dbpath),pattern = "destrand",full.names = TRUE)))    
            }
            #merge raw methylation calls together
            
            objList <- sapply(object,FUN = function(x) x@dbpath)
            dir <- dirname(object[[1]]@dbpath)
            filename <- paste(getSampleID(object),collapse = "_")
            
            if(is.null(min.per.group)) {
              
               dbpath <- mergeTabix(tabixList = objList ,dir = dir,filename = filename,mc.cores = 1) 
               
               dbpath <- tools::file_path_sans_ext(dbpath)
               
            } else {
              # if the the min.per.group argument is supplied, remove the rows that doesn't have enough coverage
              
              # keep rows with no matching in all samples  
              tmpPath <- mergeTabix(tabixList = objList ,dir = dir,filename = paste0(filename,".tmp"),mc.cores = 1,all=TRUE) 
              tmpPath <- tools::file_path_sans_ext(tmpPath)
              
              # get indices of coverage in the data frame 
              coverage.ind=seq(5,by=3,length.out=length(object))

              filter <- function(df, coverage.ind, treatment,min.per.group){
                
                df <- as.data.table(df)
                
                for(i in unique(treatment) ){
                  
                  my.ind=coverage.ind[treatment==i]
                  ldat = !is.na(df[,my.ind,with=FALSE])
                  
                  if(  is.null(dim(ldat))  ){  # if there is only one dimension
                    df=df[ldat>=min.per.group,]
                  }else{
                    df=df[rowSums(ldat)>=min.per.group,]
                  }
                  
                }
                return(as.data.frame(df))
              }

              dbpath <- applyTbxByChunk(tbxFile = tmpPath, dir = dir, filename = filename, 
                                        return.type = "tabix", FUN = filter,treatment=object@treatment,
                                        coverage.ind=coverage.ind,min.per.group=min.per.group)
              
              unlink(list.files(dirname(tmpPath),pattern = basename(tools::file_path_sans_ext(tmpPath)),full.names = TRUE))
              }
            
            
            coverage.ind=seq(5,by=3,length.out=length(object))
            numCs.ind   =coverage.ind+1
            numTs.ind   =coverage.ind+2
            
            readMethylBaseDB(dbpath = dbpath,dbtype = object[[1]]@dbtype,
                             sample.ids = getSampleID(object),assembly = object[[1]]@assembly,
                             context = object[[1]]@context,resolution = object[[1]]@resolution,
                             treatment = object@treatment,coverage.index = coverage.ind,
                             numCs.index = numCs.ind,numTs.index = numTs.ind,
                             destranded = destrand)
            
          }
)           



#' @rdname getCorrelation-methods
#' @aliases getCorrelation,methylBaseDB-method
setMethod("getCorrelation", "methylBaseDB",
          function(object,method,plot,nrow=2e6){
            
            if(is.null(nrow)){ 
              
              meth.fun <- function(data, numCs.index, numTs.index){
                
                data[, numCs.index]/( data[,numCs.index] + data[,numTs.index] )
                
              }
              meth.mat = applyTbxByChunk(object@dbpath,return.type = "data.frame",
                                         FUN = meth.fun, numCs.index = object@numCs.index,
                                         numTs.index = object@numTs.index)
              
            }else{
              data = headTabix(object@dbpath,nrow=nrow,return.type = "data.frame")
            
              meth.mat = data[, object@numCs.index]/( data[,object@numCs.index] + data[,object@numTs.index] )
              
              }
            
            names(meth.mat)=object@sample.ids


            print( cor(meth.mat,method=method) )
            
            
            panel.cor.pearson <- function(x, y, digits=2, prefix="", cex.cor, ...)
            {
              usr <- par("usr"); on.exit(par(usr))
              par(usr = c(0, 1, 0, 1))
              r <- abs(cor(x, y,method="pearson"))
              txt <- format(c(r, 0.123456789), digits=digits)[1]
              txt <- paste(prefix, txt, sep="")
              if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
              text(0.5, 0.5, txt, cex = cex.cor * r)
            }
            
            panel.cor.kendall <- function(x, y, digits=2, prefix="", cex.cor, ...)
            {
              usr <- par("usr"); on.exit(par(usr))
              par(usr = c(0, 1, 0, 1))
              r <- abs(cor(x, y,method="kendall"))
              txt <- format(c(r, 0.123456789), digits=digits)[1]
              txt <- paste(prefix, txt, sep="")
              if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
              text(0.5, 0.5, txt, cex = cex.cor * r)
            }
            
            panel.cor.spearman <- function(x, y, digits=2, prefix="", cex.cor, ...)
            {
              usr <- par("usr"); on.exit(par(usr))
              par(usr = c(0, 1, 0, 1))
              r <- abs(cor(x, y,method="spearman"))
              txt <- format(c(r, 0.123456789), digits=digits)[1]
              txt <- paste(prefix, txt, sep="")
              if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
              text(0.5, 0.5, txt, cex = cex.cor * r)
            }
            
            
            
            panel.my.smooth2<-function(x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = "darkgreen", span = 2/3, iter = 3, ...) 
            {
              par(new = TRUE)    #par(usr = c(usr[1:2], 0, 1.5) )
              smoothScatter(x, y,colramp=colorRampPalette(topo.colors(100)), bg = bg)
              ok <- is.finite(x) & is.finite(y)
              if (any(ok)) 
                lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = col.smooth, ...)
              abline(lm(y[ok]~x[ok]), col="red")
            }
            
            panel.my.smooth<-function(x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 0.3, col.smooth = "green", span = 2/3, iter = 3, ...) 
            {
              points(x, y, pch = 20, col = densCols(x,y,colramp=colorRampPalette(topo.colors(20))), bg = bg, cex = 0.1)
              ok <- is.finite(x) & is.finite(y)
              if (any(ok)){
                lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = col.smooth, ...);
                abline(lm(y[ok]~x[ok]), col="red")}
            }
            panel.hist <- function(x, ...)
            {
              usr <- par("usr"); on.exit(par(usr))
              par(usr = c(usr[1:2], 0, 1.5) )
              h <- hist(x, plot = FALSE)
              breaks <- h$breaks; nB <- length(breaks)
              y <- h$counts; y <- y/max(y)
              rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
            }
            
            if(plot)
            {  
              
              if(method=="spearman")
              { pairs(meth.mat, 
                      lower.panel=panel.my.smooth2, 
                      upper.panel=panel.cor.spearman,
                      diag.panel=panel.hist,main=paste(object@context, object@resolution ,method,"cor.") )
              }
              if(method=="kendall")
              { pairs(meth.mat, 
                      lower.panel=panel.my.smooth2, 
                      upper.panel=panel.cor.kendall,
                      diag.panel=panel.hist,main=paste(object@context, object@resolution ,method,"cor.") )
              }
              if(method=="pearson")
              { pairs(meth.mat, 
                      lower.panel=panel.my.smooth2, 
                      upper.panel=panel.cor.pearson,
                      diag.panel=panel.hist,main=paste(object@context, object@resolution ,method,"cor.") )
              }
              
              
            }
          }  
)




#' @rdname percMethylation-methods
#' @aliases percMethylation,methylBaseDB-method
setMethod("percMethylation", "methylBaseDB",
          function(methylBase.obj,rowids=FALSE,save.txt){
            
            meth.fun <- function(data, numCs.index, numTs.index){
              
              100 * data[, numCs.index]/( data[,numCs.index] + data[,numTs.index] )
              
            }
            if (save.txt) {
              
              filename <- paste(basename(tools::file_path_sans_ext(methylBase.obj@dbpath)),"methMath.txt",sep = "_")
              
              meth.mat = applyTbxByChunk(methylBase.obj@dbpath,return.type = "text",
                                         dir = dirname(methylBase.obj@dbpath), filename = filename,
                                         FUN = meth.fun, numCs.index = methylBase.obj@numCs.index,
                                         numTs.index = methylBase.obj@numTs.index)
              
              return(meth.mat)
              
            } else {
              
              meth.mat = applyTbxByChunk(methylBase.obj@dbpath,return.type = "data.frame",
                                         FUN = meth.fun, numCs.index = methylBase.obj@numCs.index,
                                         numTs.index = methylBase.obj@numTs.index)
              
              names(meth.mat)=methylBase.obj@sample.ids
              if(rowids){
                rownames(meth.mat)=as.character(paste(x[,1],x[,2],x[,3],sep=".") )
              }
              return(as.matrix(meth.mat))
              
            }
            })



#' @rdname reconstruct-methods
#' @aliases reconstruct,methylBaseDB-method
setMethod("reconstruct",signature(mBase="methylBaseDB"), function(methMat,mBase){
  
  if(is.matrix(methMat) ) {
      
    # check if indeed methMat is percent methylation matrix
    if(max(methMat)<=1){
      warning("\nmake sure 'methMat' is percent methylation matrix (values between 0-100) \n")
    }
    
    # check if methMat is percent methylation matrix fitting to mBase  
    if(nrow(methMat) != mBase@num.records | ncol(methMat) != length(mBase@numCs.index) ){
      stop("\nmethMat dimensions do not match number of samples\n",
           "and number of bases in methylBase object\n")
    }
    
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    matFile=paste(rndFile,"methMat.txt",sep="_")
    write.table(methMat,matFile,quote=FALSE,col.names=FALSE,row.names=FALSE,
                sep="\t")
    methMat = matFile
    
  } else {
    
    mat <- read.table(methMat,header = F)
    
    # check if indeed methMat is percent methylation matrix
    if(max(mat)<=1){
      warning("\nmake sure 'methMat' is percent methylation matrix (values between 0-100) \n")
    }
    
    # check if methMat is percent methylation matrix fitting to mBase  
    if(nrow(mat) != mBase@num.records | ncol(mat) != length(mBase@numCs.index) ){
      stop("\nmethMat dimensions do not match number of samples\n",
           "and number of bases in methylBase object\n")
    }
    
  }
  
  
  reconstr <- function(data, methMat, chunk, numCs.index, numTs.index) {
    
    mat=data[,numCs.index]+data[,numTs.index]
    
    methMat = read.table(methMat,header = F, nrows = chunk)
    
    # get new unmethylated and methylated counts
    numCs=round(methMat*mat/100)
    numTs=round((100-methMat)*mat/100)
    
    data[,numCs.index]=numCs
    data[,numTs.index]=numTs
    
    return(data)
  }
  
  dir <- dirname(mBase@dbpath)
  filename <- paste(basename(tools::file_path_sans_ext(mBase@dbpath)),"reconstruct",sep="_")
  chunk.size <- 1e6
  con <- file(methMat,open = "r") 
  
  newdbpath <- applyTbxByChunk(tbxFile = mBase@dbpath,chunk.size = chunk.size, dir=dir,filename = filename, 
                               return.type = "tabix", FUN = reconstr, con,chunk.size,mBase@numCs.index,mBase@numTs.index)
  
  close(con)
  if(file.exists(matFile)) {unlink(matFile)}
  
  readMethylBaseDB(dbpath = newdbpath,dbtype = mBase@dbtype,
                   sample.ids = mBase@sample.ids,assembly = mBase@assembly,
                   context = mBase@context,resolution = mBase@resolution,
                   treatment = mBase@treatment,coverage.index = mBase@coverage.index,
                   numCs.index = mBase@numCs.index,numTs.index = mBase@numTs.index,
                   destranded = mBase@destranded)
  
  
  
}
)



# methylDiffDB -------------------------------------------------------


valid.methylDiffDB <- function(object) {
  
  
  check1=( (object@resolution == "base") | (object@resolution == "region") )
  check2=file.exists(object@dbpath)
  check3=( length(object@sample.ids) != length(object@treatment) )
  
  if(check1 & check2 ){
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
  else if(! check3){
    message("The number of samples is different from the number of treatments, check the length of 'treatment'")
    FALSE
  }
  
}

#' An S4 class that holds differential methylation information as flat file database
#'
#' This class is designed to hold statistics and locations for differentially 
#' methylated regions/bases as flat file database.
#' \code{\link[methylKit]{calculateDiffMeth}} function returns an object 
#' with \code{methylDiffDB} class.
#'          
#' @section Slots:\describe{
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
#'
#' }
#' 
#' @section Details:
#' \code{methylDiffDB} class has the same functionality as \code{\link{methylDiff}} class, 
#' but the data is saved in a flat database file and therefore allocates less space in memory.
#' 
#' 
#' @section Subsetting:
#'  In the following code snippets, \code{x} is a \code{methylDiffDB}.
#'  Subsetting by \code{x[i,]} will produce a new object if subsetting is done 
#'  on rows. Column subsetting is not directly allowed to prevent errors in the 
#'  downstream analysis. see ?methylKit[ .
#' 
#' @section Coercion:
#'   \code{methylDiffDB} object can be coerced to \code{\link[GenomicRanges]{GRanges}} object via \code{\link{as}} function.
#' 
#' @section Accessors: 
#' The following functions provides access to data slots of methylDiffDB:
#' \code{\link[methylKit]{getData}},\code{\link[methylKit]{getAssembly}},\code{\link[methylKit]{getContext}}
#' 
#' @examples
#' data(methylKit)
#' library(GenomicRanges)
#' my.gr=as(methylDiffDB.obj,"GRanges")
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
                           resolution,treatment,destranded){
  
  # new tabix file is named by concatenation of sample.ids, works for now
  filepath=paste0(dbpath,"/",paste0(sample.ids,collapse = "_"),"_diff.txt")
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
  
  new("methylDiffDB",dbpath=normalizePath(dbpath),num.records=num.records,
      sample.ids = sample.ids, assembly = assembly,context=context,
      resolution=resolution,dbtype=dbtype,treatment=treatment,
      destranded=destranded)
}


setMethod("calculateDiffMeth", "methylBaseDB",
          function(.Object,covariates,overdispersion=c("none","MN","shrinkMN"),
                   adjust=c("SLIM","holm","hochberg","hommel","bonferroni","BH","BY","fdr","none","qvalue"),
                   effect=c("wmean","mean","predicted"),parShrinkNM=list(),
                   test=c("F","Chisq"),mc.cores=1,slim=TRUE,weighted.mean=TRUE){
            
            if(length(.Object@treatment)<2 ){
              stop("can not do differential methylation calculation with less than two samples")
            }
            
            if(length(unique(.Object@treatment))<2 ){
              stop("can not do differential methylation calculation when there is no control\n
                   treatment option should have 0 and 1 designating treatment and control samples")
            }
            
            if(length(unique(.Object@treatment))>2 ){
              stop("can not do differential methylation calculation when there are more than\n
                   two groups, treatment vector indicates more than two groups")
            }
            
            #### check if covariates+intercept+treatment more than replicates ####
            if(!is.null(covariates)){if(ncol(covariates)+2 >= length(.Object@numTs.index)){stop("Too many covariates/too few replicates.")}}
            
            # add backwards compatibility with old parameters
            if(slim==FALSE) adjust="BH" else adjust=adjust
            if(weighted.mean==FALSE) effect="mean" else effect=effect

                        
            # function to apply the test to data
            diffMeth <- function(data,Ccols,Tcols,formula,vars,treatment,overdispersion,effect,
                                 parShrinkNM,test,mc.cores){
              
              cntlist=split(as.matrix(data[,c(Ccols,Tcols)]),1:nrow(data))
              
              tmp=simplify2array(
                mclapply(cntlist,logReg,formula,vars,treatment=.Object@treatment,overdispersion=overdispersion,effect=effect,
                         parShrinkNM=parShrinkNM,test=test,mc.cores=mc.cores))
              tmp <- as.data.frame(t(tmp))
              #print(head(tmp))
              x=data.frame(data[,1:4],tmp$p.value,p.adjusted(tmp$q.value,method=adjust),meth.diff=tmp$meth.diff.1,stringsAsFactors=FALSE)
              
              return(x)
              
              
            }
              
            dir <- dirname(.Object@dbpath)
            filename <- paste(basename(tools::file_path_sans_ext(.Object@dbpath)),"diffMeth",sep="_")
            chunk.size <- 1e4
            
            dbpath <- applyTbxByChunk(.Object@dbpath,dir = dir,chunk.size = chunk.size,  filename = filename, return.type = "tabix", FUN = diffMeth,
                                   Ccols = .Object@numCs.index,Tcols = .Object@numTs.index,formula=formula,vars=covariates,
                                   treatment=.Object@treatment,overdispersion=overdispersion,effect=effect,
                                   parShrinkNM=parShrinkNM,test=test,mc.cores=mc.cores)
            
            
            obj=readMethylDiffDB(dbpath = dbpath,dbtype = .Object@dbtype, sample.ids=.Object@sample.ids,assembly=.Object@assembly,context=.Object@context,
                    destranded=.Object@destranded,treatment=.Object@treatment,resolution=.Object@resolution)
            obj
            }
)

#' @aliases get.methylDiff,methylDiffDB-method
#' @rdname get.methylDiff-methods
setMethod(f="get.methylDiff", signature="methylDiffDB", 
          definition=function(.Object,difference,qvalue,type) {
            
            if(!( type %in% c("all","hyper","hypo") )){
              stop("Wrong 'type' argument supplied for the function, it can be 'hypo', 'hyper' or 'all' ")
            }
              
            #function applied to data
            f <- function(data,difference,qv,type){
              
              data <- data.table(data)
              .setMethylDBNames(data,methylDBclass = "methylDiffDB")
              
              if(type=="all"){
                data <- data[(qvalue < qv) & (abs(meth.diff) > difference)]
              }else if(type=="hyper"){
                data <- data[(qvalue < qv) & (meth.diff > difference)]
              }else if(type=="hypo"){
                data <- data[(qvalue < qv) & (meth.diff < -1*difference)]
              }
              return(data)
            }
            
            dir <- dirname(.Object@dbpath)
            filename <- paste(basename(tools::file_path_sans_ext(.Object@dbpath)),type,sep="_")
            
            dbpath <- applyTbxByChunk(.Object@dbpath,dir = dir, filename = filename, return.type = "tabix", FUN = f,
                                      difference = difference, qv = qvalue, type = type)
            
            
            obj=readMethylDiffDB(dbpath = dbpath,dbtype = .Object@dbtype, sample.ids=.Object@sample.ids,assembly=.Object@assembly,context=.Object@context,
                                 destranded=.Object@destranded,treatment=.Object@treatment,resolution=.Object@resolution)
            return(obj)

            }) 



# bedgraph methods --------------------------------------------------------

#' @rdname bedgraph-methods
#' @aliases bedgraph,methylRawDB-method
setMethod("bedgraph", signature(methylObj="methylRawDB"),
          function(methylObj,file.name,col.name,unmeth,log.transform,negative,add.on){
            if(!col.name %in%  c('coverage', 'numCs','numTs','perc.meth') ){
              stop("col.name argument is not one of 'coverage', 'numCs','numTs','perc.meth'")
            }
           
            bedgr <- function(data,col.name,file.name,unmeth,log.transform,negative,add.on,sample.id){
              
              data <- as.data.table(data)
              .setMethylDBNames(df = data,methylDBclass = "methylRawDB")

              if(col.name=="perc.meth"){
                df= data[,.(chr,start=start-1,end,score=100*numCs/coverage )]
              }else{
                df=data[,.(chr,start=start-1,end,score=get(col.name) )]
                if(log.transform){
                  df=df[,.(chr,start,end,score=log10(score))]
                }
                if(negative){
                  df=df[,.(chr,start,end,score=-1*score )]
                }
              }
              
              if(is.null(file.name)){
                return(as.data.frame(df))
              }else{
              
                if(unmeth & col.name=="perc.meth")
                {
                  # write meth data to single file
                  write.table(df,file=paste0(file.name,"_meth"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t",append=TRUE)
                  
                  # write unmeth data to single file
                  dfu=df[,.(chr,start,end,score=100-score)]
                  write.table(dfu,file=paste0(file.name,"_unmeth"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t",append=TRUE)
                  
                }else{
                  write.table(df,file=file.name,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t",append=TRUE)
                }
              }
            }
            
            
            if(is.null(file.name)){
              
              applyTbxByChunk(methylObj@dbpath,chunk.size = 1e6, return.type = "data.frame", FUN = bedgr,
                              col.name = col.name, file.name = file.name, unmeth = unmeth, log.transform=log.transform, negative= negative,
                              add.on = add.on, sample.id = methylObj@sample.id)
              
            } else {
              
              rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
              filename2=paste(file.name,rndFile,sep="_")
              
              applyTbxByChunk(methylObj@dbpath,chunk.size = 1e6, return.type = "data.frame", FUN = bedgr,
                                    col.name = col.name, file.name= filename2, unmeth = unmeth, log.transform=log.transform, negative= negative,
                                    add.on = add.on, sample.id = methylObj@sample.id)
              
              
              if(unmeth & col.name=="perc.meth")
              {
                # combine single files produced by bedgr function
                track.line=paste(
                  "track type=bedGraph name='",methylObj@sample.id," METH Cs","' description='",methylObj@sample.id," METH Cs",
                  "' visibility=full color=255,0,0 maxHeightPixels=80:80:11 ",add.on,sep="")                        
                cat(track.line,"\n",file=file.name)
                file.append(file.name,paste0(filename2,"_meth"))
                track.line2=paste(
                  "track type=bedGraph name='",methylObj@sample.id," UNMETH Cs","' description='",methylObj@sample.id," UNMETH Cs",
                  "' visibility=full color=0,0,255 maxHeightPixels=80:80:11 ",add.on,sep="")
                cat(track.line2,"\n",file=file.name,append=TRUE)
                file.append(file.name,paste0(filename2,"_unmeth"))
              
              }else{
                
                track.line=paste(
                  "track type=bedGraph name='",methylObj@sample.id," ",col.name,"' description='",methylObj@sample.id," ",col.name,
                  "' visibility=full color=255,0,0 maxHeightPixels=80:80:11 ",add.on,sep="")
                cat(track.line,"\n",file=file.name)
                file.append(file.name,filename2)
                
              }
              # tidy up
              unlink(list.files(path = dirname(file.name),pattern = rndFile,full.names = T))
            
            }
          
          }
          
)


#' @rdname bedgraph-methods
#' @aliases bedgraph,methylRawListDB-method
setMethod("bedgraph", signature(methylObj="methylRawListDB"),
          function(methylObj,file.name,col.name,unmeth,log.transform,negative,add.on){
            if(!col.name %in%  c('coverage', 'numCs','numTs','perc.meth') ){
              stop("col.name argument is not one of 'coverage', 'numCs','numTs','perc.meth' options")
            }
            
            
            
            if( is.null(file.name) ){
              
              result=list()
              for(i in 1:length(methylObj))
              {
                result[[ methylObj[[i]]@sample.id ]]=bedgraph(methylObj[[i]],file.name=NULL,
                                                              col.name=col.name,unmeth=FALSE,log.transform=log.transform,
                                                              negative=negative)
              }
              
              return(result)
            }else{

                bedgraph(methylObj[[1]],file.name=file.name,
                         col.name=col.name,unmeth=unmeth,log.transform=log.transform,
                         negative=negative)
                
                for(i in 2:length(methylObj))
                {
                  bedgraph(methylObj[[i]],file.name=paste(file.name,i,sep = "_"),
                           col.name=col.name,unmeth=unmeth,log.transform=log.transform,
                           negative=negative)
                  
                  file.append(file.name,paste(file.name,i,sep = "_"))
                }
                
                unlink(list.files(path = dirname(file.name),pattern = paste0(basename(file.name),"_[[:digit:]]+"),full.names = T))
                
            }
          })


#' @rdname bedgraph-methods
#' @aliases bedgraph,methylDiffDB-method
setMethod("bedgraph", signature(methylObj="methylDiffDB"),
          function(methylObj,file.name,col.name,log.transform,negative,add.on){
            
            if(! col.name %in% c('pvalue','qvalue', 'meth.diff') )
            {
              stop("col.name argument is not one of 'pvalue','qvalue', 'meth.diff'")
            }
            
            bedgr <- function(data,col.name,file.name,log.transform,negative,add.on,sample.id){
              
              data <- as.data.table(data)
              .setMethylDBNames(df = data,methylDBclass = "methylDiffDB")
              
              df=data[,.(chr,start=start-1,end,score=get(col.name) )]
              if(log.transform){
                df=df[,.(chr,start,end,score=log10(score))]
              }
              if(negative){
                df=df[,.(chr,start,end,score=-1*score )]
              }
              if(is.null(file.name)){
                return(as.data.frame(df))
              }else{
                  write.table(df,file=file.name,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t",append=TRUE)
                }
             }
            
            if(is.null(file.name)){
              
              applyTbxByChunk(methylObj@dbpath,chunk.size = 1e6, return.type = "data.frame", FUN = bedgr,
                              col.name = col.name, file.name= file.name, log.transform=log.transform, negative= negative,
                              add.on = add.on, sample.id = methylObj@sample.id)
            } else {
              
              rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
              filename2=paste(file.name,rndFile,sep="_")
              
              txtPath <- applyTbxByChunk(methylObj@dbpath,chunk.size = 1e6, return.type = "data.frame", FUN = bedgr,
                                         col.name = col.name,file.name= filename2, log.transform=log.transform, negative= negative,
                                         add.on = add.on, sample.id = methylObj@sample.id)
              
              track.line=paste(
                "track type=bedGraph name='",file.name,"' description='",col.name,
                "' visibility=full color=255,0,255 altColor=102,205,170 maxHeightPixels=80:80:11 ",add.on,sep="")
              cat(track.line,"\n",file=file.name)
              file.append(file.name,filename2)
              
              unlink(filename2)
              
            }
        }
)

# coercion functions ------------------------------------------------------

setAs("methylRawDB", "GRanges", function(from)
{
  from2=getData(from)
  GRanges(seqnames=as.character(from2$chr),ranges=IRanges(start=from2$start, end=from2$end),
          strand=from2$strand, 
          coverage=from2$coverage,
          numCs   =from2$numCs,
          numTs  =from2$numTs                                
  )
  
})

setAs("methylBaseDB", "GRanges", function(from)
{
  from=getData(from)
  GRanges(seqnames=as.character(from$chr),ranges=IRanges(start=from$start, end=from$end),
          strand=from$strand, 
          data.frame(from[,5:ncol(from)])
  )
  
})

setAs("methylDiffDB", "GRanges", function(from)
{
  
  from = getData(from)
  GRanges(seqnames=as.character(from$chr),ranges=IRanges(start=from$start, end=from$end),
          strand=from$strand, 
          qvalue=from$qvalue,
          meth.diff=from$meth.diff
  )
  
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
  df <- headTabix(tbxFile = x@dbpath, nrow = x@num.records, return.type = "data.frame")
  .setMethylDBNames(df,"methylRawDB")
  
  return(df)
})

#' @rdname getData-methods
#' @aliases getData,methylBaseDB-method
setMethod("getData", signature="methylBaseDB", definition=function(x) {
  df <- headTabix(tbxFile = x@dbpath, nrow = x@num.records, return.type = "data.frame")
 .setMethylDBNames(df,"methylBaseDB")
  
  return(df)
})

#' @rdname getData-methods
#' @aliases getData,methylDiffDB-method
setMethod(f="getData", signature="methylDiffDB", definition=function(x) {
  df <- headTabix(tbxFile = x@dbpath, nrow = x@num.records, return.type = "data.frame")
  .setMethylDBNames(df,"methylDiffDB")
  
  return(df)
})


# show functions ----------------------------------------------------------

#' @rdname show-methods
#' @aliases show,methylRawDB
setMethod("show", "methylRawDB", function(object) {
  
  cat("methylRawDB object with",object@num.records,"rows\n--------------\n")
  print(.setMethylDBNames(headTabix(object@dbpath,nrow = 6,return.type = "data.frame"),methylDBclass = "methylRawDB"))
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
  print(.setMethylDBNames(headTabix(object@dbpath,nrow = 6,return.type = "data.frame"),methylDBclass = "methylBaseDB"))
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
  
  cat("methylDiffDB object with",nrow(object@num.records),"rows\n--------------\n")
  print(.setMethylDBNames(headTabix(object@dbpath,nrow = 6,return.type = "data.frame"),methylDBclass = "methylDiffDB"))
  cat("--------------\n")
  cat("sample.ids:",object@sample.ids,"\n")
  cat("destranded",object@destranded,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("treament:", object@treatment,"\n")
  cat("resolution:", object@resolution,"\n")
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


#' @aliases select,methylBaseDB-method
#' @rdname select-methods
setMethod("select", "methylBaseDB",
          function(x, i)
          {
            
            new("methylBase",getData(x)[i,],
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
            
            new("methylDiff",getData(x)[i,],
                sample.ids = x@sample.ids,
                assembly = x@assembly,
                context = x@context,
                treatment=x@treatment,
                destranded=x@destranded,
                resolution=x@resolution)
          }
)


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

#' @aliases [,methylBaseDB-method
#' @rdname extract-methods
setMethod("[",signature(x="methylBaseDB", i = "ANY", j="ANY"), 
          function(x,i,j,drop){
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
#' @aliases [,methylDiffDB-method
setMethod("[","methylDiffDB", 
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


#' get treatment vector of the methylRawListDB or methylBaseDB object
#' 
#' The function returns the treatment vector stored in any of the \code{\link{methylBaseDB}},\code{\link{methylRawListDB}} objects 
#' 
#' @param x an \code{\link{methylBaseDB}} or \code{\link{methylRawListDB}} object
#' @usage getTreatment(x)
#' @examples
#' 
#' data(methylKit)
#' 
#' getTreatment(methylBaseDB.obj)
#' getTreatment(methylRawListDB.obj)
#' 
#' 
#' @return the treatment vector for the object
#' @export
#' @docType methods
#' @rdname getTreatment-methods
setGeneric("getTreatment", def=function(x) standardGeneric("getTreatment"))
#' @rdname getTreatment-methods
#' @aliases getTreatment,methylRawListDB-method
setMethod("getTreatment", signature = "methylRawListDB", function(x) {
    return(x@treatment)
})

#' @rdname getTreatment-methods
#' @aliases getTreatment,methylBaseDB-method
setMethod("getTreatment", signature = "methylBaseDB", function(x) {
  return(x@treatment)
})

#' @rdname getTreatment-methods
#' @aliases getTreatment,methylDiffDB-method
setMethod("getTreatment", signature = "methylDiffDB", function(x) {
  return(x@treatment)
})


#' set treatment vector of the methylRawListDB or methylBaseDB object
#' 
#' The function directly modifies the treatment vector of the given \code{\link{methylBaseDB}},\code{\link{methylRawListDB}} object.
#' If you rather want to create a new object use 
#' 
#' @param x an \code{\link{methylBaseDB}} or \code{\link{methylRawListDB}} object
#' @usage setTreatment(x) <- c(1,1,0,0)
#' @examples
#' 
#' data(methylKit)
#' 
#' setTreatment(methylBaseDB.obj) <- c(2,2,1,1)
#' setTreatment(methylRawListDB.obj) <- c(2,2,1,1)
#' 
#' 
#' @return \code{\link{methylBaseDB}} or \code{\link{methylRawListDB}} object
#' @export
#' @docType methods
#' @rdname setTreatment-methods
setGeneric("setTreatment<-", def=function(x, value="numeric") {standardGeneric("setTreatment<-")})
#' @rdname setTreatment-methods
#' @aliases setTreatment,methylRawListDB-method
setReplaceMethod("setTreatment", signature = "methylRawListDB", function(x, value) {
  
  if(! ( length(x@treatment) == length(value) ) ){
    message("The new treatment vector is not valid, check the length of input")
    return(x)
  } else {
    x@treatment <- value
    return(x)
  }
  
})

#' @rdname setTreatment-methods
#' @aliases setTreatment,methylBaseDB-method
setReplaceMethod("setTreatment", signature = "methylBaseDB", function(x, value) {
  
  if(! ( length(x@treatment) == length(value) ) ){
    message("The new treatment vector is not valid, check the length of input")
    return(x)
  } else {
    x@treatment <- value
    return(x)
  }
  
})



#' get sample ids of the methylRawListDB or methylBaseDB object
#' 
#' The function returns the sample ids stored in any of the \code{\link{methylBaseDB}},\code{\link{methylRawListDB}} objects 
#' 
#' @param x an \code{\link{methylBaseDB}} or \code{\link{methylRawListDB}} object
#' @usage getSampleID(x)
#' @examples
#' 
#' data(methylKit)
#' 
#' getSampleID(methylBaseDB.obj)
#' getSampleID(methylRawListDB.obj)
#' 
#' 
#' @return the sample ids stored in the object
#' @export
#' @docType methods
#' @rdname getSampleID-methods
setGeneric("getSampleID", def=function(x) standardGeneric("getSampleID"))
#' @rdname getSampleID-methods
#' @aliases getSampleID,methylRawListDB-method
setMethod("getSampleID", signature = "methylRawListDB", function(x) {
  names <- vapply(x,function(z) z@sample.id,FUN.VALUE = "character")
  return(names)
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylBaseDB-method
setMethod("getSampleID", signature = "methylBaseDB", function(x) {
  return(x@sample.ids)
})


#' set sample ids of the methylRawListDB or methylBaseDB object
#' 
#' The function returns the \code{\link{methylBaseDB}},\code{\link{methylRawListDB}} objects with the new sample ids stored in it.
#' 
#' @param x an \code{\link{methylBaseDB}} or \code{\link{methylRawListDB}} object
#' @usage setSampleID(x) <- c("id.1","id.2","id.3",...,"id.n")
#' @examples
#' 
#' data(methylKit)
#' 
#' setSampleID(methylBaseDB.obj) <- c("t1","t2","c1","c2")
#' setSampleID(methylRawListDB.obj) <- c("test-one","test-two","control-one","control-two")
#' 
#' 
#' @return the object with new sample ids
#' @export
#' @docType methods
#' @rdname setSampleID-methods
setGeneric("setSampleID<-", def=function(x, value="numeric") {standardGeneric("setSampleID<-")})
#' @rdname setSampleID-methods
#' @aliases setSampleID,methylRawDB-method
setReplaceMethod("setSampleID", signature = "methylRawDB", function(x, value) {
  
  if(! ( length(value) == 1 ) ){
    message("The vector of new sample ids is not valid, check the length of input")
    return(x)
  } else {
    x@sample.id <- value
    return(x)
  }
  
})

#' @rdname setSampleID-methods
#' @aliases setSampleID,methylRawListDB-method
setReplaceMethod("setSampleID", signature = "methylRawListDB", function(x, value) {
  
  if(! ( length(getSampleID(x)) == length(value) ) ){
    message("The vector of new sample ids is not valid, check the length of input")
    return(x)
  } else {
    treatment <- x@treatment
    x <- mapply(`setSampleID<-`, x, value)
    x <- new("methylRawListDB",x,treatment=treatment)
    return(x)
  }
  
})

#' @rdname setSampleID-methods
#' @aliases setSampleID,methylBaseDB-method
setReplaceMethod("setSampleID", signature = "methylBaseDB", function(x, value) {
  
  if(! ( length(x@sample.ids) == length(value) ) ){
    message("The vector of new sample ids is not valid, check the length of input")
    return(x)
  } else {
    x@sample.ids <- value
    return(x)
  }
  
})


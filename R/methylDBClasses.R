
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
          function(methylObj,lo.count,lo.perc,hi.count,hi.perc, save.db){
            if( is.null(lo.count) & is.null(lo.perc) & is.null(hi.count) & is.null(hi.perc) ){return(methylObj)}
            
            data=getData(methylObj) # get the data part
            
            #figure out which cut-offs to use, maybe there is more elagent ways, quick&dirty works for now
            if(is.numeric(lo.count) ){lo.count=lo.count}
            if(is.numeric(lo.perc)){lo.count=quantile(data$coverage,lo.perc/100)}
            if(is.numeric(hi.count)){hi.count=hi.count}
            if(is.numeric(hi.perc)){hi.count=quantile(data$coverage,hi.perc/100)}
            
            if(is.numeric(lo.count)){data=data[data$coverage>=lo.count,]}
            if(is.numeric(hi.count)){data=data[data$coverage<hi.count,]}
            
            if(!save.db){
              
              new("methylRaw",data,sample.id=methylObj@sample.id,
                  assembly=methylObj@assembly,
                  context=methylObj@context,resolution=methylObj@resolution)
              
            }else{
              #sample name to not overwrite the original tabix file
              sample.id = paste(methylObj@sample.id,
                                paste(sample(c(0:9, letters, LETTERS),3,replace = T),collapse = ""),sep = "_")
              #save dataframe as methylbaseDB object and return object
              makeMethylRawDB(df = data,dbpath = getwd(),
                              dbtype = "tabix",sample.id = sample.id,
                              assembly = methylObj@assembly,context = methylObj@context,
                              resolution = methylObj@resolution)
              
            }
            
            
            
          })

#' @aliases filterByCoverage,methylRawListDB-method
#' @rdname filterByCoverage-methods
setMethod("filterByCoverage", signature(methylObj="methylRawListDB"),
          function(methylObj,lo.count,lo.perc,hi.count,hi.perc,save.db){
            if(!save.db){
              new.list=lapply(methylObj,filterByCoverage,lo.count,lo.perc,hi.count,hi.perc)
              new("methylRawList", new.list,treatment=methylObj@treatment)
            }else{
              new.list=lapply(methylObj,filterByCoverage,lo.count,lo.perc,hi.count,hi.perc,save.db)
              new("methylRawListDB", new.list,treatment=methylObj@treatment)
            }
            
          })



# methylBaseDB ------------------------------------------------------------


valid.methylBaseDB <- function(object) {
  
  
  #data=getData(object,nrow=5)
  check1=( (object@resolution == "base") | (object@resolution == "region") )
  check2=file.exists(object@dbpath)
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
  
}




#' An S4 class for methylation events sampled in multiple experiments
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
#' my.gr=as(methylBase.obj,"GRanges")
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
# it is called from modRead function or whenever this functionality is needed
makeMethylBaseDB<-function(df,dbpath,dbtype,
                           sample.ids, assembly ,context,
                           resolution,treatment,coverage.index,
                           numCs.index,numTs.index,destranded){
  
  # new tabix file is named by concatenation of sample.ids, works for now
  filepath=paste0(dbpath,"/",paste0(sample.ids,collapse = "_"),".txt")
  df <- df[with(df,order(chr,start,end)),]
  df2tabix(df,filepath)
  num.records=Rsamtools::countTabix(paste0(filepath,".bgz"))[[1]] ## 
  
  new("methylBaseDB",dbpath=paste0(filepath,".bgz"),num.records=num.records,
      sample.ids = sample.ids, assembly = assembly,context=context,
      resolution=resolution,dbtype=dbtype,treatment=treatment,
      coverage.index=coverage.index,numCs.index=numCs.index,numTs.index=numTs.index,
      destranded=destranded)
}


#' @rdname unite-methods
#' @aliases unite,methylRawListDB-method
setMethod("unite", "methylRawListDB",
          function(object,destrand,min.per.group,save.db=FALSE,dbdir=getwd()){
            
            
            
            
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
            
            #merge raw methylation calls together
            df=getData(object[[1]])
            if(destrand & (object[[1]]@resolution == "base") ){df=.CpG.dinuc.unify(df)}
            df=data.table(df,key=c("chr","start","end","strand"))
            sample.ids=c(object[[1]]@sample.id)
            assemblies=c(object[[1]]@assembly)
            contexts  =c(object[[1]]@context)
            for(i in 2:length(object))
            {
              df2=getData(object[[i]])
              if(destrand & (object[[1]]@resolution == "base") ){df2=.CpG.dinuc.unify(df2)}
              #
              
              if( is.null(min.per.group) ){
                df2=data.table(df2,key=c("chr","start","end","strand"))
                df=merge(df,df2,by=c("chr","start","end","strand"),suffixes=c(as.character(i-1),as.character(i) ) ) # merge the dat to a data.frame
                #df=df[df2, nomatch=FALSE]
              }else{
                df2=data.table(df2,key=c("chr","start","end","strand") )
                # using hacked data.table merge called merge2: temporary fix
                df=merge(df,df2,by=c("chr","start","end","strand"),suffixes=c(as.character(i-1),as.character(i) ) ,all=TRUE)
                #setkeyv(X,c("chr","start","end","strand"))
                #df=df[df2, nomatch=FALSE]
              }
              sample.ids=c(sample.ids,object[[i]]@sample.id)
              contexts=c(contexts,object[[i]]@context)
            }
            
            # stop if the assembly of object don't match
            if( length( unique(assemblies) ) != 1 ){stop("assemblies of methylrawList elements should be same\n")}
            
            
            if(  ! is.null(min.per.group) ){
              # if the the min.per.group argument is supplied, remove the rows that doesn't have enough coverage
              
              # get indices of coverage,numCs and numTs in the data frame 
              coverage.ind=seq(5,by=3,length.out=length(object))
              numCs.ind   =coverage.ind+1
              numTs.ind   =coverage.ind+2
              start.ind   =2 # will be needed to weed out NA values on chr/start/end/strand
              
              for(i in unique(object@treatment) ){
                my.ind=coverage.ind[object@treatment==i]
                ldat = !is.na(df[,my.ind,with=FALSE])
                if(  is.null(dim(ldat))  ){  # if there is only one dimension
                  df=df[ldat>=min.per.group,]
                }else{
                  df=df[rowSums(ldat)>=min.per.group,]
                }
              }

            }
            df=as.data.frame(df)
            # get indices of coverage,numCs and numTs in the data frame 
            coverage.ind=seq(5,by=3,length.out=length(object))
            numCs.ind   =coverage.ind+1
            numTs.ind   =coverage.ind+2
            
            # change column names
            names(df)[coverage.ind]=paste(c("coverage"),1:length(object),sep="" )
            names(df)[numCs.ind]   =paste(c("numCs"),1:length(object),sep="" )
            names(df)[numTs.ind]   =paste(c("numTs"),1:length(object),sep="" )
            
            if(!save.db){
              #make methylbase object and return the object
              obj=new("methylBase",(df),sample.ids=sample.ids,
                      assembly=unique(assemblies),context=unique(contexts),
                      treatment=object@treatment,coverage.index=coverage.ind,
                      numCs.index=numCs.ind,numTs.index=numTs.ind,destranded=destrand,resolution=object[[1]]@resolution )
            }else{
              dbdir <- .check.dbdir(dbdir)
              #save dataframe as methylbaseDB object and return object
              obj=makeMethylBaseDB(df = df,dbpath = dbdir,dbtype = "tabix",sample.ids=sample.ids,
                                   assembly=unique(assemblies),context=unique(contexts),
                                   treatment=object@treatment,coverage.index=coverage.ind,
                                   numCs.index=numCs.ind,numTs.index=numTs.ind,destranded=destrand,
                                   resolution=object[[1]]@resolution)
            }

            obj
          }
)           



#' @rdname getCorrelation-methods
#' @aliases getCorrelation,methylBaseDB-method
setMethod("getCorrelation", "methylBaseDB",
          function(object,method,plot){
            
            data = getData(object)
            meth.mat = data[, object@numCs.index]/
              ( data[,object@numCs.index] + data[,object@numTs.index] )                                      
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

#' @rdname getCoverageStats-methods
#' @aliases getCoverageStats,methylRawDB-method
setMethod("getCoverageStats", "methylRawDB",
          function(object,plot,both.strands,labels,...){
            
            if(!plot){
              qts=seq(0,0.9,0.1) # get quantiles
              qts=c(qts,0.95,0.99,0.995,0.999,1)                          
              
              if(both.strands){
                plus.cov=object[object[][["strand"]]=="+",][["coverage"]]
                mnus.cov=object[object[][["strand"]]=="-",][["coverage"]]
                
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
                
                all.cov=object[][["coverage"]]
                
                cat("read coverage statistics per base\n")
                cat("summary:\n")
                print( summary( all.cov ) )
                cat("percentiles:\n")
                print(quantile( all.cov,p=qts ))
                cat("\n")
              }
              
            }else{
              if(both.strands){   
                plus.cov=object[object[][["strand"]]=="+",][["coverage"]]
                mnus.cov=object[object[][["strand"]]=="-",][["coverage"]]
                
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
                all.cov= object[][["coverage"]]
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
            
            plus.met=100* object[object[][["strand"]]=="+",][["numCs"]]/object[object[][["strand"]]=="+",][["coverage"]]
            mnus.met=100* object[object[][["strand"]]=="-",][["numCs"]]/object[object[][["strand"]]=="-",][["coverage"]]
            all.met =100* object[][["numCs"]]/object[][["coverage"]]
            
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

#' @rdname getAssembly-methods
#' @aliases getAssembly,methylBaseDB-method
setMethod("getAssembly", signature="methylBaseDB", definition=function(x) {
  return(x@assembly)
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

#' @rdname getContext-methods
#' @aliases getContext,methylBaseDB-method
setMethod("getContext", signature="methylBaseDB", definition=function(x) {
  return(x@context)
})

#' @rdname getData-methods
#' @aliases getData,methylRawDB-method
setMethod("getData", signature="methylRawDB", definition=function(x) {
  df <- headTabix(tbxFile = x@dbpath, nrow = x@num.records, return.type = "data.frame")
  data.table::setnames(x = df,old = names(df), new = c("chr","start","end","strand","coverage","numCs","numTs"))
  return(df)
})

#' @rdname getData-methods
#' @aliases getData,methylRawListDB-method
setMethod("getData", signature="methylRawListDB", definition=function(x) {
  return(lapply(x,getData))
})

#' @rdname getData-methods
#' @aliases getData,methylBaseDB-method
setMethod("getData", signature="methylBaseDB", definition=function(x) {
  df <- headTabix(tbxFile = x@dbpath, nrow = x@num.records, return.type = "data.frame")
  data.table::setnames(x = df,old = c("V1","V2","V3","V4"), new = c("chr","start","end","strand"))
  # get indices of coverage,numCs and numTs in the data frame 
  numsamples = (length(df)-4)/3
  coverage.ind=seq(5,by=3,length.out=numsamples)
  numCs.ind   =coverage.ind+1
  numTs.ind   =coverage.ind+2
  
  # change column names
  names(df)[coverage.ind]=paste(c("coverage"),1:numsamples,sep="" )
  names(df)[numCs.ind]   =paste(c("numCs"),1:numsamples,sep="" )
  names(df)[numTs.ind]   =paste(c("numTs"),1:numsamples,sep="" )
  
  return(df)
})


# show functions ----------------------------------------------------------

setMethod("show", "methylRawDB", function(object) {
  
  cat("methylRawDB object with",object@num.records,"rows\n--------------\n")
  print(head(getData(object),n = 6))
  cat("--------------\n")
  cat("sample.id:",object@sample.id,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("resolution:", object@resolution,"\n")
  cat("dbtype:", object@dbtype,"\n")
  #cat("dbpath:",object@dbpath,"\n")
  cat("\n")
})

setMethod("show", "methylRawListDB", function(object) {
  
  cat("methylRawListDB object with",length(object),"methylRawDB objects\n\n")
  
  lapply(object,show)
  cat("treament:", object@treatment,"\n")
  
})

#' @rdname show-methods
#' @aliases show,methylBaseDB
setMethod("show", "methylBaseDB", function(object) {
  
  cat("methylBaseDB object with",object@num.records,"rows\n--------------\n")
  print(head(getData(object),n = 6))
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


#' set treatment vector of the methylRawListDB or methylBaseDB object
#' 
#' The function returns the \code{\link{methylBaseDB}},\code{\link{methylRawListDB}} objects with the new treatment vector stored in it.
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
#' @return the object with new treatment vector
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
#' @usage getSampleNames(x)
#' @examples
#' 
#' data(methylKit)
#' 
#' getSampleNames(methylBaseDB.obj)
#' getSampleNames(methylRawListDB.obj)
#' 
#' 
#' @return the sample ids stored in the object
#' @export
#' @docType methods
#' @rdname getSampleNames-methods
setGeneric("getSampleNames", def=function(x) standardGeneric("getSampleNames"))
#' @rdname getSampleNames-methods
#' @aliases getSampleNames,methylRawListDB-method
setMethod("getSampleNames", signature = "methylRawListDB", function(x) {
  names <- vapply(x,function(z) z@sample.id,FUN.VALUE = "character")
  return(names)
})

#' @rdname getSampleNames-methods
#' @aliases getSampleNames,methylBaseDB-method
setMethod("getSampleNames", signature = "methylBaseDB", function(x) {
  return(x@sample.ids)
})


#' set sample ids of the methylRawListDB or methylBaseDB object
#' 
#' The function returns the \code{\link{methylBaseDB}},\code{\link{methylRawListDB}} objects with the new sample ids stored in it.
#' 
#' @param x an \code{\link{methylBaseDB}} or \code{\link{methylRawListDB}} object
#' @usage setSampleNames(x) <- c("id.1","id.2","id.3",...,"id.n")
#' @examples
#' 
#' data(methylKit)
#' 
#' setSampleNames(methylBaseDB.obj) <- c("t1","t2","c1","c2")
#' setSampleNames(methylRawListDB.obj) <- c("test-one","test-two","control-one","control-two")
#' 
#' 
#' @return the object with new sample ids
#' @export
#' @docType methods
#' @rdname setSampleNames-methods
setGeneric("setSampleNames<-", def=function(x, value="numeric") {standardGeneric("setSampleNames<-")})
#' @rdname setSampleNames-methods
#' @aliases setSampleNames,methylRawDB-method
setReplaceMethod("setSampleNames", signature = "methylRawDB", function(x, value) {
  
  if(! ( length(value) == 1 ) ){
    message("The vector of new sample ids is not valid, check the length of input")
    return(x)
  } else {
    x@sample.id <- value
    return(x)
  }
  
})

#' @rdname setSampleNames-methods
#' @aliases setSampleNames,methylRawListDB-method
setReplaceMethod("setSampleNames", signature = "methylRawListDB", function(x, value) {
  
  if(! ( length(getSampleNames(x)) == length(value) ) ){
    message("The vector of new sample ids is not valid, check the length of input")
    return(x)
  } else {
    treatment <- x@treatment
    x <- mapply(`setSampleNames<-`, x, value)
    x <- new("methylRawListDB",x,treatment=treatment)
    return(x)
  }
  
})

#' @rdname setSampleNames-methods
#' @aliases setSampleNames,methylBaseDB-method
setReplaceMethod("setSampleNames", signature = "methylBaseDB", function(x, value) {
  
  if(! ( length(x@sample.ids) == length(value) ) ){
    message("The vector of new sample ids is not valid, check the length of input")
    return(x)
  } else {
    x@sample.ids <- value
    return(x)
  }
  
})



# MethylRawDB and MethylRawListDB -----------------------------------------

#' @aliases filterByCoverage,methylRawDB-method
#' @rdname filterByCoverage-methods
setMethod("filterByCoverage", signature(methylObj="methylRawDB"),
          function(methylObj,lo.count,lo.perc,hi.count,hi.perc,chunk.size,
                   save.db=TRUE,...){
  if( is.null(lo.count) & is.null(lo.perc) & is.null(hi.count) & 
      is.null(hi.perc) ){return(methylObj)}
  
  if(save.db) {
    
    filter <- function(data,lo.count,lo.perc,hi.count,hi.perc) {
      
      .setMethylDBNames(data,"methylRawDB")
      
      #figure out which cut-offs to use, maybe there is more elagent 
      # ways, quick&dirty works for now
      if(is.numeric(lo.count) ){lo.count=lo.count}
      if(is.numeric(lo.perc)){lo.count=quantile(data$coverage,lo.perc/100)}
      if(is.numeric(hi.count)){hi.count=hi.count}
      if(is.numeric(hi.perc)){hi.count=quantile(data$coverage,hi.perc/100)}
      
      if(is.numeric(lo.count)){data=data[data$coverage>=lo.count,]}
      if(is.numeric(hi.count)){data=data[data$coverage<hi.count,]}
      
      return(data)
      
    }
    
    # catch additional args 
    args <- list(...)
    dir <- dirname(methylObj@dbpath)
    
    if( "dbdir" %in% names(args) ){
      if( !(is.null(args$dbdir)) ){
        dir <- .check.dbdir(args$dbdir) 
      }
    } 
    
    if(!( "suffix" %in% names(args) ) ){
      suffix <- paste0("_","filtered")
    } else { 
      suffix <- paste0("_",args$suffix)
    }
    
    filename <- paste0(paste(methylObj@sample.id,collapse = "_"),suffix,".txt")
    #filename <- paste0(basename(gsub(".txt.bgz",replacement = "",methylObj@dbpath)),
    # suffix,".txt")
    
    #print(filename)
    
    newdbpath <- applyTbxByChunk(methylObj@dbpath,chunk.size = chunk.size, 
                                 dir=dir,filename = filename, 
                                 return.type = "tabix", FUN = filter, 
                                 lo.count=lo.count, lo.perc=lo.perc, 
                                 hi.count=hi.count, hi.perc=hi.perc)
    
    readMethylRawDB(dbpath = newdbpath,dbtype = "tabix",
                    sample.id = methylObj@sample.id,
                    assembly = methylObj@assembly, context = methylObj@context,
                    resolution = methylObj@resolution)
    
  } else {
    
    methylObj <- methylObj[]
    #print(class(methylObj))
    filterByCoverage(methylObj,lo.count,lo.perc,hi.count,hi.perc,
                     chunk.size,save.db=FALSE,...)
    
  }
            
})

#' @aliases filterByCoverage,methylRawListDB-method
#' @rdname filterByCoverage-methods
setMethod("filterByCoverage", signature(methylObj="methylRawListDB"),
          function(methylObj,lo.count,lo.perc,hi.count,hi.perc,chunk.size,
                   save.db=TRUE,...){
            
    if(save.db){
      args <- list(...)
      
      if( !( "dbdir" %in% names(args)) ){
        dbdir <- NULL
      } else { dbdir <- basename(.check.dbdir(args$dbdir)) }
      
      new.list=lapply(methylObj,filterByCoverage,lo.count,lo.perc,hi.count,
                      hi.perc,chunk.size,save.db,dbdir=dbdir,...)
      new("methylRawListDB", new.list,treatment=methylObj@treatment)
      
    } else {
      new.list=lapply(methylObj,filterByCoverage,lo.count,lo.perc,hi.count,
                      hi.perc,chunk.size,save.db=FALSE,...)
      new("methylRawList", new.list,treatment=methylObj@treatment)
    }
    
})


#' @rdname getCoverageStats-methods
#' @aliases getCoverageStats,methylRawDB-method
setMethod("getCoverageStats", "methylRawDB",
          function(object,plot,both.strands,labels,...,chunk.size){
            
            tmp = applyTbxByChunk(object@dbpath,chunk.size = chunk.size,
                                  return.type = "data.table", 
                                  FUN = function(x) { 
                                    .setMethylDBNames(x,"methylRawDB");
                                    return(x[,c("strand","coverage"),with=FALSE])} )
            
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
                  a=hist(log10(plus.cov),plot=FALSE)
                  my.labs=as.character(round(100*a$counts/length(plus.cov),1))
                }else{my.labs=FALSE}
                
                hist(log10(plus.cov),col="chartreuse4",
                     xlab=paste("log10 of read coverage per",object@resolution),
                     main=paste("Histogram of", object@context, 
                                "coverage: Forward strand"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
                if(labels){
                  a=hist(log10(mnus.cov),plot=FALSE)
                  my.labs=as.character(round(100*a$counts/length(mnus.cov),1))
                }else{my.labs=FALSE}
                a=hist(log10(mnus.cov),plot=FALSE)
                hist(log10(mnus.cov),col="chartreuse4",
                     xlab=paste("log10 of read coverage per",object@resolution),
                     main=paste("Histogram of", object@context, 
                                "coverage: Reverse strand"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
              }else{
                all.cov=tmp[,coverage]
                if(labels){
                  a=hist(log10(all.cov),plot=FALSE)
                  my.labs=as.character(round(100*a$counts/length(all.cov),1))
                }else{my.labs=FALSE}                          
                
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
          function(object,plot,both.strands,labels,...,chunk.size){
            
            numCs=coverage=strand=.=NULL
            tmp = applyTbxByChunk(object@dbpath,chunk.size = chunk.size,
                                  return.type = "data.table", 
                                  FUN = function(x) { 
                                    .setMethylDBNames(x,"methylRawDB"); 
                                    return(x[,.(strand,coverage,numCs)])} )
            
            
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
                  a=hist((plus.met),plot=FALSE,...)
                  my.labs=as.character(round(100*a$counts/length(plus.met),1))
                }else{my.labs=FALSE}
                hist((plus.met),col="cornflowerblue",
                     xlab=paste("% methylation per",object@resolution),
                     main=paste("Histogram of %", object@context,
                                "methylation: Forward strand"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
                if(labels){                          
                  a=hist((mnus.met),plot=FALSE,...)
                  my.labs=as.character(round(100*a$counts/length(mnus.met),1))
                }
                else{my.labs=FALSE}
                
                hist((mnus.met),col="cornflowerblue",
                     xlab=paste("% methylation per",object@resolution),
                     main=paste("Histogram of %", object@context,
                                "methylation: Reverse strand"),
                     labels=my.labs,...)
                mtext(object@sample.id, side = 3)
                
              }else{
                if(labels){                          
                  
                  a=hist((all.met),plot=FALSE,...)
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


#' @rdname adjustMethylC
#' @aliases adjustMethylC,methylRawDB,methylRawDB-method
setMethod("adjustMethylC", c("methylRawDB","methylRawDB"),
          function(mc,hmc,save.db=TRUE,...,chunk.size){
  
  if(save.db) {
    
    lst=new("methylRawListDB",list(mc,hmc),treatment=c(1,0))
    base=unite(lst)
    
    adjust <- function(data) {
      
      .setMethylDBNames(data,"methylBaseDB")
      diff=(data$numCs1)-round(data$coverage1*(data$numCs2/data$coverage2))
      diff[diff<0]=0
      data$numCs1=diff
      data$numTs1=data$coverage1-data$numCs1
      return(data[1:7])
      
    }
    
    # catch additional args 
    args <- list(...)
    dir <- dirname(mc@dbpath)
    
    if( "dbdir" %in% names(args) ){
      if( !(is.null(args$dbdir)) ){
        dir <- .check.dbdir(args$dbdir) 
      }
    } 
    
    if(!( "suffix" %in% names(args) ) ){
      suffix <- paste0("_","adjusted")
    } else { 
      suffix <- paste0("_",args$suffix)
    }
    
    
    filename <- paste0(paste(mc@sample.id,collapse = "_"),suffix,".txt")
    #filename <- paste0(basename(gsub(".txt.bgz",replacement = "",mc@dbpath)),suffix,".txt")
    
    newdbpath <- applyTbxByChunk(base@dbpath,chunk.size = chunk.size, 
                                 dir=dir,filename = filename, 
                                 return.type = "tabix", FUN = adjust)
    
    unlink(list.files(dirname(base@dbpath),
                      pattern = basename(gsub(".txt.bgz","",base@dbpath)),
                      full.names = TRUE))
    
    readMethylRawDB(dbpath = newdbpath,sample.id=mc@sample.id,
                    assembly=mc@assembly, context =mc@context,
                    resolution=mc@resolution,
                    dbtype = mc@dbtype)
    
  } else {
    
    mc.tmp <- mc[]
    hmc.tmp <- hmc[]
    adjustMethylC(mc.tmp,hmc.tmp,save.db=FALSE,...)
    
    
  }
  
  
})


#' @rdname adjustMethylC
#' @aliases adjustMethylC,methylRawListDB,methylRawListDB-method
setMethod("adjustMethylC", c("methylRawListDB","methylRawListDB"),
          function(mc,hmc,save.db=TRUE,...,chunk.size){
  
  # check lengths equal if not give error
  if(length(mc) != length(hmc)){stop("lengths of methylRawList objects should be same\n")}
  
  if(save.db){
    
    args <- list(...)
    
    if( !( "dbdir" %in% names(args)) ){
      dbdir <- NULL
    } else { dbdir <- .check.dbdir(args$dbdir) }
    
    my.list=list()
    for(i in 1:length(mc)){
      my.list[[i]]=adjustMethylC(mc[[i]],hmc[[i]],save.db,
                                  dbdir=basename(dbdir),...,chunk.size)
    }
    new("methylRawListDB",my.list,treatment=mc@treatment )
    
  } else {
    
    my.list=list()
    for(i in 1:length(mc)){
      my.list[[i]]=adjustMethylC(mc[[i]],hmc[[i]],save.db=FALSE,...)
    }
    new("methylRawList",my.list,treatment=mc@treatment )
    
  }
  
})

#' @rdname normalizeCoverage-methods
#' @aliases normalizeCoverage,methylRawListDB-method
setMethod("normalizeCoverage", "methylRawListDB",
          function(obj,method,chunk.size,save.db=TRUE,...){
            
            if( !(method %in% c("median","mean") ) )
            {
              stop("method option should be either 'mean' or 'median'\n")
            }
            
            if(save.db) { 
              
              normCov <- function(data,method) {
                
                .setMethylDBNames(data)
                if(method=="median"){
                  x=median(data$coverage)
                }else {
                  x=mean(data$coverage)
                }
                
                sc.fac=max(x)/x #get scaling factor
                
                all.cov=data$coverage
                fCs    =data$numCs/all.cov
                fTs    =data$numT/all.cov
                data$coverage=round(sc.fac*data$coverage)
                data$numCs   =round(data$coverage*fCs)
                data$numTs   =round(data$coverage*fTs)
                
                return(data)
              }
              
              # catch additional args 
              args <- list(...)
              
              if( ( "dbdir" %in% names(args))   ){
                #                 if(!(is.null(args$dbdir))) {
                dir <- .check.dbdir(args$dbdir)
                #                 }
              } else { 
                dir <- dirname(obj[[1]]@dbpath)
              }
              if(!( "suffix" %in% names(args) ) ){
                suffix <- paste0("_","normed")
              } else { 
                suffix <- paste0("_",args$suffix)
              }
              
              outList <- list()
              
              for(i in 1:length(obj)){
                
                filename <- paste0(paste(obj[[i]]@sample.id,collapse = "_"),
                                   suffix,".txt")
                #filename <- paste0(basename(gsub(".txt.bgz",replacement = "",
                #obj[[i]]@dbpath)),suffix,".txt")
                
                newdbpath <- applyTbxByChunk(obj[[i]]@dbpath,
                                             chunk.size = chunk.size, 
                                             dir=dir,filename = filename, 
                                             return.type = "tabix", 
                                             FUN = normCov,method = method)
                
                outList[[i]] <- readMethylRawDB(dbpath = newdbpath,
                                                sample.id=obj[[i]]@sample.id,
                                                assembly=obj[[i]]@assembly, 
                                                context =obj[[i]]@context,
                                                resolution=obj[[i]]@resolution,
                                                dbtype = obj[[i]]@dbtype)
                
              }
              
              new("methylRawListDB",outList,treatment=obj@treatment)
              
            } else {
              
              # coerce methylRawDB to methylRaw
              tmp.list <- lapply(obj,function(x) x[])
              tmp <- new("methylRawList",tmp.list,treatment=obj@treatment)
              normalizeCoverage(tmp,method,save.db=FALSE)
            }
})

#' @rdname reorganize-methods
#' @aliases reorganize,methylRawListDB-method
setMethod("reorganize", signature(methylObj="methylRawListDB"),
          function(methylObj,sample.ids,treatment,chunk.size,save.db=TRUE,...){
            
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
  
  if(save.db) {
    
    # catch additional args 
    args <- list(...)
    
    outList=list()
    # do not rename per default
    suffix <- NULL
    
    if( ("dbdir" %in% names(args)) | ( "suffix" %in% names(args) ) ) {
      if( "suffix" %in% names(args) ) { suffix <- paste0("_",args$suffix) }
      if( ( "dbdir" %in% names(args)) ){ 
        dir <- .check.dbdir(args$dbdir) 
        
        for(i in 1:length(sample.ids)){
          obj <- methylObj[[ col.ord[i]  ]]
          filename <- paste0(dir,"/",paste(obj@sample.ids,collapse = "_"),
                             suffix,".txt.bgz")
          #filename <- paste0(dir,"/",basename(gsub(".txt.bgz",replacement = 
          # "",obj@dbpath)),suffix,".txt.bgz")
          file.copy(obj@dbpath,filename)
          
          outList[[i]]=readMethylRawDB(dbpath = filename,
                                       dbtype = obj@dbtype,
                                       sample.id = obj@sample.id,
                                       assembly = obj@assembly, 
                                       context = obj@context, 
                                       resolution = obj@resolution)
        }
      } else {
        
        for(i in 1:length(sample.ids)){
          obj <- methylObj[[ col.ord[i]  ]]
          filename <- paste0(paste(obj@sample.ids,collapse = "_")
                             ,suffix,".txt.bgz")
          #filename <- paste0(gsub(".txt.bgz",replacement = "",obj@dbpath),
          # suffix,".txt.bgz")
          file.copy(obj@dbpath,filename)
          
          outList[[i]]=readMethylRawDB(dbpath = filename,
                                       dbtype = obj@dbtype,
                                       sample.id = obj@sample.id,
                                       assembly = obj@assembly, 
                                       context = obj@context, 
                                       resolution = obj@resolution)
        }
      }
    } else {
      
      for(i in 1:length(sample.ids)){
        outList[[i]]=methylObj[[ col.ord[i]  ]]
      }
    }
    
    new("methylRawListDB",outList,treatment=treatment)
    
  } else {
    
    outList=list()    
    for(i in 1:length(sample.ids)){
      outList[[i]]=methylObj[[ col.ord[i]  ]][]
      
    }
    
    new("methylRawList",outList,treatment=treatment)
  }
  
})


# MethylBaseDB ------------------------------------------------------------

#' @rdname unite-methods
#' @aliases unite,methylRawListDB-method
setMethod("unite", "methylRawListDB",
          function(object,destrand,min.per.group,chunk.size,mc.cores,save.db=TRUE,...){
            
if(save.db) {
  
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
  if( (!is.null(min.per.group)) &  ( ! is.integer( min.per.group ) )  )
  {
    stop("min.per.group should be an integer\ntry providing integers as 1L, 2L,3L etc.\n")
  }
  
  if(Sys.info()['sysname']=="Windows") {mc.cores = 1}
  # destrand single objects contained in methylRawListDB
  if(destrand) { 
    
    destrandFun <- function(obj){
      
      if(obj@resolution == "base") {
        dir <- dirname(obj@dbpath)
        filename <- paste(basename(gsub(".txt.bgz","",obj@dbpath)),
                          "destrand.txt",sep="_")
        # need to use .CpG.dinuc.unifyOld because output needs to be ordered
        newdbpath <- applyTbxByChunk(obj@dbpath,
                                     chunk.size = chunk.size, 
                                     dir=dir,filename = filename,
                                     return.type = "tabix", 
                                     FUN = function(x) { 
                                       .CpG.dinuc.unifyOld(.setMethylDBNames(x,
                                                            "methylRawDB") )})
        
        readMethylRawDB(dbpath = newdbpath,dbtype = "tabix",
                        sample.id = obj@sample.id,
                        assembly = obj@assembly, context = obj@context,
                        resolution = obj@resolution)
      }else {obj}
      
    }
    new.list=lapply(object,destrandFun)
    object <- new("methylRawListDB", new.list,treatment=object@treatment)
    
    on.exit(unlink(list.files(dirname(dbpath),pattern = "destrand",
                              full.names = TRUE)))    
  }
  #merge raw methylation calls together
  
  objList <- sapply(object,FUN = function(x) x@dbpath)
  
  args <- list(...)
  #print(args)
  if( ( "dbdir" %in% names(args)) )
  { 
    dir <- .check.dbdir(args$dbdir) 
  } else {
    dir <- dirname(object[[1]]@dbpath)
  }
  
  filename <- paste0(paste(getSampleID(object),collapse = "_"),".txt")
  
  if(is.null(min.per.group)) {
    
    dbpath <- mergeTabix(tabixList = objList ,dir = dir,
                         filename = filename,mc.cores = mc.cores) 
    dbpath <- gsub(".tbi","",dbpath)
    
  } else {
    # if the the min.per.group argument is supplied, remove the rows
    # that doesn't have enough coverage
    
    # keep rows with no matching in all samples  
    tmpPath <- mergeTabix(tabixList = objList ,dir = dir,
                          filename = paste0(filename,".tmp"),
                          mc.cores = mc.cores,all=TRUE) 
    tmpPath <- gsub(".tbi","",tmpPath)
    
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
    
    dbpath <- applyTbxByChunk(tbxFile = tmpPath,chunk.size = chunk.size, dir = dir, filename = filename, 
                              return.type = "tabix", FUN = filter,treatment=object@treatment,
                              coverage.ind=coverage.ind,min.per.group=min.per.group)
    
    unlink(list.files(dirname(tmpPath),pattern = basename(gsub(".txt.bgz","",tmpPath)),full.names = TRUE))
  }
  readMethylBaseDB(dbpath = dbpath,dbtype = object[[1]]@dbtype,
                   sample.ids = getSampleID(object),assembly = object[[1]]@assembly,
                   context = object[[1]]@context,resolution = object[[1]]@resolution,
                   treatment = object@treatment,destranded = destrand)
} else {
  obj <- object
  class(obj) <- "methylRawList"
  unite(obj,destrand,min.per.group,chunk.size,mc.cores,save.db=FALSE,...)
}
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
    meth.mat = data[, object@numCs.index]/( data[,object@numCs.index] + 
                                              data[,object@numTs.index] )
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
  
  panel.my.smooth2<-function(x, y, col = par("col"), 
                             bg = NA, pch = par("pch"), cex = 1, 
                             col.smooth = "darkgreen", span = 2/3, 
                             iter = 3, ...) 
  {
    par(new = TRUE)    #par(usr = c(usr[1:2], 0, 1.5) )
    smoothScatter(x, y,colramp=colorRampPalette(topo.colors(100)), bg = bg)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = col.smooth, ...)
    abline(lm(y[ok]~x[ok]), col="red")
  }
  
  panel.my.smooth<-function(x, y, col = par("col"), bg = NA, 
                            pch = par("pch"), cex = 0.3, 
                            col.smooth = "green", span = 2/3, iter = 3, ...) 
  {
    points(x, y, pch = 20, col = densCols(x,y,
                                          colramp=colorRampPalette(topo.colors(20))), 
           bg = bg, cex = 0.1)
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
            diag.panel=panel.hist,main=paste(object@context, object@resolution ,
                                             method,"cor.") )
    }
    if(method=="kendall")
    { pairs(meth.mat, 
            lower.panel=panel.my.smooth2, 
            upper.panel=panel.cor.kendall,
            diag.panel=panel.hist,main=paste(object@context, object@resolution ,
                                             method,"cor.") )
    }
    if(method=="pearson")
    { pairs(meth.mat, 
            lower.panel=panel.my.smooth2, 
            upper.panel=panel.cor.pearson,
            diag.panel=panel.hist,main=paste(object@context, object@resolution ,
                                             method,"cor.") )
    }
  }
}  
)

#' @rdname reconstruct-methods
#' @aliases reconstruct,methylBaseDB-method
setMethod("reconstruct",signature(mBase="methylBaseDB"), 
          function(methMat,mBase,chunk.size,save.db=TRUE,...){
  
  if(save.db){
    
    if(is.matrix(methMat) ) {
      
      # check if indeed methMat is percent methylation matrix
      if(max(methMat)<=1){
        warning("\nmake sure 'methMat' is percent methylation matrix (values between 0-100) \n")
      }
      
      # check if methMat is percent methylation matrix fitting to mBase  
      if(nrow(methMat) != mBase@num.records | ncol(methMat) != 
         length(mBase@numCs.index) ){
        stop("\nmethMat dimensions do not match number of samples\n",
             "and number of bases in methylBase object\n")
      }
      
      rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
      matFile=paste(rndFile,"methMat.txt",sep="_")
      write.table(methMat,matFile,quote=FALSE,col.names=FALSE,row.names=FALSE,
                  sep="\t")
      methMat = matFile
      
    } else {
      
      # methMat can also be a file containing a percent methylation matrix
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
      rm(mat)
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
    
    # catch additional args 
    args <- list(...)
    
    if( !( "dbdir" %in% names(args)) ){
      dir <- dirname(mBase@dbpath)
    } else { dir <- .check.dbdir(args$dbdir) }
    if(!( "suffix" %in% names(args) ) ){
      suffix <- "_reconstructed"
    } else { 
      suffix <- paste0("_",args$suffix)
    }
    
    filename <- paste0(paste(mBase@sample.ids,collapse = "_"),suffix,".txt")
    #filename <- paste0(basename(gsub(".txt.bgz","",mBase@dbpath)),suffix,".txt")
    
    con <- file(methMat,open = "r") 
    
    newdbpath <- applyTbxByChunk(tbxFile = mBase@dbpath,chunk.size = chunk.size, 
                                 dir=dir,filename = filename, 
                                 return.type = "tabix", FUN = reconstr, 
                                 con,chunk.size,mBase@numCs.index,
                                 mBase@numTs.index)
    
    close(con)
    if(file.exists(matFile)) {unlink(matFile)}
    
    readMethylBaseDB(dbpath = newdbpath,dbtype = mBase@dbtype,
                     sample.ids = mBase@sample.ids,assembly = mBase@assembly,
                     context = mBase@context,resolution = mBase@resolution,
                     treatment = mBase@treatment,destranded = mBase@destranded)
  } else {
    
    tmp <- mBase[]
    reconstruct(methMat,tmp,save.db=FALSE)
    
  }
}
)

#' @rdname removeComp-methods
#' @aliases removeComp,methylBaseDB-method
setMethod("removeComp",signature(mBase="methylBaseDB"), 
          function(mBase,comp,chunk.size,save.db=TRUE,...){
  if(is.na(comp) || is.null(comp)){
    stop("no component to remove\n")
  }
  
  if(any(comp > length(mBase@sample.ids) )){
    stop("'comp' elements can only take values between 1 and number of samples\n")
  }
  
  scale=TRUE
  center=TRUE
  mat=percMethylation(mBase,chunk.size=chunk.size)  
  mat=scale(mat,scale=scale,center=center)
  centers <- attr(mat,'scaled:center') 
  scales <- attr(mat,'scaled:scale') 
  pr=prcomp(mat,scale.=FALSE,center=FALSE)
  
  pr$rotation[,comp]=0
  res=pr$x %*% t(pr$rotation)
  
  res=(scale(res,center=(-centers/scales),scale=1/scales))
  attr(res,"scaled:center")<-NULL 
  attr(res,"scaled:scale")<-NULL 
  res[res>100]=100
  res[res<0]=0
  reconstruct(res,mBase,chunk.size,save.db = save.db,...=...)
}
)

#' @rdname percMethylation-methods
#' @aliases percMethylation,methylBaseDB-method
setMethod("percMethylation", "methylBaseDB",
          function(methylBase.obj,rowids=FALSE,save.txt,chunk.size){
            
            meth.fun <- function(data, numCs.index, numTs.index){
              
              dat=100 * data[, numCs.index]/( data[,numCs.index] + 
                                                data[,numTs.index] )
              rownames(dat)=paste(as.character(data[,1]),
                                  data[,2],data[,3],sep=".")
            }
            if (save.txt) {
              
              filename <- paste0(basename(gsub(".txt.bgz","",
                                      methylBase.obj@dbpath)),"_methMath.txt")
              
              meth.mat = applyTbxByChunk(methylBase.obj@dbpath,
                                         return.type = "text", 
                                         chunk.size = chunk.size,
                                         dir = dirname(methylBase.obj@dbpath), 
                                         filename = filename,
                                         FUN = meth.fun, 
                                         numCs.index = methylBase.obj@numCs.index,
                                         numTs.index = methylBase.obj@numTs.index)
              return(meth.mat)
              
            } else {
              
              meth.mat = applyTbxByChunk(methylBase.obj@dbpath,
                                         return.type = "data.frame", 
                                         chunk.size = chunk.size,
                                         FUN = meth.fun, 
                                         numCs.index = methylBase.obj@numCs.index,
                                         numTs.index = methylBase.obj@numTs.index)
              
              names(meth.mat)=methylBase.obj@sample.ids
              if(!rowids){
                rownames(meth.mat)=NULL 
              }
              return(as.matrix(meth.mat))
              
            }
})

#' @rdname clusterSamples-methods
#' @aliases clusterSamples,methylBaseDB-method
setMethod("clusterSamples", "methylBaseDB",
          function(.Object, dist, method ,sd.filter, sd.threshold, 
                   filterByQuantile, plot,chunk.size)
          {
            
            getMethMat <- function(mat,numCs.index,numTs.index,sd.filter, 
                                   sd.threshold, filterByQuantile){
              
              # remove rows containing NA values, they might be 
              # introduced at unite step
              mat      =mat[ rowSums(is.na(mat))==0, ] 
              
              meth.mat = mat[, numCs.index]/
                (mat[,numCs.index] + mat[,numTs.index] )                                      
              
              # if Std. Dev. filter is on remove rows with low variation
              if(sd.filter){
                if(filterByQuantile){
                  sds=rowSds(as.matrix(meth.mat))
                  cutoff=quantile(sds,sd.threshold)
                  meth.mat=meth.mat[sds>cutoff,]
                }else{
                  meth.mat=meth.mat[rowSds(as.matrix(meth.mat))>sd.threshold,]
                }
              }
            }
            
            meth.mat <- applyTbxByChunk(.Object@dbpath,chunk.size = chunk.size,
                                        return.type = "data.frame",
                                        FUN=getMethMat,
                                        numCs.ind=.Object@numCs.index,
                                        numTs.ind=.Object@numTs.index,
                                        sd.filter=sd.filter, 
                                        sd.threshold=sd.threshold, 
                                        filterByQuantile=filterByQuantile)
            
            names(meth.mat)=.Object@sample.ids
            
            .cluster(meth.mat, dist.method=dist, hclust.method=method, 
                     plot=plot, treatment=.Object@treatment,
                     sample.ids=.Object@sample.ids,
                     context=.Object@context)
            
          }
)

#' @rdname PCASamples-methods
#' @aliases PCASamples,methylBaseDB-method
setMethod("PCASamples", "methylBaseDB",
          function(.Object, screeplot, adj.lim,scale,center,comp,
                   transpose,sd.filter, sd.threshold, 
                   filterByQuantile,obj.return,chunk.size)
{
  
  getMethMat <- function(mat,numCs.index,numTs.index,sd.filter, 
                         sd.threshold, filterByQuantile){
    
    # remove rows containing NA values, they might be introduced at unite step
    mat      =mat[ rowSums(is.na(mat))==0, ] 
    
    meth.mat = mat[, numCs.index]/
      (mat[,numCs.index] + mat[,numTs.index] )                                      
    
    
    # if Std. Dev. filter is on remove rows with low variation
    if(sd.filter){
      if(filterByQuantile){
        sds=rowSds(as.matrix(meth.mat))
        cutoff=quantile(sds,sd.threshold)
        meth.mat=meth.mat[sds>cutoff,]
      }else{
        meth.mat=meth.mat[rowSds(as.matrix(meth.mat))>sd.threshold,]
      }
    }
    
  }
  
  meth.mat <- applyTbxByChunk(.Object@dbpath,chunk.size = chunk.size,
                              return.type = "data.frame",FUN=getMethMat,
                              numCs.ind=.Object@numCs.index,
                              numTs.ind=.Object@numTs.index,
                              sd.filter=sd.filter, sd.threshold=sd.threshold, 
                              filterByQuantile=filterByQuantile)
  
  names(meth.mat)=.Object@sample.ids
  
  if(transpose){
    .pcaPlotT(meth.mat,comp1=comp[1],comp2=comp[2],screeplot=screeplot, 
              adj.lim=adj.lim, 
              treatment=.Object@treatment,sample.ids=.Object@sample.ids,
              context=.Object@context
              ,scale=scale,center=center,obj.return=obj.return)
    
  }else{
    .pcaPlot(meth.mat,comp1=comp[1],comp2=comp[2],screeplot=screeplot, 
             adj.lim=adj.lim, 
             treatment=.Object@treatment,sample.ids=.Object@sample.ids,
             context=.Object@context,
             scale=scale,center=center,  obj.return=obj.return)
  }
  
}      
)

#' @rdname reorganize-methods
#' @aliases reorganize,methylBaseDB-method
setMethod("reorganize", signature(methylObj="methylBaseDB"),
          function(methylObj,sample.ids,treatment,chunk.size,save.db=TRUE,...){
            
  #sample.ids length and treatment length should be equal
  if(length(sample.ids) != length(treatment) ){
    stop("length of sample.ids should be equal to treatment")
  }
  
  if( ! all(sample.ids %in% methylObj@sample.ids) ){
    stop("provided sample.ids is not a subset of the sample ids of the object")
  }
  
  if(save.db) {
    
    temp.id = methylObj@sample.ids # get the subset of ids
    # get the column order in the original matrix
    col.ord = order(match(temp.id,sample.ids))[1:length(sample.ids)] 
    
    # make a matrix indices for easy access 
    ind.mat=rbind(methylObj@coverage.index[col.ord],  
                  methylObj@numCs.index[col.ord],
                  methylObj@numTs.index[col.ord])
    
    
    getSub <- function(data,ind.mat) {
      
      newdat =data[,1:4]
      for(i in 1:ncol(ind.mat))
      {
        newdat=cbind(newdat,data[,ind.mat[,i]])
      }
      return(newdat)
    }
    
    # catch additional args 
    args <- list(...)
    
    if( ( "dbdir" %in% names(args))   ){
       dir <- .check.dbdir(args$dbdir) 
    } else { dir <- dirname(methylObj@dbpath) }
    
    if(!( "suffix" %in% names(args) ) ){
      suffix <- NULL
    } else { 
      suffix <- paste0("_",args$suffix)
    }
    
    filename <- paste0(paste(sample.ids,collapse = "_"),suffix,".txt")
    # filename <- paste0(basename(gsub(".txt.bgz",replacement = "",methylObj@dbpath)),suffix,".txt")
    
    newdbpath <- applyTbxByChunk(tbxFile = methylObj@dbpath,
                                 chunk.size = chunk.size, dir=dir,
                                 filename = filename, 
                                 return.type = "tabix", FUN = getSub, 
                                 ind.mat=ind.mat) 
    
    readMethylBaseDB(dbpath = newdbpath,dbtype = methylObj@dbtype,
                     sample.ids=sample.ids,
                     assembly=methylObj@assembly,context=methylObj@context,
                     treatment=treatment,destranded=methylObj@destranded, 
                     resolution=methylObj@resolution )
  } else {
    
    obj <- methylObj[]
    reorganize(obj,sample.ids,treatment,save.db=FALSE,...)
    
  }
})

#' @rdname pool-methods
#' @aliases pool,methylBaseDB-method
setMethod("pool", "methylBaseDB",
          function(obj,sample.ids,chunk.size,save.db,...){
            
  if(save.db) {
    
    mypool <- function(df,treatment,numCs.index){
      
      treat=unique(treatment)
      res=df[,1:4]
      for(i in 1:length(treat) ){
        
        # get indices
        setCs=numCs.index[treatment==treat[i]]
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
      }
      return(res)
    }
    
    treat = unique(obj@treatment)
    coverage.ind=3*(1:length(treat)) + 2
    
    # catch additional args 
    args <- list(...)
    
    if( ( "dbdir" %in% names(args))   ){
      if( !(is.null(args$dbdir)) ) { 
        dir <- .check.dbdir(args$dbdir) }
    } else { dir <- dirname(obj@dbpath) }
    if(!( "suffix" %in% names(args) ) ){
      suffix <- NULL
    } else { 
      suffix <- paste0("_",args$suffix)
    }
    
    filename <- paste0(paste(sample.ids,collapse = "_"),suffix,".txt")
    
    newdbpath <- applyTbxByChunk(tbxFile = obj@dbpath,chunk.size = chunk.size,
                                 dir=dir,filename = filename, 
                                 return.type = "tabix", FUN = mypool, 
                                 treatment = obj@treatment,
                                 numCs.index = obj@numCs.index) 
    
    readMethylBaseDB(dbpath = newdbpath,dbtype = obj@dbtype,
                     sample.ids=sample.ids,
                     assembly=obj@assembly,context=obj@context,
                     treatment=treat,destranded=obj@destranded,
                     resolution=obj@resolution )
  } else {
    
    tmp <- obj[]
    pool(tmp,sample.ids,save.db=FALSE)
  }
})




# MethylDiffDB ------------------------------------------------------------

#' @aliases calculateDiffMeth,methylBaseDB-method
#' @rdname calculateDiffMeth-methods
setMethod("calculateDiffMeth", "methylBaseDB",
          function(.Object,covariates,overdispersion=c("none","MN","shrinkMN"),
                   adjust=c("SLIM","holm","hochberg","hommel","bonferroni",
                            "BH","BY","fdr","none","qvalue"),
                   effect=c("wmean","mean","predicted"),parShrinkMN=list(),
                   test=c("F","Chisq"),mc.cores=1,slim=TRUE,weighted.mean=TRUE,
                   chunk.size,save.db=TRUE,...){
            
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
    if(!is.null(covariates)){if(ncol(covariates)+2 >= length(.Object@numTs.index))
    {
      stop("Too many covariates/too few replicates.")}
    }
    
    if(save.db) {
      
      # add backwards compatibility with old parameters
      if(slim==FALSE) adjust="BH" else adjust=adjust
      if(weighted.mean==FALSE) effect="mean" else effect=effect
      
      vars <- covariates
      
      # function to apply the test to data
      diffMeth <- function(data,Ccols,Tcols,formula,vars,treatment,
                           overdispersion,effect,
                           parShrinkMN,test,adjust,mc.cores){
        
        cntlist=split(as.matrix(data[,c(Ccols,Tcols)]),1:nrow(data))
        
        tmp=simplify2array(
          mclapply(cntlist,logReg,formula,vars,treatment=treatment,
                   overdispersion=overdispersion,effect=effect,
                   parShrinkMN=parShrinkMN,test=test,mc.cores=mc.cores))
        tmp <- as.data.frame(t(tmp))
        #print(head(tmp))
        x=data.frame(data[,1:4],tmp$p.value,
                     p.adjusted(tmp$q.value,method=adjust),
                     meth.diff=tmp$meth.diff.1,stringsAsFactors=FALSE)
        
        return(x)
      }
      
      # catch additional args 
      args <- list(...)
      dir <- dirname(.Object@dbpath)
      
      if( ( "dbdir" %in% names(args))   ){
        if( !(is.null(args$dbdir)) ) { 
          dir <- .check.dbdir(args$dbdir) 
        }
      }
      if(!( "suffix" %in% names(args) ) ){
        suffix <- "_diffMeth"
      } else { 
        suffix <- paste0("_",args$suffix)
      }
      
      filename <- paste0(paste(.Object@sample.ids,collapse = "_"),suffix,".txt")
      
      #filename <- paste(basename(gsub(".txt.bgz","",.Object@dbpath)),suffix,".txt")
      
      dbpath <- applyTbxByChunk(.Object@dbpath,dir = 
                                  dir,chunk.size = chunk.size,  
                                filename = filename, 
                                return.type = "tabix", FUN = diffMeth,
                                Ccols = .Object@numCs.index,
                                Tcols = .Object@numTs.index,
                                formula=formula,vars=vars,
                                treatment=.Object@treatment,
                                overdispersion=overdispersion,effect=effect,
                                parShrinkMN=parShrinkMN,test=test,
                                adjust=adjust,mc.cores=mc.cores)
      
      obj=readMethylDiffDB(dbpath = dbpath,dbtype = .Object@dbtype, 
                           sample.ids=.Object@sample.ids,
                           assembly=.Object@assembly,context=.Object@context,
                           destranded=.Object@destranded,
                           treatment=.Object@treatment,
                           resolution=.Object@resolution)
      obj
      
    } else {
      
      tmp <- .Object[]
      calculateDiffMeth(tmp,covariates,overdispersion=overdispersion,
                        adjust=adjust,
                        effect=effect,parShrinkMN=parShrinkMN,
                        test=test,mc.cores=mc.cores,
                        slim=slim,weighted.mean=weighted.mean,
                        save.db=FALSE)
      
    }
}
)

#' @aliases getMethylDiff,methylDiffDB-method
#' @rdname getMethylDiff-methods
setMethod(f="getMethylDiff", signature="methylDiffDB", 
          definition=function(.Object,difference,qvalue,type,
                              chunk.size,save.db=TRUE,...) {
            
  if(!( type %in% c("all","hyper","hypo") )){
    stop("Wrong 'type' argument supplied for the function, ",
         "it can be 'hypo', 'hyper' or 'all' ")
  }
  meth.diff=NULL
  if(save.db) {
    
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
    
    # catch additional args 
    args <- list(...)
    dir <- dirname(.Object@dbpath)
    
    if( ( "dbdir" %in% names(args))   ){
      if( !(is.null(args$dbdir)) ) { 
        dir <- .check.dbdir(args$dbdir) 
      } 
    }
    
    if(!( "suffix" %in% names(args) ) ){
      suffix <- paste0("_diffMeth_",type)
    } else { 
      suffix <- paste0("_",args$suffix)
    }
    
    
    filename <- paste0(paste(.Object@sample.ids,collapse = "_"),suffix,".txt")
    
    #filename <- paste(basename(gsub(".txt.bgz","",.Object@dbpath)),suffix,".txt")
    
    dbpath <- applyTbxByChunk(.Object@dbpath,chunk.size = chunk.size, 
                              dir = dir, filename = filename, 
                              return.type = "tabix", FUN = f,
                              difference = difference, qv = qvalue, type = type)
    
    
    obj=readMethylDiffDB(dbpath = dbpath,dbtype = .Object@dbtype, 
                         sample.ids=.Object@sample.ids,
                         assembly=.Object@assembly,context=.Object@context,
                         destranded=.Object@destranded,
                         treatment=.Object@treatment,resolution=.Object@resolution)
    return(obj)
    
  } else {
    
    tmp <- .Object[]
    getMethylDiff(tmp,difference,qvalue,type,save.db=FALSE)
    
  }
    
}) 

#' @aliases diffMethPerChr,methylDiffDB-method
#' @rdname  diffMethPerChr-methods
setMethod("diffMethPerChr", signature(x = "methylDiffDB"),
          function(x,plot,qvalue.cutoff, meth.cutoff,exclude,...){
            

diffMeth <- function(data,qvalue.cutoff, meth.cutoff){
  
  .setMethylDBNames(data,"methylDiffDB")
  
  temp.hyper=data[data$qvalue < qvalue.cutoff & data$meth.diff >= meth.cutoff,]
  temp.hypo =data[data$qvalue < qvalue.cutoff & data$meth.diff <= -meth.cutoff,]
  
  # plot barplot for percentage of DMCs per chr
  dmc.hyper.chr=merge(as.data.frame(table(temp.hyper$chr)), 
                      as.data.frame(table(data$chr)),by="Var1")
  dmc.hyper.chr=cbind(dmc.hyper.chr,perc=100*dmc.hyper.chr[,2]/
                        dmc.hyper.chr[,3])
  
  dmc.hypo.chr=merge(as.data.frame(table(temp.hypo$chr)), 
                     as.data.frame(table(data$chr)),by="Var1")
  dmc.hypo.chr=cbind(dmc.hypo.chr,perc=100*dmc.hypo.chr[,2]/dmc.hypo.chr[,3])
  
  dmc.hyper.hypo=merge(dmc.hyper.chr[,c(1,2,4)],dmc.hypo.chr[,c(1,2,4)],
                       by="Var1",all=TRUE) # merge hyper hypo per chromosome
  
  
  names(dmc.hyper.hypo)=c("chr","number.of.hypermethylated",
                          "percentage.of.hypermethylated",
                          "number.of.hypomethylated",
                          "percentage.of.hypomethylated")
  
  return(dmc.hyper.hypo)
}

res <- applyTbxByChr(x@dbpath,return.type = "data.frame",FUN = diffMeth,qvalue.cutoff=qvalue.cutoff, meth.cutoff=meth.cutoff)

dmc.hyper=100*sum(res$number.of.hypermethylated)/x@num.records # get percentages of hypo/ hyper
dmc.hypo =100*sum(res$number.of.hypomethylated)/x@num.records

all.hyper.hypo=data.frame(number.of.hypermethylated=sum(res$number.of.hypermethylated),
                          percentage.of.hypermethylated=dmc.hyper,
                          number.of.hypomethylated=sum(res$number.of.hypomethylated),
                          percentage.of.hypomethylated=dmc.hypo)

dmc.hyper.hypo=res[order(as.numeric(sub("chr","",res$chr))),] # order the chromosomes

if(plot){
  
  if(!is.null(exclude)){dmc.hyper.hypo=
    dmc.hyper.hypo[! dmc.hyper.hypo$chr %in% exclude,]}
  
  barplot(
    t(as.matrix(data.frame(hyper=dmc.hyper.hypo[,3],
                           hypo=dmc.hyper.hypo[,5],
                           row.names=dmc.hyper.hypo[,1]) ))
    ,las=2,horiz=TRUE,col=c("magenta","aquamarine4"),
    main=paste("% of hyper & hypo methylated regions per chromosome",sep=""),
    xlab="% (percentage)",...)
  mtext(side=3,paste("qvalue<",qvalue.cutoff," & methylation diff. >=",
                     meth.cutoff," %",sep="") )
  legend("topright",legend=c("hyper","hypo"),fill=c("magenta","aquamarine4"))
}else{
  
  list(diffMeth.per.chr=dmc.hyper.hypo,diffMeth.all=all.hyper.hypo)
  
}
            
})

#' @aliases annotateWithGenicParts,methylDiffDB,GRangesList-method
#' @rdname annotateWithGenicParts-methods
setMethod("annotateWithGenicParts", 
          signature(target = "methylDiffDB",GRangesList.obj="GRangesList"),
          function(target,GRangesList.obj,strand){
            gr=as(target,"GRanges")
            annotateWithGenicParts(gr,GRangesList.obj,strand)
          })

#' @aliases annotateWithFeatureFlank,methylDiffDB,GRanges,GRanges-method
#' @rdname annotateWithFeatureFlank-methods
setMethod("annotateWithFeatureFlank", 
          signature(target= "methylDiffDB",
                                                feature="GRanges",flank="GRanges"),
          function(target, feature, flank,feature.name,flank.name,strand){
  gr=as(target,"GRanges")
  annotateWithFeatureFlank(gr,feature, flank,feature.name,flank.name,strand)
})

#' @aliases annotateWithFeature,methylDiffDB,GRanges-method
#' @rdname annotateWithFeature-methods
setMethod("annotateWithFeature", signature(target = "methylDiffDB",feature="GRanges"),
          function(target, feature, strand,extend,feature.name){                      
            gr=as(target,"GRanges")
            annotateWithFeature(gr, feature, strand,extend,feature.name)
          })

# Regionalize methods ----------------------------------------------------

#' @rdname regionCounts
#' @aliases regionCounts,methylRawDB,GRanges-method
setMethod("regionCounts", signature(object="methylRawDB",regions="GRanges"),
          function(object,regions,cov.bases,
                   strand.aware,chunk.size,save.db=TRUE,...){
            
  if(save.db) {
    
    # sort regions
    regions <- sortSeqlevels(regions)
    regions <- sort(regions,ignore.strand=TRUE)
    
    getCounts <- function(data,regions,cov.bases,strand.aware){
      .setMethylDBNames(data)
      # overlap data with regions
      # convert data to GRanges without metacolumns
      g.meth=with(data, GRanges(chr, IRanges(start=start, end=end),strand =strand))
      if(!strand.aware){
        strand(g.meth)="*"
        mat=IRanges::as.matrix( findOverlaps(regions,g.meth ) )
        #mat=matchMatrix( findOverlaps(regions,g.meth ) )
        
      }else{
        mat=IRanges::as.matrix( findOverlaps(regions,g.meth) )
        #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
        
      }
      
      # create a temporary data.table row ids from regions and counts from object
      temp.dt=data.table(id = mat[, 1], data[mat[, 2], c(5, 6, 7)])
      
      coverage=NULL
      numCs=NULL
      numTs=NULL
      id=covered=NULL
      
      # use data.table to sum up counts per region
      sum.dt=temp.dt[,list(coverage=sum(coverage),
                           numCs   =sum(numCs),
                           numTs   =sum(numTs),covered=length(numTs)),by=id] 
      sum.dt=sum.dt[covered>=cov.bases,]
      temp.df=as.data.frame(regions) # get regions to a dataframe
      
      #create a new methylRaw object to return
      new.data=data.frame(#id      =new.ids,
        chr     =temp.df[sum.dt$id,"seqnames"],
        start   =temp.df[sum.dt$id,"start"],
        end     =temp.df[sum.dt$id,"end"],
        strand  =temp.df[sum.dt$id,"strand"],
        coverage=sum.dt$coverage,
        numCs   =sum.dt$numCs,
        numTs   =sum.dt$numTs)
    }
    
    # catch additional args 
    args <- list(...)
    dir <- dirname(object@dbpath)
    
    if( "dbdir" %in% names(args) ){
      if( !(is.null(args$dbdir)) ){
        dir <- .check.dbdir(args$dbdir) 
      }
    } 
    if(!( "suffix" %in% names(args) ) ){
      suffix <- paste0("_","regions")
    } else { 
      suffix <- paste0("_",args$suffix)
    }
    
    filename <- paste0(paste(object@sample.id,collapse = "_"),suffix,".txt")
    
    newdbpath <- applyTbxByOverlap(object@dbpath,chunk.size = chunk.size, 
                                   ranges=regions,  dir=dir,filename = filename, 
                                   return.type = "tabix", FUN = getCounts,
                                   regions,cov.bases,strand.aware)
    
    readMethylRawDB(dbpath = newdbpath,sample.id=object@sample.id,
                    assembly=object@assembly, context =object@context,
                    resolution="region",
                    dbtype = object@dbtype)
    
  } else {
    tmp <- object[]
    regionCounts(tmp,regions,cov.bases,strand.aware,save.db=FALSE)
  }
}
)

#' @rdname regionCounts
#' @aliases regionCounts,methylRawDB,GRangesList-method
# assume that each name of the element in the GRangesList is unique and 
setMethod("regionCounts", signature(object="methylRawDB",regions="GRangesList"),
          function(object,regions,cov.bases,strand.aware,chunk.size,save.db=TRUE,...){
            
            # combine and sort GRanges from List
            regions <- unlist(regions)
            regions <- sortSeqlevels(regions)
            regions <- sort(regions,ignore.strand=TRUE)
            regions <- unique(regions)
            
            if(save.db) {
              
              getCounts <- function(data,regions,cov.bases,strand.aware){
                
                .setMethylDBNames(data)
                
                # overlap data with regions
                # convert data to GRanges without metacolumns
                g.meth=with(data, GRanges(chr, IRanges(start=start, end=end),strand =strand))
                
                if(!strand.aware){
                  strand(g.meth)="*"
                  mat=IRanges::as.matrix( findOverlaps(regions,g.meth ) )
                  #mat=matchMatrix( findOverlaps(regions,g.meth ) )
                  
                }else{
                  mat=IRanges::as.matrix( findOverlaps(regions,g.meth) )
                  #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
                  
                }
                
                # create a temporary data.table row ids from regions and counts from object
                temp.dt=data.table(id = mat[, 1], data[mat[, 2], c(5, 6, 7)])
                
                coverage=NULL
                numCs=NULL
                numTs=NULL
                id=covered=NULL
                
                # use data.table to sum up counts per region
                sum.dt=temp.dt[,list(coverage=sum(coverage),
                                     numCs   =sum(numCs),
                                     numTs   =sum(numTs),covered=length(numTs)),by=id] 
                sum.dt=sum.dt[covered>=cov.bases,]
                temp.df=as.data.frame(regions) # get regions to a dataframe
                
                #create a new methylRaw object to return
                new.data=data.frame(#id      =new.ids,
                  chr     =temp.df[sum.dt$id,"seqnames"],
                  start   =temp.df[sum.dt$id,"start"],
                  end     =temp.df[sum.dt$id,"end"],
                  strand  =temp.df[sum.dt$id,"strand"],
                  coverage=sum.dt$coverage,
                  numCs   =sum.dt$numCs,
                  numTs   =sum.dt$numTs)
              }
              
              
              
              # catch additional args 
              args <- list(...)
              dir <- dirname(object@dbpath)
              
              if( "dbdir" %in% names(args) ){
                if( !(is.null(args$dbdir)) ){
                  dir <- .check.dbdir(args$dbdir) 
                }
              } 
              if(!( "suffix" %in% names(args) ) ){
                suffix <- paste0("_","regions")
              } else { 
                suffix <- paste0("_",args$suffix)
              }
              
              
              filename <- paste0(paste(object@sample.id,collapse = "_"),suffix,".txt")
              
              newdbpath <- applyTbxByOverlap(object@dbpath,chunk.size = chunk.size, ranges=regions,  dir=dir,filename = filename, 
                                             return.type = "tabix", FUN = getCounts,regions,cov.bases,strand.aware)
              
              
              readMethylRawDB(dbpath = newdbpath,sample.id=object@sample.id,
                              assembly=object@assembly, context =object@context,resolution="region",
                              dbtype = object@dbtype)
              
            } else {
              tmp <- object[]
              regionCounts(tmp,regions,cov.bases,strand.aware,save.db=FALSE)
            }
            
          }
)

#' @rdname regionCounts
#' @aliases regionCounts,methylRawListDB,GRanges-method
setMethod("regionCounts", signature(object="methylRawListDB",regions="GRanges"),
          function(object,regions,cov.bases,strand.aware,chunk.size,save.db=TRUE,...){
            
            if(save.db){
              args <- list(...)
              
              if( !( "dbdir" %in% names(args)) ){
                dbdir <- NULL
              } else { dbdir <- basename(.check.dbdir(args$dbdir)) }
              
              outList = lapply(object,regionCounts,regions,cov.bases, strand.aware,
                               chunk.size,save.db,dbdir=dbdir,...)
              new("methylRawListDB", outList,treatment=object@treatment)
              
            } else {
              outList = lapply(object,regionCounts,regions,cov.bases,strand.aware,
                               chunk.size,save.db=FALSE)
              new("methylRawList", outList,treatment=object@treatment)
            }
          }
)

#' @rdname regionCounts
#' @aliases regionCounts,methylRawListDB,GRangesList-method
setMethod("regionCounts", signature(object="methylRawListDB",
                                    regions="GRangesList"),
          function(object,regions,cov.bases,strand.aware,chunk.size,save.db=TRUE,...){
            
            if(save.db){
              args <- list(...)
              
              if( !( "dbdir" %in% names(args)) ){
                dbdir <- NULL
              } else { dbdir <- basename(.check.dbdir(args$dbdir)) }
              
              outList = lapply(object,regionCounts,regions,cov.bases, strand.aware,
                               chunk.size,save.db,dbdir=dbdir,...)
              new("methylRawListDB", outList,treatment=object@treatment)
              
            } else {
              outList = lapply(object,regionCounts,regions,cov.bases,strand.aware,
                               chunk.size,save.db=FALSE)
              new("methylRawList", outList,treatment=object@treatment)
            }
          }
)

#' @rdname regionCounts
#' @aliases regionCounts,methylBaseDB,GRanges-method
setMethod("regionCounts", signature(object="methylBaseDB",regions="GRanges"),
          function(object,regions,cov.bases,strand.aware,chunk.size,
                   save.db=TRUE,...){
            
    if(save.db) {
      
      # sort regions
      regions <- sortSeqlevels(regions)
      regions <- sort(regions,ignore.strand=TRUE)
      
      getCounts <- function(data,regions,cov.bases,strand.aware){
        .setMethylDBNames(data)
        # overlap data with regions
        # convert data to GRanges without metacolumns
        g.meth=with(data, GRanges(chr, IRanges(start=start, end=end),
                                  strand =strand))
        if(!strand.aware){
          strand(g.meth)="*"
          mat=IRanges::as.matrix( findOverlaps(regions,g.meth ) )
          #mat=matchMatrix( findOverlaps(regions,g.meth ) )
          
        }else{
          mat=IRanges::as.matrix( findOverlaps(regions,g.meth) )
          #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
          
        }
        
        #require(data.table)
        # create a temporary data.table row ids from regions and counts 
        #from object
        temp.dt=data.table(id = mat[, 1], data[mat[, 2], 5:ncol(data)])
        
        coverage=NULL
        numCs=NULL
        numTs=NULL
        id=NULL
        .SD=NULL
        numTs1=covered=NULL
        # use data.table to sum up counts per region
        sum.dt=temp.dt[,c(lapply(.SD,sum),covered=length(numTs1)),by=id] 
        sum.dt=sum.dt[covered>=cov.bases,]
        temp.df=as.data.frame(regions) # get regions to a dataframe
        
        #create a new methylRaw object to return
        new.data=data.frame(#id      =new.ids,
          chr     =temp.df[sum.dt$id,"seqnames"],
          start   =temp.df[sum.dt$id,"start"],
          end     =temp.df[sum.dt$id,"end"],
          strand  =temp.df[sum.dt$id,"strand"],
          as.data.frame(sum.dt[,c(2:(ncol(sum.dt)-1)),with=FALSE]),
          stringsAsFactors=FALSE)
      }
      
      
      # catch additional args 
      args <- list(...)
      dir <- dirname(object@dbpath)
      
      if( "dbdir" %in% names(args) ){
        dir <- .check.dbdir(args$dbdir) 
      } 
      if(!( "suffix" %in% names(args) ) ){
        suffix <- paste0("_","regions")
      } else { 
        suffix <- paste0("_",args$suffix)
      }
      
      filename <- paste0(paste(object@sample.ids,collapse = "_"),suffix,".txt")
      
      #filename <- paste(basename(gsub(".txt.bgz","",object@dbpath)),suffix,".txt")
      
      newdbpath <- applyTbxByOverlap(object@dbpath,chunk.size = chunk.size, 
                                     ranges=regions,  dir=dir,
                                     filename = filename, 
                                     return.type = "tabix", 
                                     FUN = getCounts,regions,cov.bases,
                                     strand.aware)
      
      if(strand.aware & !(object@destranded) ){destranded=FALSE}else{
        destranded=TRUE}
      readMethylBaseDB(dbpath = newdbpath,dbtype = object@dbtype,
                       sample.ids=object@sample.ids,
                       assembly=object@assembly,
                       context=object@context,treatment=object@treatment,
                       destranded=destranded,resolution="region")
      
    } else {
      
      tmp <- object[]
      regionCounts(tmp,regions,cov.bases,strand.aware,save.db=FALSE)
      
    }
  }
)

#' @rdname regionCounts
#' @aliases regionCounts,methylBaseDB,GRangesList-method
setMethod("regionCounts", signature(object="methylBaseDB",regions="GRangesList"),
          function(object,regions,cov.bases,strand.aware,chunk.size,
                   save.db=TRUE,...){
            
            # combine and sort GRanges from List
            regions <- unlist(regions)
            regions <- sortSeqlevels(regions)
            regions <- sort(regions,ignore.strand=TRUE)
            regions <- unique(regions)
            
            if(save.db) {
              
              getCounts <- function(data,regions,cov.bases,strand.aware){
                .setMethylDBNames(data)
                # overlap data with regions
                # convert data to GRanges without metacolumns
                g.meth=with(data, GRanges(chr, IRanges(start=start, end=end),
                                          strand =strand))
                if(!strand.aware){
                  strand(g.meth)="*"
                  mat=IRanges::as.matrix( findOverlaps(regions,g.meth ) )
                  #mat=matchMatrix( findOverlaps(regions,g.meth ) )
                  
                }else{
                  mat=IRanges::as.matrix( findOverlaps(regions,g.meth) )
                  #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
                  
                }
                
                #require(data.table)
                # create a temporary data.table row ids from regions and counts 
                #from object
                temp.dt=data.table(id = mat[, 1], data[mat[, 2], 5:ncol(data)])
                #dt=data.table::data.table(dt)
                #dt=data.table(id=mat[,1],data[mat[,2],c(5,6,7)] ) 
                #worked with data.table 1.7.7
                
                coverage=NULL
                numCs=NULL
                numTs=NULL
                id=NULL
                .SD=NULL
                numTs1=covered=NULL
                # use data.table to sum up counts per region
                sum.dt=temp.dt[,c(lapply(.SD,sum),covered=length(numTs1)),by=id] 
                sum.dt=sum.dt[covered>=cov.bases,]
                temp.df=as.data.frame(regions) # get regions to a dataframe
                
                #create a new methylRaw object to return
                new.data=data.frame(#id      =new.ids,
                  chr     =temp.df[sum.dt$id,"seqnames"],
                  start   =temp.df[sum.dt$id,"start"],
                  end     =temp.df[sum.dt$id,"end"],
                  strand  =temp.df[sum.dt$id,"strand"],
                  as.data.frame(sum.dt[,c(2:(ncol(sum.dt)-1)),with=FALSE]),
                  stringsAsFactors=FALSE)
              }
              
              # catch additional args 
              args <- list(...)
              dir <- dirname(object@dbpath)
              
              if( "dbdir" %in% names(args) ){
                dir <- .check.dbdir(args$dbdir) 
              } 
              if(!( "suffix" %in% names(args) ) ){
                suffix <- paste0("_","regions")
              } else { 
                suffix <- paste0("_",args$suffix)
              }
              
              filename <- paste0(paste(object@sample.ids,collapse = "_"),suffix,".txt")
              
              #filename <- paste(basename(gsub(".txt.bgz","",object@dbpath)),suffix,".txt")
              
              newdbpath <- applyTbxByOverlap(object@dbpath,chunk.size = chunk.size, ranges=regions,  dir=dir,filename = filename, 
                                             return.type = "tabix", FUN = getCounts,regions,cov.bases,strand.aware)
              
              if(strand.aware & !(object@destranded) ){destranded=FALSE}else{destranded=TRUE}
              readMethylBaseDB(dbpath = newdbpath,dbtype = object@dbtype,sample.ids=object@sample.ids,
                               assembly=object@assembly,context=object@context,treatment=object@treatment,
                               destranded=destranded,resolution="region")
              
            } else {
              
              tmp <- object[]
              regionCounts(tmp,regions,cov.bases,strand.aware,save.db=FALSE)
              
            }
          }
)

#' @aliases tileMethylCounts,methylRawDB-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(object="methylRawDB"),
          function(object,win.size,step.size,cov.bases,mc.cores,
                   save.db=TRUE,...){
            
            if(save.db) {
              
              tileCount <- function(data,win.size,step.size,cov.bases,resolution) 
                {
                
                .setMethylDBNames(data,"methylRawDB")
                # overlap data with regions
                # convert data to GRanges without metacolumns
                g.meth=with(data, GRanges(chr, IRanges(start=start, end=end),
                                          strand =strand))
                
                chrs   =as.character(unique(seqnames(g.meth)))
                
                # get max length of feature covered chromosome
                max.length=max(IRanges::end(g.meth)) 
                # get start of sliding window
                min.length=min(IRanges::end(g.meth))
                win.start = floor(min.length/win.size)*win.size
                
                #get sliding windows with covered CpGs
                numTiles=floor(  (max.length-(win.size-step.size)- win.start)/step.size )+1
                
                all.wins=GRanges(seqnames=rep(chrs,numTiles),
                                 ranges=IRanges(start=win.start + 1 + 0:(numTiles-1)*step.size,
                                                width=rep(win.size,numTiles)) )
                # clean up
                rm(g.meth)
                
                tmp <- new("methylRaw",data,resolution=object@resolution)
                # clean up
                rm(data)
                
                res <- regionCounts(tmp,all.wins,cov.bases,strand.aware=FALSE)
                # clean up
                rm(tmp)
                
                return(getData(res))
                
              }
              
              # catch additional args 
              args <- list(...)
              dir <- dirname(object@dbpath)
              
              if( "dbdir" %in% names(args) ){
                if( !(is.null(args$dbdir)) ){
                  dir <- .check.dbdir(args$dbdir) 
                }
              } 
              if(!( "suffix" %in% names(args) ) ){
                suffix <- paste0("_","tiled")
              } else { 
                suffix <- paste0("_",args$suffix)
              }
              
              
              filename <- paste0(paste(object@sample.id,collapse = "_"),
                                 suffix,".txt")
              
              newdbpath <- applyTbxByChr(object@dbpath, return.type = "tabix",
                                         dir = dir,filename = filename,
                                         FUN = tileCount, win.size = win.size,
                                         step.size = step.size,
                                         cov.bases = cov.bases,
                                         resolution=object@resolution,
                                         mc.cores = mc.cores)
              
              readMethylRawDB(dbpath = newdbpath,sample.id=object@sample.id,
                              assembly=object@assembly, context =object@context,
                              resolution="region",
                              dbtype = object@dbtype)
            } else {
              
              tmp <- object[]
              tileMethylCounts(tmp,win.size,step.size,cov.bases,mc.cores)
              
            }
            
          }
)

#' @aliases tileMethylCounts,methylRawListDB-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(object="methylRawListDB"),
          function(object,win.size,step.size,cov.bases,mc.cores,save.db=TRUE,...){
            
            
            if(save.db){
              args <- list(...)
              
              if( !( "dbdir" %in% names(args)) ){
                dbdir <- NULL
              } else { dbdir <- basename(.check.dbdir(args$dbdir)) }
              
              new.list=lapply(object,tileMethylCounts,win.size,step.size,cov.bases,mc.cores,save.db,dbdir=dbdir,...) 
              new("methylRawListDB", new.list,treatment=object@treatment)
              
            } else {
              new.list=lapply(object,tileMethylCounts,win.size,step.size,cov.bases,save.db=FALSE) 
              new("methylRawList", new.list,treatment=object@treatment)
            }
            
          })

#' @aliases tileMethylCounts,methylBaseDB-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(object="methylBaseDB"),
          function(object,win.size,step.size,cov.bases,mc.cores,save.db=TRUE,...){
            
            if(save.db) {
              
              tileCount <- function(data,win.size,step.size,cov.bases,resolution,destranded) {
                
                .setMethylDBNames(data,"methylBaseDB")
                # overlap data with regions
                # convert data to GRanges without metacolumns
                g.meth=with(data, GRanges(chr, IRanges(start=start, end=end),strand =strand))
                
                chrs   =as.character(unique(seqnames(g.meth)))
                
                # get max length of feature covered chromosome
                max.length=max(IRanges::end(g.meth)) 
                # get start of sliding window
                min.length=min(IRanges::end(g.meth))
                win.start = floor(min.length/win.size)*win.size
                
                #get sliding windows with covered CpGs
                numTiles=floor(  (max.length-(win.size-step.size)- win.start)/step.size )+1
                
                all.wins=GRanges(seqnames=rep(chrs,numTiles),
                                 ranges=IRanges(start=win.start + 1 + 0:(numTiles-1)*step.size,
                                                width=rep(win.size,numTiles)) )
                
                # clean up
                rm(g.meth)
                
                tmp <- new("methylBase",data,
                           resolution=object@resolution,
                           destranded=object@destranded)
                # clean up
                rm(data)
                
                res <- regionCounts(tmp,all.wins,cov.bases,strand.aware=FALSE)
                # clean up
                rm(tmp)
                
                getData(res)
                
              }
              
              # catch additional args 
              args <- list(...)
              dir <- dirname(object@dbpath)
              
              if( "dbdir" %in% names(args) ){
                  dir <- .check.dbdir(args$dbdir) 
              } 
              if(!( "suffix" %in% names(args) ) ){
                suffix <- paste0("_","tiled")
              } else { 
                suffix <- paste0("_",args$suffix)
              }
              
              filename <- paste0(paste(object@sample.ids,collapse = "_"),suffix,".txt")
              
              newdbpath <- applyTbxByChr(object@dbpath, return.type = "tabix",dir = dir,filename = filename,
                                         FUN = tileCount, win.size = win.size,step.size = step.size,cov.bases = cov.bases,resolution=object@resolution,
                                         mc.cores = mc.cores)
              
              readMethylBaseDB(dbpath = newdbpath,dbtype = object@dbtype,sample.ids=object@sample.ids,
                               assembly=object@assembly,context=object@context,treatment=object@treatment,
                               destranded=TRUE,resolution="region")
            } else {
              
              tmp <- object[]
              tileMethylCounts(tmp,win.size,step.size,cov.bases,mc.cores)
              
            }
            
          }
)

# BedGraph methods --------------------------------------------------------

#' @rdname bedgraph-methods
#' @aliases bedgraph,methylRawDB-method
setMethod("bedgraph", signature(methylObj="methylRawDB"),
          function(methylObj,file.name,col.name,unmeth,log.transform,negative
                   ,add.on,chunk.size){
  if(!col.name %in%  c('coverage', 'numCs','numTs','perc.meth') ){
    stop("col.name argument is not one of 'coverage', 'numCs','numTs','perc.meth'")
  }
  
  bedgr <- function(data,col.name,file.name,unmeth,log.transform,
                    negative,add.on,sample.id){
    
    data <- as.data.table(data)
    .setMethylDBNames(df = data,methylDBclass = "methylRawDB")
    chr=start=end=score=numCs=coverage=.=NULL
    if(col.name=="perc.meth"){
      df= data[,.(chr,start=start-1,end,score=100*numCs/coverage )]
    }else{
      df=data[,.(chr,start=start-1,end,get(col.name) )]
      df <- as.data.frame(df)
      names(df)[4] <- col.name
      if(log.transform){
        df[,4]=log10(df[,4])
      }
      if(negative){
        df[,4]=-1*(df[,4])
      }
    }
    
    if(is.null(file.name)){
      return(df)
    }else{
      
      if(unmeth & col.name=="perc.meth")
      {
        # write meth data to single file
        write.table(df,file=paste0(file.name,"_meth"),quote=FALSE,
                    col.names=FALSE,row.names=FALSE,sep="\t",
                    append=TRUE)
        
        # write unmeth data to single file
        df[,4]=100-df[,4]
        write.table(df,file=paste0(file.name,"_unmeth"),quote=FALSE,
                    col.names=FALSE,row.names=FALSE,sep="\t",
                    append=TRUE)
        
      }else{
        write.table(df,file=file.name,quote=FALSE,col.names=FALSE,
                    row.names=FALSE,sep="\t",append=TRUE)
      }
    }
  }
  
  
  if(is.null(file.name)){
    
    applyTbxByChunk(methylObj@dbpath,chunk.size = chunk.size, 
                    return.type = "data.frame", FUN = bedgr,
                    col.name = col.name, file.name = file.name, 
                    unmeth = unmeth, log.transform=log.transform,
                    negative= negative,
                    add.on = add.on, sample.id = methylObj@sample.id)
    
  } else {
    
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(file.name,rndFile,sep="_")
    
    applyTbxByChunk(methylObj@dbpath,chunk.size = chunk.size,
                    return.type = "data.frame", FUN = bedgr,
                    col.name = col.name, file.name= filename2, 
                    unmeth = unmeth, log.transform=log.transform, 
                    negative= negative,
                    add.on = add.on, sample.id = methylObj@sample.id)
    
    
    if(unmeth & col.name=="perc.meth")
    {
      # combine single files produced by bedgr function
      track.line=paste(
        "track type=bedGraph name='",methylObj@sample.id," METH Cs",
        "' description='",methylObj@sample.id," METH Cs",
        "' visibility=full color=255,0,0 maxHeightPixels=80:80:11 ",
        add.on,sep="")                        
      cat(track.line,"\n",file=file.name)
      file.append(file.name,paste0(filename2,"_meth"))
      track.line2=paste(
        "track type=bedGraph name='",methylObj@sample.id," UNMETH Cs",
        "' description='",methylObj@sample.id," UNMETH Cs",
        "' visibility=full color=0,0,255 maxHeightPixels=80:80:11 ",
        add.on,sep="")
      cat(track.line2,"\n",file=file.name,append=TRUE)
      file.append(file.name,paste0(filename2,"_unmeth"))
      
    }else{
      
      track.line=paste(
        "track type=bedGraph name='",methylObj@sample.id," ",col.name,
        "' description='",methylObj@sample.id," ",col.name,
        "' visibility=full color=255,0,0 maxHeightPixels=80:80:11 ",
        add.on,sep="")
      cat(track.line,"\n",file=file.name)
      file.append(file.name,filename2)
      
    }
    # tidy up
    unlink(list.files(path = dirname(file.name),pattern = rndFile,
                      full.names = T))
    
  }
  
}

)

#' @rdname bedgraph-methods
#' @aliases bedgraph,methylRawListDB-method
setMethod("bedgraph", signature(methylObj="methylRawListDB"),
          function(methylObj,file.name,col.name,unmeth,log.transform,negative,add.on,chunk.size){
  if(!col.name %in%  c('coverage', 'numCs','numTs','perc.meth') ){
    stop("col.name argument is not one of 'coverage', 'numCs','numTs','perc.meth' options")
  }
  if( is.null(file.name) ){
    
    result=list()
    for(i in 1:length(methylObj))
    {
      result[[ methylObj[[i]]@sample.id ]]=bedgraph(methylObj[[i]],
                                                    file.name=NULL,
                                                    col.name=col.name,
                                                    unmeth=FALSE,
                                                    log.transform=log.transform,
                                                    negative=negative,
                                                    chunk.size=chunk.size)
    }
    return(result)
  }else{
    
    bedgraph(methylObj[[1]],file.name=file.name,
             col.name=col.name,unmeth=unmeth,log.transform=log.transform,
             negative=negative,chunk.size=chunk.size)
    
    for(i in 2:length(methylObj))
    {
      bedgraph(methylObj[[i]],file.name=paste(file.name,i,sep = "_"),
               col.name=col.name,unmeth=unmeth,log.transform=log.transform,
               negative=negative,chunk.size=chunk.size)
      
      file.append(file.name,paste(file.name,i,sep = "_"))
    }
    
    unlink(list.files(path = dirname(file.name),
                      pattern = paste0(basename(file.name),
                                       "_[[:digit:]]+"),full.names = TRUE))
    
  }
})

#' @rdname bedgraph-methods
#' @aliases bedgraph,methylDiffDB-method
setMethod("bedgraph", signature(methylObj="methylDiffDB"),
          function(methylObj,file.name,col.name,log.transform,negative,add.on,chunk.size){
            
            if(! col.name %in% c('pvalue','qvalue', 'meth.diff') )
            {
              stop("col.name argument is not one of 'pvalue','qvalue', 'meth.diff'")
            }
            
            bedgr <- function(data,col.name,file.name,log.transform,negative,add.on,sample.id){
              
              data <- as.data.table(data)
              .setMethylDBNames(df = data,methylDBclass = "methylDiffDB")
              chr=start=end=.=NULL
              df=data[,.(chr,start=start-1,end,get(col.name) )]
              df <- as.data.frame(df)
              names(df)[4] <- col.name
              if(log.transform){
                df[,4]=log10(df[,4])
              }
              if(negative){
                df[,4]=-1*(df[,4])
              }
              if(is.null(file.name)){
                return(df)
              }else{
                write.table(df,file=file.name,quote=FALSE,col.names=FALSE,
                            row.names=FALSE,sep="\t",append=TRUE)
              }
            }
            
            if(is.null(file.name)){
              
              applyTbxByChunk(methylObj@dbpath,chunk.size = chunk.size, 
                              return.type = "data.frame", FUN = bedgr,
                              col.name = col.name, file.name= file.name,
                              log.transform=log.transform, negative= negative,
                              add.on = add.on, sample.id = methylObj@sample.id)
            } else {
              
              rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUERUE),collapse="")
              filename2=paste(file.name,rndFile,sep="_")
              
              txtPath <- applyTbxByChunk(methylObj@dbpath,
                                         chunk.size = chunk.size, 
                                         return.type = "data.frame", 
                                         FUN = bedgr,
                                         col.name = col.name,
                                         file.name= filename2, 
                                         log.transform=log.transform, 
                                         negative= negative,
                                         add.on = add.on, 
                                         sample.id = methylObj@sample.id)
              
              track.line=paste(
                "track type=bedGraph name='",file.name,"' description='",col.name,
                "' visibility=full color=255,0,255 altColor=102,205,170 maxHeightPixels=80:80:11 ",add.on,sep="")
              cat(track.line,"\n",file=file.name)
              file.append(file.name,filename2)
              
              unlink(filename2)
              
            }
          }
)


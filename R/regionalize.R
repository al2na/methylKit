# functions to get regional counts

# get counts for a given regions


#
# S4 FUNCTIONS
#


#' Get regional counts for given GRanges or GRangesList object
#'
#' Convert \code{\link{methylRaw}}, \code{\link{methylRawDB}},
#'  \code{\link{methylRawList}}, 
#' \code{\link{methylRawListDB}}, \code{\link{methylBase}} or 
#' \code{\link{methylBaseDB}}  object into 
#' regional counts for a given \code{\link{GRanges}} or \code{\link{GRangesList}} object.
#' @param object a \code{\link{methylRaw}}, \code{\link{methylRawDB}}, 
#' \code{\link{methylRawList}}, 
#' \code{\link{methylRawListDB}}, \code{\link{methylBase}} or
#'  \code{\link{methylBaseDB}} object
#' @param regions a GRanges or GRangesList object. Make sure that the GRanges
#'  objects are
#'        unique in chr,start,end and strand columns.You can make them unique by 
#'        using unique() function.
#' @param cov.bases number minimum bases covered per region (Default:0). 
#' Only regions with base coverage above this threshold are returned.
#' @param strand.aware if set to TRUE only CpGs that match the strand of 
#' the region will be summarized. (default:FALSE)
#' @param chunk.size Number of rows to be taken as a chunk for processing 
#' the \code{methylDB} objects (default: 1e6)
#' @param save.db A Logical to decide whether the resulting object should be 
#' saved as flat file database or not, default: explained in Details sections  
#' @param ... optional Arguments used when save.db is TRUE
#'            
#'            \code{suffix}
#'                  A character string to append to the name of the output flat 
#'                  file database, 
#'                  only used if save.db is true, 
#'                  default actions: append \dQuote{_regions} to current filename 
#'                  if database already exists or generate new file with 
#'                  filename \dQuote{sampleID_regions} or 
#'                  \dQuote{methylBase_filtered} dependent on input object
#'                  
#'            \code{dbdir} 
#'                  The directory where flat file database(s) should be stored, 
#'                  defaults
#'                  to getwd(), working directory for newly stored databases
#'                  and to same directory for already existing database
#'                  
#'            \code{dbtype}
#'                The type of the flat file database, currently only option is "tabix"
#'                  (only used for newly stored databases)
#' 
#' @return  a new methylRaw,methylBase or methylRawList object. If \code{strand.aware} is
#'          set to FALSE (default). Even though the resulting object will have
#'          the strand information of \code{regions} it will still contain 
#'          methylation information from both strands.
#'         
#' @usage regionCounts(object,regions,cov.bases=0,strand.aware=FALSE,chunk.size,save.db,...)
#' @examples
#' data(methylKit)
#' 
#' # get the windows of interest as a GRanges object, this can be any set 
#' # of genomic locations
#' library(GenomicRanges)
#' my.win=GRanges(seqnames="chr21",
#' ranges=IRanges(start=seq(from=9764513,by=10000,length.out=20),width=5000) )
#' 
#' # getting counts per region
#' regional.methylRaw=regionCounts(object=methylRawList.obj, regions=my.win, 
#' cov.bases=0,strand.aware=FALSE)
#' 
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with 
#' \code{methylRawDB}, \code{methylBaseDB} or \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects 
#' which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most 
#' systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects as 
#' \code{methylRawDB}, \code{methylBaseDB} or \code{methylRawListDB}, 
#' while being per default FALSE for \code{methylRaw}, \code{methylBase} or
#'  \code{methylRawList}. If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the database 
#' allows the calculation in-memory, 
#' then you might want to change the value of this parameter.
#' 
#' @export
#' @docType methods
#' @rdname regionCounts
setGeneric("regionCounts", 
           function(object,regions,cov.bases=0,strand.aware=FALSE,chunk.size=1e6,save.db=FALSE,...)
           standardGeneric("regionCounts") )


# GETs regional counts for given GRanges object
# RETURNS a new methylRaw object
# @param object a \code{methylRaw} object
# @param regions a GRanges object.
#' @rdname regionCounts
#' @aliases regionCounts,methylRaw,GRanges-method
setMethod("regionCounts", signature(object="methylRaw",regions="GRanges"),
  function(object,regions,cov.bases,strand.aware,save.db=FALSE,...){
    #require(GenomicRanges)
    # sort regions
    regions <- sortSeqlevels(regions)
    regions <- sort(regions,ignore.strand=TRUE)
    # overlap object with regions
    # convert object to GRanges
    if(!strand.aware){
      g.meth=as(object,"GRanges")
      strand(g.meth)="*"
      mat=IRanges::as.matrix( findOverlaps(regions,g.meth ) )
      #mat=matchMatrix( findOverlaps(regions,g.meth ) )
      
    }else{
      mat=IRanges::as.matrix( findOverlaps(regions,as(object,"GRanges")) )
      #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
      
    }
    
    #require(data.table)
    # create a temporary data.table row ids from regions and counts from object
    coverage=numCs=numTs=id=covered=NULL
    df=data.frame(id = mat[, 1], getData(object)[mat[, 2], c(5, 6, 7)])
    dt=data.table(df)
    #dt=data.table(id=mat[,1],object[mat[,2],c(6,7,8)] ) worked with data.table 1.7.7
    

    # use data.table to sum up counts per region
    sum.dt=dt[,list(coverage=sum(coverage),
                    numCs   =sum(numCs),
                    numTs   =sum(numTs),covered=length(numTs)),by=id] 
    sum.dt=sum.dt[sum.dt$covered>=cov.bases,]
    temp.df=as.data.frame(regions) # get regions to a dataframe
                      
    # look for values with "name" in it, eg. "tx_name" or "name"
    # valuesList = names(values(regions))
    # nameid = valuesList[grep (valuesList, pattern="name")]
    
    #create id string for the new object to be returned
    #ids have to be unique and we can not assume GRanges objects will 
    #have a name attribute
    if("name" %in% names(temp.df))
    {
      new.ids=paste(temp.df[sum.dt$id,"seqnames"],temp.df[sum.dt$id,"start"],
                    temp.df[sum.dt$id,"end"],temp.df[sum.dt$id,"name"],sep=".")
      
    }else{
      new.ids=paste(temp.df[sum.dt$id,"seqnames"],temp.df[sum.dt$id,"start"],
                    temp.df[sum.dt$id,"end"],sep=".")
    }
    
    #create a new methylRaw object to return
    new.data=data.frame(#id      =new.ids,
                        chr     =temp.df[sum.dt$id,"seqnames"],
                        start   =temp.df[sum.dt$id,"start"],
                        end     =temp.df[sum.dt$id,"end"],
                        strand  =temp.df[sum.dt$id,"strand"],
                        coverage=sum.dt$coverage,
                        numCs   =sum.dt$numCs,
                        numTs   =sum.dt$numTs)
    
    if(!save.db) {
      
      new("methylRaw",new.data,sample.id=object@sample.id,
          assembly=object@assembly,context=object@context,
          resolution="region")
    
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
        suffix <- paste0("_","regions")
      } else { 
        suffix <- args$suffix
        suffix <- paste0("_",suffix)
      }
      
      # create methylRawDB
      obj <- makeMethylRawDB(df=new.data,dbpath=dbdir,dbtype="tabix",sample.id=paste0(object@sample.id,suffix),
                             assembly=object@assembly,context=object@context,resolution="region")
      obj@sample.id <- object@sample.id
      
      obj
      
    }
    
    
  }
)

#' @rdname regionCounts
#' @aliases regionCounts,methylBase,GRanges-method
setMethod("regionCounts", signature(object="methylBase",regions="GRanges"),
          function(object,regions,cov.bases,strand.aware,save.db=FALSE,...){
            #require(GenomicRanges)
            
            # sort regions
            regions <- sortSeqlevels(regions)
            regions <- sort(regions,ignore.strand=TRUE)
            
            # overlap object with regions
            # convert object to GRanges
            if(!strand.aware){
              g.meth=as(object,"GRanges")
              strand(g.meth)="*"
              mat=IRanges::as.matrix( findOverlaps(regions,g.meth ) )
              #mat=matchMatrix( findOverlaps(regions,g.meth ) )
              
            }else{
              mat=IRanges::as.matrix( findOverlaps(regions,as(object,"GRanges")) )
              #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
              
            }
            
            #require(data.table)
            # create a temporary data.table row ids from regions and counts from object
            df=data.frame(id = mat[, 1], getData(object)[mat[, 2], 5:ncol(object) ])
            dt=data.table(df)
            #dt=data.table(id=mat[,1],object[mat[,2],c(6,7,8)] ) worked with data.table 1.7.7
            
            coverage=.SD=numTs=id=numTs1=covered=NULL
            
            # use data.table to sum up counts per region
            sum.dt=dt[,c(lapply(.SD,sum, na.rm=T),covered=length(numTs1)),by=id] 
            sum.dt=sum.dt[sum.dt$covered>=cov.bases,]
            temp.df=as.data.frame(regions) # get regions to a dataframe
            
            # look for values with "name" in it, eg. "tx_name" or "name"
            # valuesList = names(values(regions))
            # nameid = valuesList[grep (valuesList, pattern="name")]
            
            #set all zero coverage tiles to missing
            for ( j in seq(2,ncol(sum.dt),by=3) ){
              data.table::set(sum.dt, which(sum.dt[[j]]==0),j:(j+2),NA)
            }


            
            #create a new methylBase object to return
            new.data=data.frame(#id      =new.ids,
              chr     =temp.df[sum.dt$id,"seqnames"],
              start   =temp.df[sum.dt$id,"start"],
              end     =temp.df[sum.dt$id,"end"],
              strand  =temp.df[sum.dt$id,"strand"],
              as.data.frame(sum.dt[,c(2:(ncol(sum.dt)-1)),with=FALSE]),stringsAsFactors=FALSE)
            
            if(strand.aware & !(object@destranded) ){destranded=FALSE}else{destranded=TRUE}
            
            if(!save.db) {
              
            new("methylBase",new.data,sample.ids=object@sample.ids,
                assembly=object@assembly,context=object@context,treatment=object@treatment,
                coverage.index=object@coverage.index,numCs.index=object@numCs.index,
                numTs.index=object@numTs.index,destranded=destranded,
                resolution="region")
           
          } else {
            
            # catch additional args 
            args <- list(...)
            
            if( !( "dbdir" %in% names(args)) ){
              dbdir <- .check.dbdir(getwd())
            } else { dbdir <- .check.dbdir(args$dbdir) }
            if(!( "suffix" %in% names(args) ) ){
              suffix <- "_regions"
            } else { 
              suffix <- paste0("_",args$suffix)
            }
            
            # create methylBaseDB
            makeMethylBaseDB(df=new.data,dbpath=dbdir,dbtype="tabix",
                             sample.ids=object@sample.ids,
                             assembly=object@assembly,context=object@context,
                             treatment=object@treatment,
                             coverage.index=object@coverage.index,
                             numCs.index=object@numCs.index,
                             numTs.index=object@numTs.index,destranded=destranded,
                             resolution="region", suffix=suffix )
          }
             
          }
)




# RETURNS a new methylRawList object
# gets regional counts for all elements in methylRawList for given regions
# MAKE SURE an element of the list, which will be a set of GRanges rows, 
#  are not on different chromosomes and strands
# Also, make sure id column of returned methylRaw object is unique
# you can add refseq id to the id column: chr.start.end.refseqid
# @param object a \code{methylRaw} object
# @param regions a GRangesList object.
#' @rdname regionCounts
#' @aliases regionCounts,methylRaw,GRangesList-method
# assume that each name of the element in the GRangesList is unique and 
setMethod("regionCounts", signature(object="methylRaw",regions="GRangesList"),
  function(object,regions,cov.bases,strand.aware,save.db=FALSE,...){
            
    #require(GenomicRanges)
    
    # combine and sort GRanges from List
    regions <- unlist(regions)
    regions <- sortSeqlevels(regions)
    regions <- sort(regions,ignore.strand=TRUE)
    regions <- unique(regions)
    
    # overlap object with regions
    # convert object to GRanges
    #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
    
    if(!strand.aware){
      g.meth=as(object,"GRanges")
      strand(g.meth)="*"
      mat=IRanges::as.matrix( findOverlaps(regions,g.meth ) )
      #mat=matchMatrix( findOverlaps(regions,g.meth ) )
      
    }else{
      mat=IRanges::as.matrix( findOverlaps(regions,as(object,"GRanges")) )
      #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
      
    }
    
    
    
    
    #require(data.table)
    # create a temporary data.table row ids from regions and counts from object
    df=data.frame(id = mat[, 1], getData(object)[mat[, 2], c(5, 6, 7)])
    dt=data.table(df)
    #dt=data.table(id=mat[,1],object[mat[,2],c(6,7,8)] ) worked with data.table 1.7.7
    
    coverage=NULL
    numCs=NULL
    numTs=NULL
    id=NULL
    # use data.table to sum up counts per region
    sum.dt=dt[,list(coverage=sum(coverage),
                    numCs   =sum(numCs),
                    numTs   =sum(numTs),covered=length(numTs)),by=id] 
    sum.dt=sum.dt[sum.dt$covered>=cov.bases,]
    temp.df=as.data.frame(regions) # get regions to a dataframe
    
    # look for values with "name" in it, eg. "tx_name" or "name"
    # valuesList = names(values(regions))
    # nameid = valuesList[grep (valuesList, pattern="name")]
    
    #create id string for the new object to be returned
    #ids have to be unique and we can not assume GRanges objects will 
    #have a name attribute
    if("name" %in% names(temp.df))
    {
      new.ids=paste(temp.df[sum.dt$id,"seqnames"],temp.df[sum.dt$id,"start"],
                    temp.df[sum.dt$id,"end"],temp.df[sum.dt$id,"name"],sep=".")
      
    }else{
      new.ids=paste(temp.df[sum.dt$id,"seqnames"],temp.df[sum.dt$id,"start"],
                    temp.df[sum.dt$id,"end"],sep=".")
    }
    
    #create a new methylRaw object to return
    new.data=data.frame(#id      =new.ids,
      chr     =temp.df[sum.dt$id,"seqnames"],
      start   =temp.df[sum.dt$id,"start"],
      end     =temp.df[sum.dt$id,"end"],
      strand  =temp.df[sum.dt$id,"strand"],
      coverage=sum.dt$coverage,
      numCs   =sum.dt$numCs,
      numTs   =sum.dt$numTs)
    
    
    if(!save.db) {
      new("methylRaw",new.data,sample.id=object@sample.id,
          assembly=object@assembly,context=object@context,
          resolution="region")   
    
    
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
        suffix <- paste0("_","regions")
      } else { 
        suffix <- args$suffix
        suffix <- paste0("_",suffix)
      }
      
      # create methylRawDB
      obj <- makeMethylRawDB(df=new.data,dbpath=dbdir,dbtype="tabix",sample.id=paste0(object@sample.id,suffix),
                             assembly=object@assembly,context=object@context,resolution="region")
      obj@sample.id <- object@sample.id
      
      obj
      
    }
  }
)
# Note: some genes do not have intron, need to take care of it.


#' @rdname regionCounts
#' @aliases regionCounts,methylBase,GRangesList-method
setMethod("regionCounts", signature(object="methylBase",regions="GRangesList"),
          function(object,regions,cov.bases,strand.aware,save.db=FALSE,...){
            #require(GenomicRanges)
            
            # combine and sort GRanges from List
            regions <- unlist(regions)
            regions <- sortSeqlevels(regions)
            regions <- sort(regions,ignore.strand=TRUE)
            regions <- unique(regions)
            
            # overlap object with regions
            # convert object to GRanges
            if(!strand.aware){
              g.meth=as(object,"GRanges")
              strand(g.meth)="*"
              mat=IRanges::as.matrix( findOverlaps(regions,g.meth ) )
              #mat=matchMatrix( findOverlaps(regions,g.meth ) )
              
            }else{
              mat=IRanges::as.matrix( findOverlaps(regions,as(object,"GRanges")) )
              #mat=matchMatrix( findOverlaps(regions,as(object,"GRanges")) )
              
            }
            
            #require(data.table)
            # create a temporary data.table row ids from regions and counts from object
            df=data.frame(id = mat[, 1], getData(object)[mat[, 2], 5:ncol(object) ])
            dt=data.table(df)
            #dt=data.table(id=mat[,1],object[mat[,2],c(6,7,8)] ) worked with data.table 1.7.7
            
            id=.SD=numTs1=NULL
            # use data.table to sum up counts per region
            # treat missing values as they had zero coverage
            sum.dt=dt[,c(lapply(.SD,sum,na.rm=T),covered=length(numTs1)),by=id] 
            sum.dt=sum.dt[sum.dt$covered>=cov.bases,]
            temp.df=as.data.frame(regions) # get regions to a dataframe
            
            # look for values with "name" in it, eg. "tx_name" or "name"
            # valuesList = names(values(regions))
            # nameid = valuesList[grep (valuesList, pattern="name")]
            
            #set all zero coverage tiles to missing
            for ( j in seq(2,ncol(sum.dt),by=3) ){
              data.table::set(sum.dt, which(sum.dt[[j]]==0),j:(j+2),NA)
            }
            
            #create a new methylBase object to return
            new.data=data.frame(#id      =new.ids,
              chr     =temp.df[sum.dt$id,"seqnames"],
              start   =temp.df[sum.dt$id,"start"],
              end     =temp.df[sum.dt$id,"end"],
              strand  =temp.df[sum.dt$id,"strand"],
              as.data.frame(sum.dt[,c(2:(ncol(sum.dt)-1)),with=FALSE]),stringsAsFactors=FALSE)
            
            if(strand.aware & !(object@destranded) ){destranded=FALSE}else{destranded=TRUE}
            
            if(!save.db) {
              
              new("methylBase",new.data,sample.ids=object@sample.ids,
                  assembly=object@assembly,context=object@context,treatment=object@treatment,
                  coverage.index=object@coverage.index,numCs.index=object@numCs.index,
                  numTs.index=object@numTs.index,destranded=destranded,
                  resolution="region")
              
            } else {
              
              # catch additional args 
              args <- list(...)
              
              if( !( "dbdir" %in% names(args)) ){
                dbdir <- .check.dbdir(getwd())
              } else { dbdir <- .check.dbdir(args$dbdir) }
              if(!( "suffix" %in% names(args) ) ){
                suffix <- "_regions"
              } else { 
                suffix <- paste0("_",args$suffix)
              }
              
              # create methylBaseDB
              makeMethylBaseDB(df=new.data,dbpath=dbdir,dbtype="tabix",sample.ids=object@sample.ids,
                               assembly=object@assembly,context=object@context,treatment=object@treatment,
                               coverage.index=object@coverage.index,numCs.index=object@numCs.index,
                               numTs.index=object@numTs.index,destranded=destranded,
                               resolution="region", suffix=suffix )
            }
            
          }
)

# GETs regional counts for given GRanges object
# RETURNS a new methylRawList object
# @param object a \code{methylRawList} object
# @param regions a GRanges object.
#' @rdname regionCounts
#' @aliases regionCounts,methylRawList,GRanges-method
setMethod("regionCounts", signature(object="methylRawList",regions="GRanges"),
  function(object,regions,cov.bases,strand.aware,save.db=FALSE,...){
    
    if(!save.db) {
      outList=list()
      for(i in 1:length(object))
      {
        obj = regionCounts(object = object[[i]],
                           regions=regions,
                           cov.bases,strand.aware)
        outList[[i]] = obj
      }
      myobj=new("methylRawList", outList,treatment=object@treatment)
      myobj
      
    } else {
     
      args <- list(...)
      if( !( "dbdir" %in% names(args)) ){
        dbdir <- .check.dbdir(getwd())
      } else { dbdir <- .check.dbdir(args$dbdir) }
    
          outList = lapply(object,regionCounts, regions=regions,
                                 cov.bases,strand.aware,save.db = TRUE, dbdir = basename(dbdir) ,...)
      
        new("methylRawListDB", outList,treatment=object@treatment)
    }
  }
)

# GETs regional counts for given GRangesList object
# RETURNS a new methylRawList object
# @param object a \code{methylRawList} object
# @param regions a GRangesList object.
#' @rdname regionCounts
#' @aliases regionCounts,methylRawList,GRangesList-method
setMethod("regionCounts", signature(object="methylRawList",
                                    regions="GRangesList"),
  function(object,regions,cov.bases,strand.aware,save.db=FALSE,...){
    
    if(!save.db) {
      outList=list()
      for(i in 1:length(object))
      {
        obj = regionCounts(object = object[[i]],
                           regions=regions,
                           cov.bases,strand.aware)
        outList[[i]] = obj
      }
      myobj=new("methylRawList", outList,treatment=object@treatment)
      myobj
      
    } else {
      
      args <- list(...)
      if( !( "dbdir" %in% names(args)) ){
        dbdir <- .check.dbdir(getwd())
      } else { dbdir <- .check.dbdir(args$dbdir) }
      
      outList = lapply(object,regionCounts, regions=regions,
                       cov.bases,strand.aware,save.db = TRUE, dbdir = basename(dbdir) ,...)
      
      new("methylRawListDB", outList,treatment=object@treatment)
    }
  }
)


#' Get methylated/unmethylated base counts for tilling windows 
#'
#' The function summarizes methylated/unmethylated base counts over tilling 
#' windows accross genome. This function can be used when differential 
#' methylated analysis is preferable to tilling windows instead of base pairs.
#'
#' @param object \code{\link{methylRaw}}, \code{\link{methylRawDB}},
#'   \code{\link{methylRawList}}, \code{\link{methylRawListDB}},
#'   \code{\link{methylBase}} or \code{\link{methylBaseDB}} object containing
#'   base pair resolution methylation information
#' @param win.size an integer for the size of the tiling windows
#' @param step.size an integer for the step size of tiling windows
#' @param cov.bases minimum number of bases to be covered in a given window
#' @param mc.cores number of cores to use when processing \code{methylDB}
#'   objects, default: 1, but always 1 for Windows)
#' @param save.db A Logical to decide whether the resulting object should be 
#'   saved as flat file database or not, default: explained in Details sections
#' @param ... optional Arguments used when save.db is TRUE
#'            
#'            \code{suffix}
#'                  A character string to append to the name of the output flat 
#'                  file database, 
#'                  only used if save.db is true, 
#'                  default actions: append \dQuote{_tiled} to current filename 
#'                  if database already exists or generate new file with 
#'                  filename \dQuote{sampleID_tiled} or 
#'                  \dQuote{methylBase_tiled} dependent on input object
#'                  
#'            \code{dbdir} 
#'                  The directory where flat file database(s) should be stored,
#'                  defaults to getwd(), working directory for newly stored
#'                  databases and to same directory for already existing
#'                  database
#'                  
#            \code{dbtype}
#                  The type of the flat file database, currently only option is
#                  "tabix" (only used for newly stored databases)
#'
#' @usage tileMethylCounts(object,win.size=1000,step.size=1000,cov.bases=0,mc.cores=1,save.db,...)
#' @return \code{methylRaw},\code{methylBase} or \code{methylRawList} object
#' @export
#' @examples
#' data(methylKit)
#' 
#' tiled.methylRaw=tileMethylCounts(object=methylRawList.obj,win.size=1000,
#'                                  step.size=1000,cov.bases=0)
#' 
#' 
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with 
#' \code{methylRawDB}, \code{methylBaseDB} or \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects 
#' which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most 
#' systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects as 
#' \code{methylRawDB}, \code{methylBaseDB} or \code{methylRawListDB}, 
#' while being per default FALSE for \code{methylRaw}, \code{methylBase} or
#'  \code{methylRawList}. If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the database 
#' allows the calculation in-memory, 
#' then you might want to change the value of this parameter.
#' 
#' @docType methods
#' @rdname tileMethylCounts-methods
setGeneric("tileMethylCounts", 
           function(object,win.size=1000,step.size=1000,cov.bases=0,mc.cores=1,save.db=FALSE,...)
             standardGeneric("tileMethylCounts") )

#' @aliases tileMethylCounts,methylRaw-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(object="methylRaw"),
  function(object,win.size,step.size,cov.bases,save.db=FALSE,...){
    
    g.meth =as(object,"GRanges")
    #chrs   =IRanges::levels(seqnames(g.meth))
    chrs   =as.character(unique(seqnames(g.meth)))
    #widths =seqlengths(g.meth) # this doesn't work with BioC 3.0
    widths =sapply(chrs,function(x,y) max(end(y[seqnames(y)==x,])),g.meth  )# lengths of max bp in each chr
    all.wins=GRanges()
    for(i in 1:length(chrs))
    {
      # get max length of feature covered chromosome
      max.length=max(IRanges::end(g.meth[seqnames(g.meth)==chrs[i],])) 
      
      #get sliding windows with covered CpGs
      numTiles=floor(  (max.length-(win.size-step.size) )/step.size )+1
      numTiles=ifelse(numTiles<1, 1,numTiles)
      temp.wins=GRanges(seqnames=rep(chrs[i],numTiles),
                        ranges=IRanges(start=1+0:(numTiles-1)*step.size,
                                       width=rep(win.size,numTiles)) )
      all.wins=suppressWarnings(c(all.wins,temp.wins))
    }
    #catch additional args
    args <- list(...)
    if( !( "suffix" %in% names(args)) ){
      suffix <- "tiled"
    } else { suffix <- args$suffix }
    
    regionCounts(object,all.wins,cov.bases,strand.aware=FALSE,save.db=save.db,suffix=suffix,... = ...)
  }
)

#' @aliases tileMethylCounts,methylRawList-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(object="methylRawList"),
  function(object,win.size,step.size,cov.bases,save.db=FALSE,...){
    
    if (save.db) {
    
    #catch additional args
    args <- list(...)
    if( !( "dbdir" %in% names(args)) ){
      dbdir <- .check.dbdir(getwd())
    } else { dbdir <- .check.dbdir(args$dbdir) }
    
    new.list=lapply(object,tileMethylCounts,win.size,step.size,cov.bases,save.db=save.db,dbdir=basename(dbdir),... = ...) 
    new("methylRawListDB", new.list,treatment=object@treatment)
    
    } else {
      
      new.list=lapply(object,tileMethylCounts,win.size,step.size,cov.bases) 
      new("methylRawList", new.list,treatment=object@treatment)
      
    }
    
    
    
})



#' @aliases tileMethylCounts,methylBase-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(object="methylBase"),
          function(object,win.size,step.size,cov.bases,save.db=FALSE,...){
            
            g.meth =as(object,"GRanges")
            #chrs   =IRanges::levels(seqnames(g.meth))
            chrs   =as.character(unique(seqnames(g.meth)))
            
            #widths =seqlengths(g.meth)
            widths =sapply(chrs,function(x,y) max(end(y[seqnames(y)==x,])),g.meth  )# lengths of max bp in each chr
            
            all.wins=GRanges()
            for(i in 1:length(chrs))
            {
              # get max length of feature covered chromosome
              max.length=max(IRanges::end(g.meth[seqnames(g.meth)==chrs[i],])) 
              
              #get sliding windows with covered CpGs
              numTiles=floor(  (max.length-(win.size-step.size) )/step.size )+1
              numTiles=ifelse(numTiles<1, 1,numTiles)
              temp.wins=GRanges(seqnames=rep(chrs[i],numTiles),
                                ranges=IRanges(start=1+0:(numTiles-1)*step.size,
                                               width=rep(win.size,numTiles)) )
              all.wins=suppressWarnings(c(all.wins,temp.wins))
            }
            #catch additional args
            args <- list(...)
            if( !( "suffix" %in% names(args)) ){
              suffix <- "tiled"
            } else { suffix <- args$suffix }
            
            regionCounts(object,all.wins,cov.bases,strand.aware=FALSE,save.db=save.db,suffix=suffix,... = ...)
          }
)


# back up code:
# annotation of methylation states by genes
# including promoter, intron, exon, nearest cpgi and cpgi shore.


# valPerIndex=function(query, subject, col="cov"){
# cov.per=value.per.location(query,subject,col=col)
# all.ind=1:length(query)
# rm=data.frame(ind=all.ind[!all.ind %in% (cov.per[,1])],
#  val=rep(0,length(all.ind[!all.ind %in% (cov.per[,1])]) ),
#  numRows=rep(0,length(all.ind[!all.ind %in% (cov.per[,1])]) ) )
# cov.per2=rbind(cov.per,rm)
# cov.per3=data.frame(cov.per2, tx_id=names(query)[cov.per2[,1]])
# return(cov.per3)
# }
# 
# meth.ratio=function(meth.idx, cov.idx){
#     if(identical(meth.idx[,1], cov.idx[,1])){
#     ratio = meth.idx[,2]/cov.idx[,2]
#     names(ratio) = meth.idx[,4]
#     return(ratio)
#     }
#     else stop("index does not match")
# }
# 
# # using exons
# meth.idx.exons=valPerIndex(exonRanges, g.meth, col="methCs")
# cov.idx.exons=valPerIndex(exonRanges, g.meth, col="cov")
# 
# exon.meth.ratio=meth.ratio(meth.idx.exons, cov.idx.exons)


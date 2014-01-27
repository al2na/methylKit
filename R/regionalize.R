# functions to get regional counts

# get counts for a given regions


#
# S4 FUNCTIONS
#


#' Get regional counts for given GRanges or GRangesList object
#'
#' Convert \code{\link{methylRaw}}, \code{\link{methylRawList}} or \code{\link{methylBase}}  object into 
#' regional counts for a given \code{\link{GRanges}} or \code{\link{GRangesList}} object.
#' @param object a \code{methylRaw} or \code{methlRawList} object
#' @param regions a GRanges or GRangesList object. Make sure that the GRanges objects are
#'        unique in chr,start,end and strand columns.You can make them unique by 
#'        using unique() function.
#' @param cov.bases number minimum bases covered per region (Default:0). 
#' Only regions with base coverage above this threshold are returned.
#' @param strand.aware if set to TRUE only CpGs that match the strand of 
#' the region will be summarized. (default:FALSE)
#' 
#' @return  a new methylRaw,methylBase or methylRawList object. If \code{strand.aware} is
#'          set to FALSE (default). Even though the resulting object will have
#'          the strand information of \code{regions} it will still contain 
#'          methylation information from both strands.
#'         
#' @usage regionCounts(object,regions,cov.bases=0,strand.aware=FALSE)
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
#' 
#' @export
#' @docType methods
#' @rdname regionCounts
setGeneric("regionCounts", 
           function(object,regions,cov.bases=0,strand.aware=FALSE)
           standardGeneric("regionCounts") )


# GETs regional counts for given GRanges object
# RETURNS a new methylRaw object
# @param object a \code{methylRaw} object
# @param regions a GRanges object.
#' @rdname regionCounts
#' @aliases regionCounts,methylRaw,GRanges-method
setMethod("regionCounts", signature(object="methylRaw",regions="GRanges"),
  function(object,regions,cov.bases,strand.aware){
    #require(GenomicRanges)
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
    df=data.frame(id = mat[, 1], getData(object)[mat[, 2], c(5, 6, 7)])
    dt=data.table::data.table(df)
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
    
    new("methylRaw",new.data,sample.id=object@sample.id,
        assembly=object@assembly,context=object@context,
        resolution="region")
    
  }
)

#' @rdname regionCounts
#' @aliases regionCounts,methylBase,GRanges-method
setMethod("regionCounts", signature(object="methylBase",regions="GRanges"),
          function(object,regions,cov.bases,strand.aware){
            #require(GenomicRanges)
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
            dt=data.table::data.table(df)
            #dt=data.table(id=mat[,1],object[mat[,2],c(6,7,8)] ) worked with data.table 1.7.7
            
            # use data.table to sum up counts per region
            sum.dt=dt[,c(lapply(.SD,sum),covered=length(numTs1)),by=id] 
            sum.dt=sum.dt[sum.dt$covered>=cov.bases,]
            temp.df=as.data.frame(regions) # get regions to a dataframe
            
            # look for values with "name" in it, eg. "tx_name" or "name"
            # valuesList = names(values(regions))
            # nameid = valuesList[grep (valuesList, pattern="name")]
            


            
            #create a new methylRaw object to return
            new.data=data.frame(#id      =new.ids,
              chr     =temp.df[sum.dt$id,"seqnames"],
              start   =temp.df[sum.dt$id,"start"],
              end     =temp.df[sum.dt$id,"end"],
              strand  =temp.df[sum.dt$id,"strand"],
              as.data.frame(sum.dt[,,c(2:(ncol(sum.dt)-1))]),stringsAsFactors=FALSE)
            
            if(strand.aware & !(object@destranded) ){destranded=FALSE}else{destranded=TRUE}
            new("methylBase",new.data,sample.ids=object@sample.ids,
                assembly=object@assembly,context=object@context,treatment=object@treatment,
                coverage.index=object@coverage.index,numCs.index=object@numCs.index,
                numTs.index=object@numTs.index,destranded=destranded,
                resolution="region")
            
          }
)

# RETURNS a new methylRawList object
# gets regional counts for all elements in methylRawList for given regions
# MAKE SURE an element of the list, which will be a set of GRanges rows, 
#  are not on different chromsomes and strands
# Also, make sure id column of returned methylRaw object is unique
# you can add refseq id to the id column: chr.start.end.refseqid
# @param object a \code{methylRaw} object
# @param regions a GRangesList object.
#' @rdname regionCounts
#' @aliases regionCounts,methylRaw,GRangesList-method
# assume that each name of the element in the GRangesList is unique and 
setMethod("regionCounts", signature(object="methylRaw",regions="GRangesList"),
  function(object,regions,cov.bases,strand.aware){
            
    #require(GenomicRanges)
    
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
    dt=data.table::data.table(id=mat[,1],getData(object)[mat[,2],c(5,6,7)] )
    
    # use data.table to sum up counts per region
    sum.dt=dt[,list(coverage=sum(coverage),
                    numCs   =sum(numCs),
                    numTs   =sum(numTs),covered=length(numTs)),by=id] 
    sum.dt=sum.dt[sum.dt$covered>=cov.bases,]
    
    #temp.df=as.data.frame(regions) # get regions to a dataframe
    temp.dt=data.table(as.data.frame(regions))
    first.element = function(x) x[1]
    sum.temp.dt=temp.dt[, list(chr = first.element(seqnames), 
                               start = min(start), 
                               end = max(end),
                               strand = first.element(strand)),
                        by=element]
    # if it is intron, some of the symbols may not have any intron GRanges
    sum.id=names(regions)[sum.dt$id]
    #create a new methylRaw object to return
    new.data=data.frame(
      #id      =sum.id,
      chr     =sum.temp.dt$chr[sum.temp.dt$element %in% sum.id],
      start   =sum.temp.dt$start[sum.temp.dt$element %in% sum.id],
      end     =sum.temp.dt$end[sum.temp.dt$element %in% sum.id],
      strand  =sum.temp.dt$strand[sum.temp.dt$element %in% sum.id],
      coverage=sum.dt$coverage,
      numCs   =sum.dt$numCs,
      numTs   =sum.dt$numTs)
    
    new("methylRaw",new.data,sample.id=object@sample.id,
        assembly=object@assembly,context=object@context,
        resolution="region")
    
  }
)
# Note: some genes do not have intron, need to take care of it.


# GETs regional counts for given GRanges object
# RETURNS a new methylRawList object
# @param object a \code{methylRawList} object
# @param regions a GRanges object.
#' @rdname regionCounts
#' @aliases regionCounts,methylRawList,GRanges-method
setMethod("regionCounts", signature(object="methylRawList",regions="GRanges"),
  function(object,regions,cov.bases,strand.aware){
    
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
  function(object,regions,cov.bases,strand.aware){
    
    outList=list()
    
    for(i in 1:length(object))
    {
      obj = regionCounts(object = object[[i]],regions=regions,
                         cov.bases,strand.aware)
      outList[[i]] = obj
    }
    
    myobj=new("methylRawList", outList,treatment=object@treatment)
    myobj
  }
)


#' Get methylated/unmethylated base counts for tilling windows 
#'
#' The function summarizes methylated/unmethylated base counts over tilling 
#' windows accross genome. This function can be used when differential 
#' methylated analysis is preferable to tilling windows instead of base pairs.
#'
#' @param object \code{\link{methylRaw}}, \code{\link{methylRawList}} or \code{\link{methylBase}} 
#' object containing base pair resolution methylation information
#' @param win.size an integer for the size of the tiling windows
#' @param step.size an integer for the step size of tiling windows
#' @param cov.bases minimum number of bases to be covered in a given window
#' @usage tileMethylCounts(object,win.size=1000,step.size=1000,cov.bases=0)
#' @return \code{methylRaw},\code{methylBase} or \code{methylRawList} object
#' @export
#' @examples
#' data(methylKit)
#' 
#' tiled.methylRaw=tileMethylCounts(object=methylRawList.obj,win.size=1000,
#'                                  step.size=1000,cov.bases=0)
#' 
#' 
#' @docType methods
#' @rdname tileMethylCounts-methods
setGeneric("tileMethylCounts", 
           function(object,win.size=1000,step.size=1000,cov.bases=0)
             standardGeneric("tileMethylCounts") )

#' @aliases tileMethylCounts,methylRaw-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(object="methylRaw"),
  function(object,win.size,step.size,cov.bases){
    
    g.meth =as(object,"GRanges")
    #chrs   =IRanges::levels(seqnames(g.meth))
    chrs   =Ias.character(unique(seqnames(g.meth)))
    widths =seqlengths(g.meth)
    all.wins=GRanges()
    for(i in 1:length(chrs))
    {
      # get max length of feature covered chromosome
      max.length=max(IRanges::end(g.meth[seqnames(g.meth)==chrs[i],])) 
      
      #get sliding windows with covered CpGs
      numTiles=floor(  (max.length-(win.size-step.size) )/step.size )+1
      temp.wins=GRanges(seqnames=rep(chrs[i],numTiles),
                        ranges=IRanges(start=1+0:(numTiles-1)*step.size,
                                       width=rep(win.size,numTiles)) )
      all.wins=suppressWarnings(c(all.wins,temp.wins))
    }
    regionCounts(object,all.wins,cov.bases,strand.aware=FALSE)
  }
)

#' @aliases tileMethylCounts,methylRawList-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(object="methylRawList"),
  function(object,win.size,step.size,cov.bases){
    
    new.list=lapply(object,tileMethylCounts,win.size,step.size,cov.bases) 
    new("methylRawList", new.list,treatment=object@treatment)
    
})



#' @aliases tileMethylCounts,methylBase-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(object="methylBase"),
          function(object,win.size,step.size,cov.bases){
            
            g.meth =as(object,"GRanges")
            chrs   =IRanges::levels(seqnames(g.meth))
            widths =seqlengths(g.meth)
            all.wins=GRanges()
            for(i in 1:length(chrs))
            {
              # get max length of feature covered chromosome
              max.length=max(IRanges::end(g.meth[seqnames(g.meth)==chrs[i],])) 
              
              #get sliding windows with covered CpGs
              numTiles=floor(  (max.length-(win.size-step.size) )/step.size )+1
              temp.wins=GRanges(seqnames=rep(chrs[i],numTiles),
                                ranges=IRanges(start=1+0:(numTiles-1)*step.size,
                                               width=rep(win.size,numTiles)) )
              all.wins=suppressWarnings(c(all.wins,temp.wins))
            }
            regionCounts(object,all.wins,cov.bases,strand.aware=FALSE)
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


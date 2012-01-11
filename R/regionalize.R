# functions to get regional counts

# get counts for a given regions


#
# S4 FUNCTIONS
#


#' GETs regional counts for given GRanges or GRangesList object
#'
#' Convert \code{methylRaw} or \code{methylRawList} object into regional counts for a given GRanges or GRangesList object.
#' @param methylObj a \code{methylRaw} or \code{methlRawList} object
#' @param regions a GRanges or GRangesList object.
#' @param cov.bases number minimum bases covered per region (Default:0). Only regions with base coverage above this threshold are returned.
#' 
#' @return RETURNS a new methylRaw or methylRawList object
#' @usage regionCounts(methylObj,regions,cov.bases=0)
#' @export
#' @docType methods
#' @rdname regionCounts-methods
setGeneric("regionCounts", function(methylObj,regions,cov.bases=0) standardGeneric("regionCounts") )


# GETs regional counts for given GRanges object
# RETURNS a new methylRaw object
# @param methylObj a \code{methylRaw} object
# @param regions a GRanges object.
#' @rdname regionCounts-methods
#' @aliases regionCounts,methylRaw,GRanges-method
setMethod("regionCounts", signature(methylObj="methylRaw",regions="GRanges"),
                    function(methylObj,regions,cov.bases){
                      #require(GenomicRanges)
                      # overlap methylObj with regions
                      # convert methylObj to GRanges
                      mat=matchMatrix( findOverlaps(regions,as(methylObj,"GRanges")) )
                      
                      # create a temporary data.table row ids from regions and counts from methylObj
                      dt=data.table(id=mat[,1],methylObj[mat[,2],c(6,7,8)] )
                      
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
                      #ids have to be unique and we can not assume GRanges objects will have a name attribute
                      if("name" %in% names(temp.df))
                      {
                        new.ids=paste(temp.df[sum.dt$id,"seqnames"],temp.df[sum.dt$id,"start"],
                                      temp.df[sum.dt$id,"end"],temp.df[sum.dt$id,"name"],sep=".")

                      }else{
                        new.ids=paste(temp.df[sum.dt$id,"seqnames"],temp.df[sum.dt$id,"start"],temp.df[sum.dt$id,"end"],sep=".")
                      }
                      
                      #create a new methylRaw object to return
                      new.data=data.frame(id      =new.ids,
                                          chr     =temp.df[sum.dt$id,"seqnames"],
                                          start   =temp.df[sum.dt$id,"start"],
                                          end     =temp.df[sum.dt$id,"end"],
                                          strand  =temp.df[sum.dt$id,"strand"],
                                          coverage=sum.dt$coverage,
                                          numCs   =sum.dt$numCs,
                                          numTs   =sum.dt$numTs)
                    
                    new("methylRaw",new.data,sample.id=methylObj@sample.id,assembly=methylObj@assembly,context=methylObj@context,resolution="region")
                      
                    }
)


# RETURNS a new methylRawList object
# gets regional counts for all elements in methylRawList for given regions
# MAKE SURE an element of the list, which will be a set of GRanges rows, are not on different chromsomes and strands
# Also, make sure id column of returned methylRaw object is unique
# you can add refseq id to the id column: chr.start.end.refseqid
# @param methylObj a \code{methylRaw} object
# @param regions a GRangesList object.
#' @rdname regionCounts-methods
#' @aliases regionCounts,methylRaw,GRangesList-method
# assume that each name of the element in the GRangesList is unique and 
setMethod("regionCounts", signature(methylObj="methylRaw",regions="GRangesList"),
                    function(methylObj,regions,cov.bases){
                      
                      require(GenomicRanges)
                      
                      # overlap methylObj with regions
                      # convert methylObj to GRanges
                      mat=matchMatrix( findOverlaps(regions,as(methylObj,"GRanges")) )
                      
                      # create a temporary data.table row ids from regions and counts from methylObj
                      dt=data.table(id=mat[,1],methylObj[mat[,2],c(6,7,8)] )
                      
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
                      new.data=data.frame(id      =sum.id,
                                          chr     =sum.temp.dt$chr[sum.temp.dt$element %in% sum.id],
                                          start   =sum.temp.dt$start[sum.temp.dt$element %in% sum.id],
                                          end     =sum.temp.dt$end[sum.temp.dt$element %in% sum.id],
                                          strand  =sum.temp.dt$strand[sum.temp.dt$element %in% sum.id],
                                          coverage=sum.dt$coverage,
                                          numCs   =sum.dt$numCs,
                                          numTs   =sum.dt$numTs)
                    
                    new("methylRaw",new.data,sample.id=methylObj@sample.id,assembly=methylObj@assembly,context=methylObj@context,resolution="region")
  
                    })
# Note: some genes do not have intron, need to take care of it.


# GETs regional counts for given GRanges object
# RETURNS a new methylRawList object
# @param methylObj a \code{methylRawList} object
# @param regions a GRanges object.
#' @rdname regionCounts-methods
#' @aliases regionCounts,methylRawList,GRanges-method
setMethod("regionCounts", signature(methylObj="methylRawList",regions="GRanges"),
                    function(methylObj,regions,cov.bases){
                     
                    outList=list()

                    for(i in 1:length(methylObj))
                    {
                        obj = regionCounts(methylObj = methylObj[[i]],regions=regions,cov.bases)
                        outList[[i]] = obj
                    }
                    
                    myobj=new("methylRawList", outList,treatment=methylObj@treatment)
                    myobj
                    }
)

# GETs regional counts for given GRangesList object
# RETURNS a new methylRawList object
# @param methylObj a \code{methylRawList} object
# @param regions a GRangesList object.
#' @rdname regionCounts-methods
#' @aliases regionCounts,methylRawList,GRangesList-method
setMethod("regionCounts", signature(methylObj="methylRawList",regions="GRangesList"),
                    function(methylObj,regions,cov.bases){
                     
                    outList=list()

                    for(i in 1:length(methylObj))
                    {
                        obj = regionCounts(methylObj = methylObj[[i]],regions=regions,cov.bases)
                        outList[[i]] = obj
                    }
                    
                    myobj=new("methylRawList", outList,treatment=methylObj@treatment)
                    myobj
                    }
)


#' Get methylated/unmethylated base counts for tilling windows 
#'
#' The function summarizes methylated/unmethylated base counts over tilling windows accross genome.  
#' This function can be used when differential methylated analysis is preferable to tilling windows instead of 
#' base pairs.
#'
#' @param methylObj \code{methylRaw} or \code{methylRawList} object containing base pair resolution methylation information
#' @param win.size an integer for the size of the tiling windows
#' @param step.size an integer for the step size of tiling windows
#' @param cov.bases minimum number of bases to be covered in a given window
#' @usage tileMethylCounts(methylObj,win.size=1000,step.size=1000,cov.bases=0)
#' @return \code{methylRaw} or \code{methylRawList} object
#' @export
#' @docType methods
#' @rdname tileMethylCounts-methods
setGeneric("tileMethylCounts", function(methylObj,win.size=1000,step.size=1000,cov.bases=0) standardGeneric("tileMethylCounts") )

#' @aliases tileMethylCounts,methylRaw-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(methylObj="methylRaw"),
                    function(methylObj,win.size,step.size,cov.bases){
                        
                      g.meth =as(methylObj,"GRanges")
                      chrs   =IRanges::levels(seqnames(g.meth))
                      widths =seqlengths(g.meth)
                      all.wins=GRanges()
                      for(i in 1:length(chrs))
                      {
                        max.length=max(IRanges::end(g.meth[seqnames(g.meth)==chrs[i],])) # get max length of feature covered chromosome
                        #get sliding windows with covered CpGs
                        numTiles=floor(  (max.length-(win.size-step.size) )/step.size )
                        temp.wins=GRanges(seqnames=rep(chrs[i],numTiles),ranges=IRanges(start=1+0:(numTiles-1)*step.size,width=rep(win.size,numTiles)) )
                        all.wins=c(all.wins,temp.wins)
                      }
                      regionCounts(methylObj,all.wins,cov.bases)

                      
})

#' @aliases tileMethylCounts,methylRawList-method
#' @rdname tileMethylCounts-methods
setMethod("tileMethylCounts", signature(methylObj="methylRawList"),
                    function(methylObj,win.size,step.size,cov.bases){

                    new.list=lapply(methylObj,tileMethylCounts,win.size,step.size,cov.bases) 
                    new("methylRawList", new.list,treatment=methylObj@treatment)
                    
})

# back up code:
# annotation of methylation states by genes
# including promoter, intron, exon, nearest cpgi and cpgi shore.


# valPerIndex=function(query, subject, col="cov"){
# cov.per=value.per.location(query,subject,col=col)
# all.ind=1:length(query)
# rm=data.frame(ind=all.ind[!all.ind %in% (cov.per[,1])],val=rep(0,length(all.ind[!all.ind %in% (cov.per[,1])]) ),numRows=rep(0,length(all.ind[!all.ind %in% (cov.per[,1])]) ) )
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


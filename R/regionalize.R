# functions to get regional counts

# get counts for a given regions


#
# S4 FUNCTIONS
#


#' GETs regional counts for given GRanges or GRangesList object
#'
#' Convert CpG dinucleotide methylation \code{methylRaw} or \code{methylRawList} object into regional counts for a given GRanges or GRangesList object.
#' @param methylObj a \code{methylRaw} or \code{methlRawList} object
#' @param regions a GRanges or GRangesList object.
#' 
#' @return RETURNS a new methylRaw or methylRawList object
#' 
#' @aliases regionCounts,-methods regionCounts,methylRawList-method
#' @aliases regionCounts,methylRaw,methyRawList,ANY-method
#' @export
#' @docType methods
#' @rdname regionCounts-methods
setGeneric("regionCounts", function(methylObj,regions) standardGeneric("regionCounts") )


# GETs regional counts for given GRanges object
# RETURNS a new methylRaw object
# @param methylObj a \code{methylRaw} object
# @param regions a GRanges object.
# @rdname regionCounts-methods
# @aliases regionCounts,methylRaw,ANY-method
setMethod("regionCounts", signature(methylObj="methylRaw",regions="GRanges"),
                    function(methylObj,regions){
                    
                      # overlap methylObj with regions
                      # convert methylObj to GRanges
                      mat=matchMatrix( findOverlaps(regions,as(methylObj,"GRanges")) )
                      
                      # create a temporary data.table row ids from regions and counts from methylObj
                      dt=data.table(id=mat[,1],methylObj[mat[,2],c(6,7,8)] )
                      
                      # use data.table to sum up counts per region
                      sum.dt=dt[,list(coverage=sum(coverage),
                               numCs   =sum(numCs),
                               numTs   =sum(numTs)),by=id] 
                      
                      temp.df=as.data.frame(regions) # get regions to a dataframe
                      
                      # look for values with "name" in it, eg. "tx_name" or "name"
                     # valuesList = names(values(regions))
                     # nameid = valuesList[grep (valuesList, pattern="name")]
                      
                      
                      #create a new methylRaw object to return
                      new.data=data.frame(id      =temp.df[sum.dt$id,"name"],
                                          chr     =temp.df[sum.dt$id,"seqnames"],
                                          start   =temp.df[sum.dt$id,"start"],
                                          end     =temp.df[sum.dt$id,"end"],
                                          strand  =temp.df[sum.dt$id,"strand"],
                                          coverage=sum.dt$coverage,
                                          numCs   =sum.dt$numCs,
                                          numTs   =sum.dt$numTs)
                    
                    new("methylRaw",new.data,sample.id=methylObj@sample.id,assembly=methylObj@assembly)
                      
                    }
)


# RETURNS a new methylRawList object
# gets regional counts for all elements in methylRawList for given regions
# MAKE SURE an element of the list, which will be a set of GRanges rows, are not on different chromsomes and strands
# Also, make sure id column of returned methylRaw object is unique
# you can add refseq id to the id column: chr.start.end.refseqid
# @param methylObj a \code{methylRaw} object
# @param regions a GRangesList object.
# @rdname regionCounts-methods
# @aliases regionCounts,methylRaw,ANY-method
# assume that each name of the element in the GRangesList is unique and 
setMethod("regionCounts", signature(methylObj="methylRaw",regions="GRangesList"),
                    function(methylObj,regions){
                                        
                      # overlap methylObj with regions
                      # convert methylObj to GRanges
                      mat=matchMatrix( findOverlaps(regions,as(methylObj,"GRanges")) )
                      
                      # create a temporary data.table row ids from regions and counts from methylObj
                      dt=data.table(id=mat[,1],methylObj[mat[,2],c(6,7,8)] )
                      
                      # use data.table to sum up counts per region
                      sum.dt=dt[,list(coverage=sum(coverage),
                               numCs   =sum(numCs),
                               numTs   =sum(numTs)),by=id] 
                      
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
                    
                    new("methylRaw",new.data,sample.id=methylObj@sample.id,assembly=methylObj@assembly)
  
                    })
# Note: some genes do not have intron, need to take care of it.


# GETs regional counts for given GRanges object
# RETURNS a new methylRawList object
# @param methylObj a \code{methylRawList} object
# @param regions a GRanges object.
# @rdname regionCounts-methods
# @aliases regionCounts,methylRawList,ANY-method
setMethod("regionCounts", signature(methylObj="methylRawList",regions="GRanges"),
                    function(methylObj,regions){
                     
                    outList=list()

                    for(i in 1:length(methylObj))
                    {
                        obj = regionCounts(methylObj = methylObj[[i]],regions=regions)
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
# @rdname regionCount-methods
# @aliases regionCounts,methylRawList,ANY-method
setMethod("regionCounts", signature(methylObj="methylRawList",regions="GRangesList"),
                    function(methylObj,regions){
                     
                    outList=list()

                    for(i in 1:length(methylObj))
                    {
                        obj = regionCounts(methylObj = methylObj[[i]],regions=regions)
                        outList[[i]] = obj
                    }
                    
                    myobj=new("methylRawList", outList,treatment=methylObj@treatment)
                    myobj
                    }
)


# get counts for genomic tiles
#genomeTileCounts()


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


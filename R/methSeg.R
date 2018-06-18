# segmentation functions 

#' Segment methylation or differential methylation profile
#' 
#' The function uses a segmentation algorithm (\code{\link[fastseg]{fastseg}}) 
#' to segment the methylation profiles. Following that, it uses
#' gaussian mixture modelling to cluster the segments into k components. This process
#' uses mean methylation value of each segment in the modeling phase. Each
#' component ideally indicates quantitative classification of segments, such
#' as high or low methylated regions.
#' 
#' @param obj \code{\link[GenomicRanges]{GRanges}}, \code{\link{methylDiff}}, 
#'        \code{\link{methylDiffDB}}, \code{\link{methylRaw}} or 
#'        \code{\link{methylRawDB}} . If the object is a 
#'        \code{\link[GenomicRanges]{GRanges}} it should have one meta column 
#'        with methylation scores and has to be sorted by position, 
#'        i.e. ignoring the strand information.
#' @param diagnostic.plot if TRUE a diagnostic plot is plotted. The plot shows
#'        methylation and length statistics per segment group. In addition, it 
#'        shows diagnostics from mixture modeling: the density function estimated 
#'        and BIC criterion used to decide the optimum number of components
#'        in mixture modeling.
#' @param join.neighbours if TRUE neighbouring segments that cluster to the same 
#'        seg.group are joined by extending the ranges, summing up num.marks and
#'        averaging over seg.means.
#' @param initialize.on.subset a numeric value indicating either percentage or
#'        absolute value of regions to subsample from segments before performing
#'        the mixture modeling. The value can be either between 0 and 1, 
#'        e.g. 0.1 means that 10% of data will be used for sampling, or be an 
#'        integer higher than 1 to define the number of regions to sample. 
#'        By default uses the whole dataset, which can be time consuming on 
#'        large datasets. (Default: 1)
#' @param ... arguments to \code{\link[fastseg]{fastseg}} function in fastseg 
#'        package, or to \code{\link[mclust]{densityMclust}}
#'        in Mclust package, could be used to fine tune the segmentation algorithm.
#'        E.g. Increasing "alpha" will give more segments. 
#'        Increasing "cyberWeight" will give also more segments."maxInt" controls
#'        the segment extension around a breakpoint. "minSeg" controls the minimum
#'        segment length. "G" argument
#'        denotes number of components used in BIC selection in mixture modeling.
#'        For more details see fastseg and Mclust documentation.    
#'        
#'               
#' @return A \code{\link[GenomicRanges]{GRanges}} object with segment 
#'         classification and information. 
#'        'seg.mean' column shows the mean methylation per segment.
#'        'seg.group' column shows the segment groups obtained by mixture modeling
#'               
#' @details      
#'        To be sure that the algorithm will work on your data, 
#'        the object should have at least 5000 records
#'        
#' @details After initial segmentation with fastseg(), segments are clustered 
#'        into self-similar groups based on their mean methylation profile 
#'        using mixture modeling. Mixture modeling estimates the initial 
#'        parameters of the distributions by using the whole dataset. 
#'        If "initialize.on.subset" argument set as described, the function 
#'        will use a subset of the data to initialize the parameters of the 
#'        mixture modeling prior to the Expectation-maximization (EM) algorithm 
#'        used by \code{\link{Mclust}} package.
#'  
#' @examples 
#' 
#' \donttest{
#'  download.file("https://dl.dropboxusercontent.com/u/1373164/H1.chr21.chr22.rds",
#'                destfile="H1.chr21.chr22.rds",method="curl")
#' 
#'  mbw=readRDS("H1.chr21.chr22.rds")
#' 
#'  # it finds the optimal number of componets as 6
#'  res=methSeg(mbw,diagnostic.plot=TRUE,maxInt=100,minSeg=10)
#' 
#'  # however the BIC stabilizes after 4, we can also try 4 componets
#'  res=methSeg(mbw,diagnostic.plot=TRUE,maxInt=100,minSeg=10,G=1:4)
#' 
#'  # get segments to BED file
#'  methSeg2bed(res,filename="H1.chr21.chr22.trial.seg.bed")
#' }
#' 
#' unlink(list.files(pattern="H1.chr21.chr22",full.names=TRUE))
#' 
#' @author Altuna Akalin, contributions by Arsene Wabo and Katarzyna
#' Wreczycka
#' 
#' @seealso \code{\link{methSeg2bed}}, \code{\link{joinSegmentNeighbours}} 
#' 
#' @export
#' @docType methods
#' @rdname methSeg       
methSeg<-function(obj, diagnostic.plot=TRUE, join.neighbours=FALSE,
                  initialize.on.subset=1, ...){
  
  dots <- list(...)  
  
  
  ##coerce object to granges
  if(class(obj)=="methylRaw" | class(obj)=="methylRawDB") {
    obj= as(obj,"GRanges")
    ## calculate methylation score 
    mcols(obj)$meth=100*obj$numCs/obj$coverage
    ## select only required mcol
    obj = obj[,"meth"]
  }else if (class(obj)=="methylDiff" | class(obj)=="methylDiffDB") {
    obj = as(obj,"GRanges")
    ## use methylation difference as score
    obj = obj[,"meth.diff"]
  }else if (class(obj) != "GRanges"){
    stop("only methylRaw or methylDiff objects ", 
         "or GRanges objects can be used in this function")
  }
  
  # destrand
  strand(obj) <- "*"
  
  ## check wether obj contains at least one metacol 
  if(ncol(elementMetadata(obj))<1)
    stop("GRanges does not have any meta column.")
  
  ## check wether obj contains is sorted by position
  if(is.unsorted(obj,ignore.strand=TRUE)) {
    obj <- sort(obj,ignore.strand=TRUE)
    message("Object not sorted by position, sorting now.")
  }
  
  ## check wether obj contains at least two ranges else stop
  if(length(obj)<=1)
    stop("segmentation requires at least two ranges.")
  
  # match argument names to fastseg arguments
  args.fastseg=dots[names(dots) %in% names(formals(fastseg)[-1] ) ]  
  
  # match argument names to Mclust
  args.Mclust=dots[names(dots) %in% names(formals(Mclust)[-1])  ]
  
  args.fastseg[["x"]]=obj
  
  # do the segmentation
  #seg.res=fastseg(obj)
  seg.res <- do.call("fastseg", args.fastseg)
  #seg.res <- do.call("fastseg2", args.fastseg)
  
  # stop if segmentation produced only one range
  if(length(seg.res)==1) {
    warning("segmentation produced only one range, no mixture modeling possible.")
    seg.res$seg.group <- "1"
    return(seg.res)
  }
  
  # if joining, do not show first clustering
  if(join.neighbours) {
    diagnostic.plot.old = diagnostic.plot
    diagnostic.plot = FALSE
  }
  
  if("initialization" %in% names(args.Mclust)){
    if("subset" %in% names(args.Mclust[["initialization"]])) {
      if(length(args.Mclust[["initialization"]][["subset"]]) < 9 ){
        stop("too few samples, increase the size of subset.") 
        }
      message(paste("initializing clustering with",
                    length(args.Mclust[["initialization"]][["subset"]]),
                    "segments."))
      initialize.on.subset = 1
    }
  }
  
  if(initialize.on.subset != 1 && initialize.on.subset > 0 ) {
    
    if( initialize.on.subset > 0 & initialize.on.subset < 1 )
      nbr.sample = floor(length(seg.res) * initialize.on.subset)
    
    if( initialize.on.subset > 1) 
      nbr.sample = initialize.on.subset
    
    if( nbr.sample < 9 ){stop("too few samples, increase the size of subset.") }
    
    message(paste("initializing clustering with",nbr.sample,"out of",length(seg.res),"total segments."))
    # estimate parameters for mclust
    sub <- sample(1:length(seg.res), nbr.sample,replace = FALSE)
    args.Mclust[["initialization"]]=list(subset = sub)
  }  
  
  
  # decide on number of components/groups
  args.Mclust[["score.gr"]]=seg.res
  args.Mclust[["diagnostic.plot"]]=diagnostic.plot
  dens=do.call("densityFind", args.Mclust  )
  
  
  # add components/group ids 
  mcols(seg.res)$seg.group=as.character(dens$classification)
  
  # if joining, show clustering after joining
  if(join.neighbours) {
    message("joining neighbouring segments and repeating clustering.")
    seg.res <- joinSegmentNeighbours(seg.res)
    diagnostic.plot <- diagnostic.plot.old
    
    # get the new density
    args.Mclust[["score.gr"]]=seg.res
    args.Mclust[["diagnostic.plot"]]=diagnostic.plot
    # skip second progress bar
    args.Mclust[["verbose"]]=FALSE
    dens=do.call("densityFind", args.Mclust  )
    
  }
  
  seg.res
}

# not needed
.methSeg<-function(score.gr,score.cols=NULL,...){
  #require(fastseg)
  
  
  if(!is.null(score.cols)){
    values(score.gr)=score.gr[,score.cols]
  }
  
  seg.res <- fastseg(score.gr,...)
  
}

# finds segment groups using mixture modeling
densityFind<-function(score.gr,diagnostic.plot=TRUE,...){
  dens = densityMclust(score.gr$seg.mean,... )
  
  if(diagnostic.plot){
    diagPlot(dens,score.gr)
  }
  dens
}

# Plot diagnostic graphs for segmentation results
diagPlot<-function(dens,score.gr){
  
  scores=score.gr$seg.mean
  par(mfrow=c(2,3))
  boxplot(
    lapply(1:dens$G,function(x) scores[dens$classification==x] ),
    horizontal=TRUE,main="methylation per group",xlab="methylation")
  
  boxplot(
    lapply(1:dens$G,function(x) log10(width(score.gr)[dens$classification==x]) ),
    horizontal=TRUE,main="segment length per group",
    xlab="log10(length) in bp ",outline=FALSE)
  
  
  #lapply(1:dens$G,function(x) mean(width(score.gr)[dens$classification==x] ))
  
  barplot(table(dens$classification),xlab="segment groups",
          ylab="number of segments")
  plot(dens,what="density")  
  plot(dens,what="BIC",legendArgs=list(cex=0.75))  
  
  par(mfrow=c(1,1))
  
}

#' Export segments to BED files
#' 
#' The segments are color coded based on their score (methylation or differential
#' methylation value). They are named by segment group (components in mixture modeling)
#' and the score in the BED file is obtained from 'seg.mean' column of segments
#' object.
#' 
#' @param segments \code{\link[GenomicRanges]{GRanges}} object with segment
#'        classification and information. This should be the result of 
#' \code{\link{methSeg}} function
#' @param filename name of the output data
#' @param trackLine UCSC browser trackline
#' @param colramp color scale to be used in the BED display
#'        defaults to gray,green, darkgreen scale.
#' 
#' @return A BED files with the segmented data
#' which can be visualized in the UCSC browser 
#' 
#' @seealso \code{\link{methSeg}}
#' 
#' @export
#' @docType methods
#' @rdname methSeg2bed
methSeg2bed<-function(segments,filename,
                      trackLine="track name='meth segments' description='meth segments' itemRgb=On",
                      colramp=colorRamp(c("gray","green", "darkgreen"))
){
  if(class(segments)!="GRanges"){
    stop("segments object has to be of class GRanges")
  }
  
  ## case if only one line is exported
  if(is.null(colramp) | length(segments)==1){
    trackLine <- gsub(pattern = "itemRgb=On",replacement = "",x = trackLine)
  } else {
    #require(rtracklayer)
    ramp <- colramp
    mcols(segments)$name=as.character(segments$seg.group)
    
    #scores=(segments$seg.mean-min(segments$seg.mean))/(max(segments$seg.mean))-(min(segments$seg.mean))
    scores=(segments$seg.mean-min(segments$seg.mean))/(max(segments$seg.mean)-
                                                         min(segments$seg.mean))
    
    mcols(segments)$itemRgb= rgb(ramp(scores), maxColorValue = 255)     
  }
  
  
  #strand(segments)="."
  score(segments)=segments$seg.mean
  
  if(is.null(trackLine)){
    
    export.bed(segments,filename)
  }else{
    export.bed(segments,filename,
               trackLine=as(trackLine, "BasicTrackLine"))
  }
}

# this could be used to avoid total dependece on GenomicRanges
# currently it has to be listed on depends section of DESCRIPTION
# but doing the depend thing is cleaner
#fastseg2<-function(x, type = 1, alpha = 0.05, segMedianT, minSeg = 4, 
#eps = 0, delta = 5, maxInt = 40, squashing = 0, cyberWeight = 10){
#  
#  y <- split(x, as.character(seqnames(x)))
#  nbrOfSeq <- length(y)
#  res02 <- list()
#  for (seq in seq_len(nbrOfSeq)) {
#    x <- y[[seq]]
#    res <- list()
#    for (sampleIdx in seq_len(ncol(elementMetadata(x)))) {
#      z01 <- elementMetadata(x)[[sampleIdx]]
#      sample <- names(elementMetadata(x))[sampleIdx]
#      resTmp <- fastseg:::segmentGeneral(z01, type, alpha, segMedianT, 
#                               minSeg, eps, delta, maxInt, squashing, 
#                               cyberWeight)$finalSegments
#      resTmp$sample <- sample
#      res[[sampleIdx]] <- resTmp
#    }
#    res <- do.call("rbind", res)
#    res$num.mark <- res$end - res$start
#    chrom <- as.character(seqnames(x)[1])
#    start <- start(x)[res$start]
#    end <- start(x)[res$end] + width(x)[res$end] - 1
#    resX <- data.frame(ID = res$sample, chrom = chrom, 
#                   loc.start = start, loc.end = end, num.mark = res$num.mark, 
#                   seg.mean = res$mean, startRow = res$start, endRow = res$end, 
#                   stringsAsFactors = FALSE)
#    resX <- resX[order(resX$chrom, resX$loc.start), ]
#    res02[[seq]] <- resX
#  }
#  res03 <- do.call("rbind", res02)
#  finalRes <- GRanges(seqnames = Rle(res03$chrom),
#                      ranges = IRanges(start = res03$loc.start, 
#                                      end = res03$loc.end), 
#                      ID = res03$ID, num.mark = res03$num.mark, 
#                      seg.mean = res03$seg.mean, startRow = res03$startRow, 
#                      endRow = res03$endRow)
#  return(finalRes)
#}


#' Join directly neighbouring segments produced by methSeg
#' 
#'
#' Segmentation and clustering are two distinct steps in methSeg(), 
#' leading to adjacent segments of the same class. 
#' This leads to a bias segment length distributions, 
#' which is removed by joining those neighbours.
#'
#' @param res A \code{\link[GenomicRanges]{GRanges}} object with segment 
#'         classification and information prudoced by \code{\link{methSeg}} 
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object with segment 
#'         classification and information. 
#' 
#' @seealso \code{\link{methSeg}}
#' 
#' @export
#' @importFrom data.table copy ":="
# @noRd
# @examples
joinSegmentNeighbours <- function(res) {
  
  # require(data.table)
  
  if (length(unique(seqnames(res))) > 1 ) {
    ## call recursively for multiple chromosomes
    gr <- lapply(split(res,seqnames(res)),joinSegmentNeighbours)
    gr <- do.call(c, unlist(gr,use.names = FALSE) )
    return( gr )
  } 
  else if (length(res)<=1) {
    ## return object if it has only one range
    return(res)
  }
  else{
    ## compute run-length-enconding of groups, 
    ## which is higher than 1 if neighbours are in same group
    group.neighbours <- rle(res$seg.group)
    N = length(group.neighbours$lengths)
    # res_joined <- vector(mode="list", length=N)
    
    ## k[i] is supposed to store the orginal index of the range i in res  
    k <- numeric(N)
    k[1] <- 1
    
    ## l[i] is either 0 if adjacent groups are distinct and 
    ## higher or equal to 1 if groups are the same 
    l <- group.neighbours$lengths - 1
    
    ## for some positions k[j] we have to skip l[j] index positions 
    ## since neighbouring ranges can be merged
    for ( i in 2:N ) { k[i] = k[i-1] + l[ i -1] + 1 }
    #toJoin <- l==0
    
    # print(paste("k=",k,"l=",l))
    
    res_dt <- copy(as.data.table(res))
    
    ## initialize "global variables" to silence R CMD check
    ID=endRow=num.mark=seg.group=seg.mean=startRow=NULL
    
    ## now we just merge ranges, that had a non-zero value in l
    for (i in which(l!=0)) {
      res_dt[k[i]:(k[i]+l[i]),":="(seqnames=unique(seqnames),
                                   start=min(start),
                                   end=max(end),
                                   strand = "*",
                                   width = sum(as.numeric(width)),
                                   ID = unique(ID),
                                   num.mark = sum(as.numeric(num.mark)),
                                   seg.mean = mean(seg.mean),
                                   startRow = min(startRow),
                                   endRow = max(endRow),
                                   seg.group = unique(seg.group))]
    }
    
    res_dt <- unique(res_dt)
    
    return(makeGRangesFromDataFrame(res_dt,keep.extra.columns = TRUE))
  }
}



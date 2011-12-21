# PART THAT DEALS with annotation of CpGs and differential methylation events

##############################################################################
# SECTION 1:
# reading annotation to GRanges
# makes GRanges object from a given bed6 or bed12 file to granges object
##############################################################################

#######################################
# SECTION 1: S3 functions
#######################################

# extracts exons from a bed12 file and puts them into GRanges object
# done in pure R
bed12.to.exons<-function(ref)
{
  #colClasses=c("factor","integer","integer","factor","integer","factor","integer","integer","integer","integer","factor","factor") 
  #ref         =read.table(bed.file,header=F,skip=skip,colClasses = colClasses,stringsAsFactors=F) # give colClasses for fast read-in
  ref=unique(ref)
  # apply strsplit on columns 11 and 12 to extract exon start positions and exon sizes
  b.start.size=cbind(as.integer(unlist(strsplit(as.character(ref$V12),"," ))),as.integer(unlist(strsplit(as.character(ref$V11),"," ))))
  rep.ref     =ref[rep(1:nrow(ref),ref[,10]),] # replicate rows occurs as many as its exon number
  #exon.id     =unlist(sapply(ref[,10],function(x) 1:x));rep.ref$V5=exon.id
  exon.id     =unlist( mapply( function(x,y) ifelse(x=="+",return(1:y),return(y:1) ),ref[,6],ref[,10]  ) );rep.ref$V5=exon.id

  rep.ref$V3  =rep.ref$V2+b.start.size[,1]+b.start.size[,2] # calculate exon start and ends
  rep.ref$V2  =rep.ref$V2+b.start.size[,1]
  
  #return(rep.ref[,1:6]) if you want to return as a data.frame

#  if(is.null(chrom.file))
#  {
    strands=as.character(rep.ref$V6)
    strands[strands=="."]="*"
    grange.exons=GenomicRanges::GRanges(seqnames=as.character(rep.ref$V1),
          ranges=IRanges(start=rep.ref$V2+1, end=rep.ref$V3),
          strand=strands, score=rep.ref$V5,name=rep.ref$V4)
    return(grange.exons)
#  }
  
#  #read table for chromosome info
#  chrom=read.table(chrom.file,header=F)
#  chrom=chrom[chrom[,1] %in% rep.ref$V1,1:2]
#  seqlengths=chrom[,2];names(seqlengths)=as.character(chrom[,1])
  # do adjustments so that trans and exons have the same numerical id
#  strands=as.character(rep.ref$V6)
#  strands[strands=="."]="*"
#  grange.exons=GRanges(seqnames=as.character(rep.ref$V1),
#          ranges=IRanges(start=rep.ref$V2+1, end=rep.ref$V3),
#          strand=strands, score=rep.ref$V5,name=rep.ref$V4,seqlengths=seqlengths)
#  return(grange.exons)

} 

# extracts exons from a bed12 file and puts them into GRanges object
# done in pure R
bed12.to.introns<-function(ref)
{
  #colClasses=c("factor","integer","integer","factor","integer","factor","integer","integer","integer","integer","factor","factor") 
  #ref         =read.table(bed.file,header=F,skip=skip,colClasses = colClasses,stringsAsFactors=F) # give colClasses for fast read-in

  #remove the genes with one exon only (they won't have any introns)
  ref=ref[ref[,10]>1,]
  ids=paste(ref[,1],ref[,2],ref[,3],ref[,4],sep="")
  ref=cbind(ref,id=ids)
  ref=unique(ref)
  # apply strsplit on columns 11 and 12 to extract exon start positions and exon sizes
  b.start.size=cbind(as.integer(unlist(strsplit(as.character(ref$V12),"," ))),as.integer(unlist(strsplit(as.character(ref$V11),"," ))))
  rep.ref     =ref[rep(1:nrow(ref),ref[,10]),] # replicate rows occurs as many as its exon number
  #exon.id     =unlist(sapply(ref[,10],function(x) 1:x));rep.ref$V5=exon.id
  exon.id     =unlist( mapply( function(x,y) ifelse(x=="+",return(1:y),return(y:1) ),ref[,6],ref[,10]  ) );rep.ref$V5=exon.id
  rep.ref$V3  =rep.ref$V2+b.start.size[,1]+b.start.size[,2] # calculate exon start and ends
  rep.ref$V2  =rep.ref$V2+b.start.size[,1]
  rep.ref     =rep.ref[,c(1:6,13)]
  
  # now try to get the exons by cbinding consecutive exons
  temp.ref    =cbind(rep.ref[1:(nrow(rep.ref)-1),],rep.ref[2:nrow(rep.ref),])
  temp.ref    =temp.ref[temp.ref[,7]==temp.ref[,14],]
  temp.ref[,2]=temp.ref[,3]
  temp.ref[,3]=temp.ref[,9]
  rep.ref     =temp.ref[,1:6]
  rep.ref[rep.ref[,6]=="-",5]=rep.ref[rep.ref[,6]=="-",5]-1 # subtract 1 from - strand exon numbers
    
#  if(is.null(chrom.file))
#  {
    strands=as.character(rep.ref$V6)
    strands[strands=="."]="*"
    grange.exons=GenomicRanges::GRanges(seqnames=as.character(rep.ref$V1),
          ranges=IRanges(start=rep.ref$V2+1, end=rep.ref$V3),
          strand=strands, score=rep.ref$V5,name=rep.ref$V4)
    return(grange.exons)
#  }
  
#  #read table for chromosome info
#  chrom=read.table(chrom.file,header=F)
#  chrom=chrom[chrom[,1] %in% rep.ref$V1,1:2]
#  seqlengths=chrom[,2];names(seqlengths)=as.character(chrom[,1])
#  # do adjustments so that trans and exons have the same numerical id
#  strands=as.character(rep.ref$V6)
#  strands[strands=="."]="*"
#  grange.exons=GRanges(seqnames=as.character(rep.ref$V1),
#          ranges=IRanges(start=rep.ref$V2+1, end=rep.ref$V3),
#          strand=strands, score=rep.ref$V5,name=rep.ref$V4,seqlengths=seqlengths)
#  return(grange.exons)

} 

# checks the validity of the bed data.frame if it is a legitimate bed columns
check.bed.validity<-function(bed.df,type="none")
{
  # does it have at least 3 columns
  num.col=(ncol(bed.df)>=3)
  
  # does it have 1st and 2nd column as numbers
  col1.2= (is.integer(bed.df[,2]) & is.integer(bed.df[,3]) )
  
  # does it have chr string in 1st column
  chr= sum(grepl("chr",bed.df[,1]) )

  if(type=="exons" | type=="introns")
  {
    #does it have 12>= columns
    ex=(ncol(bed.df)>=12) 
    
    return(num.col & col1.2 & chr & ex)
  }
  return(num.col & col1.2 & chr )

}

#######################################
# SECTION 1: S4 functions
#######################################

#' convert a data frame read-in from a bed file to a GRanges object
#'  
#' @param bed  a data.frame where column order and content resembles a bed file with 12 columns
#' @usage convert.bed.df(bed)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname convert.bed.df-methods
setGeneric("convert.bed.df",function(bed) standardGeneric("convert.bed.df"))

#' @aliases convert.bed.df,data.frame-method
#' @rdname convert.bed.df-methods
setMethod("convert.bed.df" ,signature(bed = "data.frame" ),
                            function(bed){
                              
                            if(! check.bed.validity(bed)){stop("this is not a valid bed file")}
                            if(ncol(bed)>=6){
                              # check if 6th column is really strand
                              if( sum( unique(bed[,6]) %in% c("+","-",".") ) != length(unique(bed[,6])) ){stop("Strand column of the bed file or data.frame is wrong")}
                              #convert to granges
                              strands=as.character(bed$V6)
                              strands[strands=="."]="*"
                              grange=GRanges(seqnames=as.character(bed$V1),
                                            ranges=IRanges(start=bed$V2+1, end=bed$V3),
                                            strand=strands, score=bed$V5,name=bed$V4)
                            }
                            if(ncol(bed)==5){
                              grange=GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2+1, end=bed$V3),strand=rep("*",nrow(bed)), score=bed$V5,name=bed$V4 )
                            }
                            if(ncol(bed)==4){
                              grange=GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2+1, end=bed$V3),strand=rep("*",nrow(bed)),name=bed$V4)
                            }    
                            if(ncol(bed)==3){
                              grange=GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2+1, end=bed$V3),strand=rep("*",nrow(bed)))
                            }                           
                            return(grange)
})

#' convert a data frame read-in from a bed file to a GRanges object for exons
#'  
#' @param bed.df  a data.frame where column order and content resembles a bed file with 12 columns
#' @usage convert.bed2exons(bed.df)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname convert.bed2exons-methods
setGeneric("convert.bed2exons",function(bed.df) standardGeneric("convert.bed2exons"))

#' @aliases convert.bed2exons,data.frame-method
#' @rdname convert.bed2exons-methods
setMethod("convert.bed2exons" ,signature(bed.df = "data.frame" ),
                            function(bed.df){
                              
                            if(! check.bed.validity(bed.df,"exon")){stop("this is not a valid bed file")}
                            bed12.to.exons(bed.df)
})

#' convert a data frame read-in from a bed file to a GRanges object for introns
#'  
#' @param bed.df  a data.frame where column order and content resembles a bed file with 12 columns
#' @usage convert.bed2introns(bed.df)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname convert.bed2introns-methods
setGeneric("convert.bed2introns",function(bed.df) standardGeneric("convert.bed2introns"))

#' @aliases convert.bed2introns,data.frame-method
#' @rdname convert.bed2introns-methods
setMethod("convert.bed2introns" ,signature(bed.df = "data.frame" ),
                            function(bed.df){
                              
                            if(! check.bed.validity(bed.df,"exon")){stop("this is not a valid bed file")}                      
                            bed12.to.introns(bed.df)
})


#' read a bed file and convert it to GRanges
#'  
#' @param location  location of the file, a character string such as: "/home/user/my.bed"
#' @param remove.unsual if TRUE(default) remove the chromomesomes with unsual names, mainly random chromsomes etc
#'
#' @usage read.bed(location,remove.unsual=T)
#' @return \code{\link{GRanges}} object
#'
#' @note one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname read.bed-methods
setGeneric("read.bed", function(location,remove.unsual=T) standardGeneric("read.bed"))

#' @aliases read.bed,character-method
#' @rdname read.bed-methods
setMethod("read.bed", signature(location = "character"),#remove.unsual="logical" ),
                    function(location,remove.unsual){
                    
                    # find out if there is a header, skip 1st line if there is a header
                    f.line=readLines(con = location, n = 1)
                    skip=0
                    if(grepl("^track",f.line))(skip=1)
                    
                    # read bed6
                    bed=.readTableFast(location,header=F,skip=skip)                    
                    if(remove.unsual){ bed=bed[grep("_", as.character(bed[,1]),invert=T),] }
                    convert.bed.df(bed)
                    
})

#' function reading exon intron, promoter structure from a given bed file
#'
#' @param location location of the bed file with 12 or more columns
#' @param remove.unsual remove the chromomesomes with unsual names, mainly random chromsomes etc
#' @param up.flank  up-stream from TSS to detect promoter boundaries
#' @param down.flank down-stream from TSS to detect promoter boundaries
#' @param unique.prom     get only the unique promoters, promoter boundaries will not have a gene name if you set this option to be TRUE
#' @usage read.transcript.features(location,remove.unsual=TRUE,up.flank=1000,down.flank=1000,unique.prom=TRUE)
#' @return a \code{\link{GRangesList}} containing locations of exon/intron/promoter/TSS
#' @note  one bed track per file is only accepted, the bed files with multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname read.transcript.features-methods
setGeneric("read.transcript.features", function(location,remove.unsual=TRUE,up.flank=1000,down.flank=1000,unique.prom=TRUE) standardGeneric("read.transcript.features"))

#' @aliases read.transcript.features,character-method
#' @rdname read.transcript.features-methods
setMethod("read.transcript.features", signature(location = "character"),#,remove.unsual="logical",up.flank="numeric",down.flank="numeric",unique.prom="logical" ),
                    function(location,remove.unsual,up.flank ,down.flank ,unique.prom){
                      
                    # find out if there is a header, skip 1st line if there is a header
                    f.line=readLines(con = location, n = 1)
                    skip=0
                    if(grepl("^track",f.line))(skip=1)
                    
                    # read bed6
                    bed=.readTableFast(location,header=F,skip=skip)                    
                    if(remove.unsual){ bed=bed[grep("_", as.character(bed[,1]),invert=T),] }
                    
                    introns=convert.bed2introns(bed)
                    exons  =convert.bed2exons(bed)
                    
                    # get the locations of TSSes
                    tss=bed
                    #  + strand
                    tss[tss$V6=="+",3]=tss[tss$V6=="+",2]
                  
                    #  - strand
                    tss[tss$V6=="-",2]=tss[tss$V6=="-",3]
                    tssg=GRanges(seqnames=as.character(tss$V1),
                            ranges=IRanges(start=tss$V2, end=tss$V3),
                            strand=as.character(tss$V6),score=rep(0,nrow(tss)),name=tss$V4)
                    
                    # get the locations of promoters
                    #  + strand
                    bed[bed$V6=="+",3]=bed[bed$V6=="+",2]+down.flank
                    bed[bed$V6=="+",2]=bed[bed$V6=="+",2]-up.flank
                  
                    #  - strand
                    bed[bed$V6=="-",2]=bed[bed$V6=="-",3]-down.flank
                    bed[bed$V6=="-",3]=bed[bed$V6=="-",3]+up.flank
                    
                    if(! unique.prom)
                    {
                      prom.df=(bed[,c(1,2,3,4,6)])
                      prom=GRanges(seqnames=as.character(prom.df$V1),
                            ranges=IRanges(start=prom.df$V2, end=prom.df$V3),
                            strand=as.character(prom.df$V6),score=rep(0,nrow(prom.df)),name=prom.df$V4)
                    }else{
                      prom.df=unique(bed[,c(1,2,3,6)])
                      prom=GRanges(seqnames=as.character(prom.df$V1),
                            ranges=IRanges(start=prom.df$V2, end=prom.df$V3),
                            strand=as.character(prom.df$V6),score=rep(0,nrow(prom.df)),name=rep(".",nrow(prom.df)) )
                      
                    }
                  
                    GRangesList(exons=exons,introns=introns,promoters=prom,TSSes=tssg)
})

#' a function to get upstream and downstream adjecent regions to a genomic feature such as CpG islands
#' 
#' @param grange GRanges object for the feature
#' @param flank  number of basepairs for the flanking regions
#' @param clean  If set to TRUE, flanks overlapping with other main features will be trimmed, and overlapping flanks will be removed
#'        this will remove multiple counts when other features overlap with flanks
#'
#' @usage getFlanks(grange,flank=2000,clean=T)
#' @return GRanges object for flanking regions
#' @export
#' @docType methods
#' @rdname getFlanks-methods
setGeneric("getFlanks", function(grange,flank=2000,clean=T) standardGeneric("getFlanks"))

#' @aliases getFlanks,GRanges-method
#' @rdname getFlanks-methods
setMethod("getFlanks", signature(grange= "GRanges"),
                    function(grange,flank=2000,clean=T){
          
                    shores=c( IRanges::flank(grange,flank),IRanges::flank(grange,flank,FALSE) )
                    if(clean){
                      shores=IRanges::reduce(IRanges::setdiff(shores, grange)) # ,erge overlapping shores remove CpG coordinates from all shores, das ist so cool!!
                    }
                    shores
})

#' a function to read-in genomic features and their upstream and downstream adjecent regions such as CpG islands and their shores
#'
#' @param location for the bed file of the feature 
#' @param flank    number of basepairs for the flanking regions
#' @param clean    If set to TRUE, flanks overlapping with other main features will be trimmed
#' @param remove.unsual  remove chromsomes with unsual names random, Un and antyhing with "_" character
#' @param feature.flank.name the names for feature and flank ranges, it should be a character vector of length 2. example: c("CpGi","shores")
#' @usage  read.feature.flank(location,remove.unsual=T,flank=2000,clean=T,feature.flank.name=NULL)
#' @return a GenomicRangesList contatining one GRanges object for flanks and one for GRanges object for the main feature.
#'   NOTE:This can not return a GRangesList at the moment because flanking regions do not have to have the same column name as the feature.
#'   GRangesList elements should resemble eachother in the column content. We can not satisfy that criteria for the flanks
#'
#' @export
#' @docType methods
#' @rdname read.feature.flank-methods
setGeneric("read.feature.flank", function(location,remove.unsual=T,flank=2000,clean=T,feature.flank.name=NULL) standardGeneric("read.feature.flank") )

#' @aliases read.feature.flank,character-method
#' @rdname read.feature.flank-methods
setMethod("read.feature.flank", signature(location = "character"),
                    function(location,remove.unsual,flank ,clean,feature.flank.name){
                    feat=read.bed(location,remove.unsual)
                    flanks=getFlanks(feat,flank=flank,clean=clean)
                    x=GenomicRangesList(features=feat,flanks=flanks)
                    if(!is.null(feature.flank.name) & length(feature.flank.name)==2)
                    {
                      names(x)=feature.flank.name
                    }
                    x
})

##############################################################################
# SECTION 2:
# annotate granges objects with annotations that read-in and converted to GRanges objects
##############################################################################

#######################################
# SECTION 2: Define new classes
#######################################

# A set of objects that will hold statistics about feature and annotation overlap

#' An S4 class that information on overlap of target features with annotation features  
#'
#' This object is desgined to hold statistics and information about genomic feature overlaps
#'          
#' @section Slots:\describe{
#'                  \item{\code{members}}{a matrix showing overlap of target features with annotation genomic features}
#'
#'                  \item{\code{annotation}}{a named vector of percentages}
#'
#'                  \item{\code{precedence}}{a named vector of percentages}
#'
#'                  \item{\code{num.hierarchica}}{vector}
#'
#'                  \item{\code{no.of.OlapFeat}}{vector}
#'
#'                  \item{\code{perc.of.OlapFeat}}{vector}
#' }
#' @name annotationByFeature-class
#' @rdname annotationByFeature-class
#' @export
setClass("annotationByFeature", representation(members         ="matrix",
                                      annotation      ="numeric",
                                      precedence    ="numeric",
                                      num.annotation  ="numeric",
                                      num.precedence="numeric",
                                      no.of.OlapFeat  ="numeric",
                                      perc.of.OlapFeat="numeric"))

#' An S4 class that information on overlap of target features with annotation features  
#'
#' This object is desgined to hold statistics and information about genomic feature overlaps
#'          
#' @section Slots:\describe{
#'                  \item{\code{members}}{a matrix showing overlap of target features with annotation genomic features}
#'
#'                  \item{\code{annotation}}{a named vector of percentages}
#'
#'                  \item{\code{precedence}}{a named vector of percentages}
#'
#'                  \item{\code{num.hierarchica}}{vector}
#'
#'                  \item{\code{no.of.OlapFeat}}{vector}
#'
#'                  \item{\code{perc.of.OlapFeat}}{vector}
#'
#'                  \item{dist.to.TSS}{a data frame showing distances to TSS and gene/TSS names and strand}
#' }
#' @name annotationByGenicParts-class
#' @rdname annotationByGenicParts-class
#' @export
setClass("annotationByGenicParts", representation(dist.to.TSS   ="data.frame"),contains="annotationByFeature")


#new.obj=new("annotationByGenicParts",
#            members=matrix(c(1,2,3,4)),annotation=c(1,2,0,3,4),
#            precedence=c(a=1,b=2,c=0,d=3,e=4),
#            num.annotation  =c(1,2,0,3,4),
#            num.precedence=c(1,2,0,3,4),
#            no.of.OlapFeat  =c(1,2,0,3,4),
#            perc.of.OlapFeat=c(1,2,0,3,4),
#            dist.to.TSS     =c(1,2,0,3,4) )
                                    
#' @rdname show-methods
#' @aliases show,annotationByGenicParts-method
setMethod("show", "annotationByGenicParts", function(object) {
  
  cat("summary of target set annotation with genic parts\n");cat(nrow(object@members));cat(" rows in target set\n--------------\n")
  cat("--------------\n")
  cat("percentage of target features overlapping with annotation :\n");print(object@annotation);cat("\n\n")
  cat("percentage of target features overlapping with annotation (with promoter>exon>intron precedence) :\n"); print(object@precedence) ;cat("\n\n")
  cat("percentage of annotation boundaries with feature overlap :\n");print(object@perc.of.OlapFeat);cat("\n\n")  
  cat("summary of distances to the nearest TSS :\n")
  print(summary(abs(object@dist.to.TSS[,2])))
})


#' @rdname show-methods
#' @aliases show,annotationByFeature-method
setMethod("show", "annotationByFeature", function(object) {
  
  cat("summary of target set annotation with feature annotation\n");cat(nrow(object@members));cat(" rows in target set\n--------------\n")
  cat("--------------\n")
  cat("percentage of target features overlapping with annotation :\n");print(object@annotation);cat("\n\n")
  cat("percentage of annotation boundaries with feature overlap :\n");print(object@perc.of.OlapFeat);cat("\n\n")  
})

#######################################
# SECTION 2: S3 FUNCTIONS 
# these shouldn't be exported
#######################################

annotate.gr.WithGenicParts<-function(gr,prom,exon,intron,strand=F)
{
  #require(GenomicRanges)
  if( ! strand){strand(gr)="*"}
  memb=data.frame(matrix(rep(0,length(gr)*3),ncol=3) )  ;colnames(memb)=c("prom","exon","intron")
  memb[countOverlaps(gr,prom)>0,1]=1
  memb[countOverlaps(gr,exon)>0,2]=1
  memb[countOverlaps(gr,intron)>0,3]=1

  annotation=c(promoter  =100*sum(memb$prom>0)/nrow(memb) ,
               exon      =100*sum(memb$exon>0)/nrow(memb) ,
               intron    =100*sum(memb$intron>0)/nrow(memb) ,
               intergenic=100*sum(rowSums(memb)==0)/nrow(memb) )
  
  num.annotation=c( promoter  =sum(memb$prom>0) ,
                    exon      =sum(memb$exon>0) ,
                    intron    =sum(memb$intron>0) ,
                    intergenic=sum(rowSums(memb)==0) )

  precedence=c(promoter  =100*sum(memb$prom>0)/nrow(memb) ,
               exon      =100*sum(memb$exon>0 & memb$prom==0)/nrow(memb) ,
               intron    =100*sum(memb$intron>0 & memb$exon==0 & memb$prom==0)/nrow(memb) ,
               intergenic=100*sum(rowSums(memb)==0)/nrow(memb) )

  num.precedence=c( promoter  =sum(memb$prom>0) ,
                      exon      =sum(memb$exon>0 & memb$prom==0),
                      intron    =sum(memb$intron>0 & memb$exon==0 & memb$prom==0) ,
                      intergenic=sum(rowSums(memb)==0) )  
  
  numberOfOlapFeat=c(promoter=sum(countOverlaps(prom,gr)>0),
    exon=sum(countOverlaps(exon,gr)>0),
    intron=sum(countOverlaps(intron,gr)>0) )
  percOfOlapFeat =100*numberOfOlapFeat/c(length(prom),length(exon),length(intron))
  return(list(members=memb,annotation=annotation,precedence=precedence,num.annotation=num.annotation,
              num.precedence=num.precedence,
              numberOfOlapFeat=numberOfOlapFeat,percOfOlapFeat=percOfOlapFeat) )
  
}

distance2nearestFeature<-function(g.idh,tss)
{
  #require(GenomicRanges)
  elementMetadata(g.idh)=DataFrame(elementMetadata(g.idh),orig.row=1:length(g.idh))
  id.col    =ncol( elementMetadata(g.idh))+5 # get the row number column

  tss.id.col=ncol( elementMetadata(g.idh))+5+7 # get the id column for tss in the merged data.frame below
  strnd.col=ncol( elementMetadata(g.idh))+5+7-2 # get the id column for tss
 
  met.tss   = .nearest.2bed(g.idh, tss) # a merged data.frame is returned by this function

  dist.col=ncol(met.tss)
  met.tss[met.tss$end<met.tss$start.y & met.tss$strand.y=="+",dist.col] = -1* met.tss[met.tss$end<met.tss$start.y & met.tss$strand.y=="+",dist.col]

  met.tss[met.tss$end>met.tss$start.y & met.tss$strand.y=="-",dist.col] = -1* met.tss[met.tss$end>met.tss$start.y & met.tss$strand.y=="-",dist.col]

  res=met.tss[order(met.tss[,id.col]),c(id.col,dist.col,tss.id.col,strnd.col)]
  names(res)=c("target.row","dist.to.feature"   ,    "feature.name"  ,   "feature.strand")

  return(res)
}


# get the nearest features between a subject and query bed file
# get grange objects in bed zormat
# g.bed is GenomicRanges object, biatch!!!
# subject is also GenomicRanges object
.nearest.2bed<-function(g.bed,subject)
{
  chrs1=IRanges::levels(seqnames(g.bed))
  chrs2=IRanges::levels(seqnames(subject))
  chrs=chrs1[chrs1 %in% chrs2]
  res.df=NULL
  for(i in 1:length(chrs))
  {
    # find the nearest for this range
    ind  =nearest(ranges(g.bed[seqnames(g.bed)==chrs[i],]),ranges(subject[seqnames(subject)==chrs[i],]))
    sbdf1 =IRanges::as.data.frame(g.bed[seqnames(g.bed)==chrs[i],] )
    sbdf2 =IRanges::as.data.frame(subject[seqnames(subject)==chrs[i],] )
    sbdf2 =sbdf2[ind,]
    names(sbdf2)=paste(names(sbdf2),".y",sep="")
    temp  =cbind(sbdf1,sbdf2)
    res.df=rbind(res.df,temp)
  }

  res.dist=rep(0,nrow(res.df))
  dist1=abs(res.df$start-res.df$end.y)+1
  dist2=abs(res.df$start.y-res.df$end)+1
  totti=res.df$width + res.df$width.y #total length
  res.dist[pmax(dist1,dist2)>totti]=pmin(dist1,dist2)[pmax(dist1, dist2) > totti] 
  res.df=cbind(res.df,dist=res.dist)
  return(res.df)
}





#######################################
# SECTION 2: S4 FUNCTIONS 
#######################################


#' function to annotate given GRanges object with promoter,exon,intron & intergenic values
#'
#' @param target      : a granges object storing chromosome locations to be annotated
#' @param GRangesList.obj  : A GRangesList object containing GRanges object for promoter,exons,introns and TSSes, or simply output of read.transcript.features function
#' @param strand           : If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
#' @usage annotate.WithGenicParts(target,GRangesList.obj,strand=F)
#' @return \code{annotationByGenicParts} object
#' 
#' @export
#' @docType methods
#' @rdname annotate.WithGenicParts-methods
setGeneric("annotate.WithGenicParts", function(target,GRangesList.obj,strand=F) standardGeneric("annotate.WithGenicParts") )

#' @aliases annotate.WithGenicParts,GRanges,GRangesList-method
#' @rdname annotate.WithGenicParts-methods
setMethod("annotate.WithGenicParts", signature(target= "GRanges",GRangesList.obj="GRangesList"),
                    function(target,GRangesList.obj,strand){
                      
                      a.list    =annotate.gr.WithGenicParts(target,GRangesList.obj$promoters,GRangesList.obj$exons,GRangesList.obj$introns,strand=strand)
                      dist2TSS  =distance2nearestFeature(target,GRangesList.obj$TSSes)

                      new("annotationByGenicParts",
                                  members         =as.matrix(a.list$members),
                                  annotation      =a.list$annotation,
                                  precedence    =a.list$precedence,
                                  num.annotation  =a.list$num.annotation,
                                  num.precedence=a.list$num.precedence,
                                  no.of.OlapFeat  =a.list$numberOfOlapFeat,
                                  perc.of.OlapFeat=a.list$percOfOlapFeat,
                                  dist.to.TSS     = dist2TSS )
})

#' @aliases annotate.WithGenicParts,methylDiff,GRangesList-method
#' @rdname annotate.WithGenicParts-methods
setMethod("annotate.WithGenicParts", signature(target = "methylDiff",GRangesList.obj="GRangesList"),
                    function(target,GRangesList.obj,strand){
                      gr=as(target,"GRanges")
                      annotate.WithGenicParts(gr,GRangesList.obj,strand)
})


#' function to annotate given GRanges object with promoter,exon,intron & intergenic values
#'  
#' @param target    a granges object storing chromosome locations to be annotated
#' @param feature   a granges object storing chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
#' @param flank     a granges object storing chromosome locations of the flanks of the feature
#' @param feature.name     string for the name of the feature
#' @param flank.name     string for the name o f the flanks
#' @param strand   If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
#' @usage annotate.WithFeature.Flank(target,feature,flank,feature.name="feat",flank.name="flank",strand=FALSE)
#' @return returns an \code{annotationByFeature} object
#' 
#' @export
#' @docType methods
#' @rdname annotate.WithFeature.Flank-methods
setGeneric("annotate.WithFeature.Flank", function(target,feature,flank,feature.name="feat",flank.name="flank",strand=FALSE) standardGeneric("annotate.WithFeature.Flank") )

#' @aliases annotate.WithFeature.Flank,GRanges,GRanges,GRanges-method
#' @rdname annotate.WithFeature.Flank-methods
setMethod( "annotate.WithFeature.Flank", signature(target = "GRanges",feature="GRanges",flank="GRanges"),
                    function(target, feature, flank,feature.name,flank.name,strand){
                      
                      if( ! strand){strand(target)="*"}
                      memb=data.frame(matrix(rep(0,length(target)*2),ncol=2) )  ;colnames(memb)=c(feature.name,flank.name)
                      memb[countOverlaps(target,feature)>0,1]=1
                      memb[countOverlaps(target,flank)>0,2]=1
                    
                      annotation=c(100*sum(memb[,1]>0)/nrow(memb) ,
                                   100*sum(memb[,2]>0)/nrow(memb) ,
                                   100*sum(rowSums(memb)==0)/nrow(memb) )
                      names(annotation)=c(feature.name,flank.name,"other")

                      num.annotation=c( sum(memb[,1]>0) , sum(memb[,2]>0) , sum(rowSums(memb)==0) )
                      names(num.annotation)=c(feature.name,flank.name,"other")                      
                      
                      
                      precedence=c(100*sum(memb[,1]>0)/nrow(memb) ,
                                     100*sum(memb[,2]>0 & memb[,1]==0)/nrow(memb) ,
                                     100*sum(rowSums(memb)==0)/nrow(memb) )
                      names(precedence)=c(feature.name,flank.name,"other")
                      
                      num.precedence=c( sum(memb[,1]>0)  , sum(memb[,2]>0 & memb[,1]==0) , sum(rowSums(memb)==0)  )
                      names(num.precedence)=c(feature.name,flank.name,"other")                      
                      
                      numberOfOlapFeat=c(sum(countOverlaps(feature,target)>0),
                                         sum(countOverlaps(flank,target)>0) )
                      names(numberOfOlapFeat)=c(feature.name,flank.name)
                      percOfOlapFeat =100*numberOfOlapFeat/c(length(feature),length(flank) )
                    
                      #return(list(members=memb,annotation=annotation,precedence=precedence,
                      #            numberOfOlapFeat=numberOfOlapFeat,percOfOlapFeat=percOfOlapFeat) )                      
                      new("annotationByFeature",
                          members         =as.matrix(memb),
                          annotation      =annotation,
                          precedence    =precedence,
                          num.annotation  =num.annotation,
                          num.precedence=num.precedence,
                          no.of.OlapFeat  =numberOfOlapFeat,
                          perc.of.OlapFeat=percOfOlapFeat)
                      
})

#' @aliases annotate.WithFeature.Flank,methylDiff,GRanges,GRanges-method
#' @rdname annotate.WithFeature.Flank-methods
setMethod("annotate.WithFeature.Flank", signature(target= "methylDiff",feature="GRanges",flank="GRanges"),
                    function(target, feature, flank,feature.name,flank.name,strand){
                      gr=as(target,"GRanges")
                      annotate.WithFeature.Flank(gr,feature, flank,feature.name,flank.name,strand)
})


#' function to annotate given GRanges object with a given genomic feature
#' 
#' @param target   a GRanges/or methylDiff object storing chromosome locations to be annotated
#' @param feature  a GRanges object storing chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
#' @param strand   If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
#' @param extend   specifiying a positive value will extend the feature on both sides as much as \code{extend}
#' @param feature.name name of the annotation feature. For example: H3K4me1,CpGisland etc.
#' @usage annotate.WithFeature(target,feature,strand=FALSE,extend=0,feature.name="feat1")
#' @return returns an \code{annotationByFeature} object
#' 
#' @export
#' @docType methods
#' @rdname annotate.WithFeature-methods
setGeneric("annotate.WithFeature", function(target,feature,strand=FALSE,extend=0,feature.name="feat1") standardGeneric("annotate.WithFeature") )

#' @aliases annotate.WithFeature,GRanges,GRanges-method
#' @rdname annotate.WithFeature-methods
setMethod("annotate.WithFeature", signature(target = "GRanges",feature="GRanges"),
                    function(target, feature, strand,extend,feature.name){
 

                      if( ! strand){strand(target)="*"}
                      memb=rep(0,length(target))
                    
                    
                      if(extend>0)
                      {
                        start(feature)=start(feature)-extend
                        end(feature)  =end(feature)+extend
                      }
                      memb[countOverlaps(target,feature)>0]=1
                      
                      annotation=c(100*sum(memb>0)/length(memb) ,
                                   100*sum((memb)==0)/length(memb) )
                      num.annotation=c(sum(memb>0) ,sum((memb)==0)  )
                      names(annotation)=c(feature.name,"other")
                      
                      numberOfOlapFeat=c(sum(countOverlaps(feature,target)>0))
                      percOfOlapFeat =100*numberOfOlapFeat/c(length(feature))
                    
                      new("annotationByFeature",
                          members         =as.matrix(memb),
                          annotation      =annotation,
                          precedence    =annotation,
                          num.annotation  =num.annotation,
                          num.precedence=num.annotation,
                          no.of.OlapFeat  =numberOfOlapFeat,
                          perc.of.OlapFeat=percOfOlapFeat)



})

#' @aliases annotate.WithFeature,methylDiff,GRanges-method
#' @rdname annotate.WithFeature-methods
setMethod("annotate.WithFeature", signature(target = "methylDiff",feature="GRanges"),
                    function(target, feature, strand,extend,feature.name){                      
                      gr=as(target,"GRanges")
                      annotate.WithFeature(gr, feature, strand,extend,feature.name)
})

# ACCESSOR FUNCTIONS
#annotationByFeature
#annotationBygenicparts

                      
#' Get the membership slot of annotationByFeature
#'
#' Membership slot defines the overlap of target features with annotation features
#'  For example, if a target feature overlaps with an exon
#' 
#' @param x a \code{annotationByFeature}  object
#' 
#' @return RETURNS a matrix showing overlap of target features with annotation features. 1 for overlap, 0 for non-overlap
#' 
#' @usage getMembers(x)
#'
#' @export
#' @docType methods
#' @rdname getMembers-methods                      
setGeneric("getMembers", def=function(x) standardGeneric("getMembers"))

#' @aliases getMembers,annotationByFeature-method
#' @rdname getMembers-methods
setMethod("getMembers", signature(x = "annotationByFeature"),
                    function(x){
                      return(x@members)
                      
})


#' Get the percentage of target features overlapping with annotation from annotationByFeature
#'
#' This function retrieves percentage/number of target features overlapping with annotation
#'  
#' @param x a \code{annotationByFeature}  object
#' @param percentage TRUE|FALSE. If TRUE percentage of target features will be returned. If FALSE, number of target features will be returned
#' @param precedence TRUE|FALSE. If TRUE there will be a hierachy of annotation features when calculating numbers (with promoter>exon>intron precedence)
#' That means if a feature overlaps with a promoter it will be counted as promoter overlapping only, or if it is overlapping with a an exon but not a promoter, 
#' it will be counted as exon overlapping only whether or not it overlaps with an intron.
#'
#' @usage getTargetAnnotationStats(x,percentage=TRUE,precedence=TRUE)
#'
#' @return RETURNS  a vector of percentages or counts showing quantity of target features overlapping with annotation
#' 
#' @export
#' @docType methods
#' @rdname getTargetAnnotationStats-methods
setGeneric("getTargetAnnotationStats", def=function(x,percentage=TRUE,precedence=TRUE) standardGeneric("getTargetAnnotationStats"))

#' @rdname getTargetAnnotationStats-methods
#' @aliases getTargetAnnotationStats,annotationByFeature-method
setMethod("getTargetAnnotationStats", signature(x = "annotationByFeature"),
                    function(x,percentage ,precedence ){                      
                      if(percentage){
                        if(precedence){return(x@precedence)
                        }else{return(x@annotation)}
                      }else{
                        if(precedence){return(x@num.precedence)
                        }else{return(x@num.annotation)}                        
                      }

})

#' Get the percentage/count of annotation features overlapping with target features from annotationByFeature
#'
#' This function retrieves percentage/number of annotation features overlapping with targets. 
#' For example, if \code{annotationByFeature}  object is containing statistics of differentially methylated 
#' regions overlapping with gene annotation. This function will return number/percentage of introns,exons and promoters
#' overlapping with differentially methylated regions.
#'  
#' @param x a \code{annotationByFeature}  object
#' @param percentage TRUE|FALSE. If TRUE percentage of annotation features will be returned. If FALSE, number of annotation features will be returned
#'
#' @return RETURNS  a vector of percentages or counts showing quantity of annotation features overlapping with target features
#' @usage getFeatsWithTargetsStats(x,percentage=TRUE)
#' @export
#' @docType methods
#' @rdname getFeatsWithTargetsStats-methods
setGeneric("getFeatsWithTargetsStats", def=function(x,percentage=TRUE) standardGeneric("getFeatsWithTargetsStats"))


#' @rdname getFeatsWithTargetsStats-methods
#' @aliases getFeatsWithTargetsStats,annotationByFeature-method
setMethod("getFeatsWithTargetsStats", signature(x = "annotationByFeature" ),
                    function( x,percentage ){                      
                      if(percentage){
                        return(x@perc.of.OlapFeat)
                      }else{
                        return(x@no.of.OlapFeat)                        
                      }
})


#' Get distance to nearest TSS and gene id from annotationByGenicParts
#'
#' This accessor function gets the nearest TSS, its distance to target feature,strand and name of TSS/gene from annotationByGenicParts object
#' @param x a \code{annotationByGenicParts}  object
#' 
#' @return RETURNS a data.frame containing row number of the target features,distance of target to nearest TSS, TSS/Gene name, TSS strand
#' @usage getAssociationWithTSS(x)
#' @aliases getAssociationWithTSS,-methods getAssociationWithTSS,annotationByGenicParts-method
#' @export
#' @docType methods
#' @rdname annotationByGenicParts-methods
setGeneric("getAssociationWithTSS", def=function(x) standardGeneric("getAssociationWithTSS"))

#' @rdname annotationByGenicParts-methods
#' @docType methods
#' @aliases getAssociationWithTSS,annotationByGenicParts-method
setMethod("getAssociationWithTSS", signature(x = "annotationByGenicParts"),
                    function(x){
                      return(x@dist.to.TSS)
                      
})

# PLOTTING FUNCTIONS

#' Plot annotation categories from annotationByGenicParts or annotationByFeature
#'
#' This function plots a pie or bar chart for showing percentages of targets annotated by genic parts or other query features
#' @param x a \code{annotationByFeature} or  \code{annotationByGenicParts} object
#' @param precedence TRUE|FALSE. If TRUE there will be a hierachy of annotation features when calculating numbers (with promoter>exon>intron precedence). 
#'  This option is only valid when x is a \code{annotationByGenicParts} object
#' @param col a vector of colors for piechart or the bar plot
#' @param ... graphical parameters to be passed to \code{pie} or \code{barplot} functions
#'
#' usage  plotTargetAnnotation(x,precedence=TRUE,col,...)
#'
#' @return plots a piechart or a barplot for percentage of the target features overlapping with annotation
#' 
#' @export
#' @docType methods
#' @rdname plotTargetAnnotation-methods
setGeneric("plotTargetAnnotation", def=function(x,precedence=TRUE,col=rainbow(length(x@annotation)),...) standardGeneric("plotTargetAnnotation"))

#' @rdname plotTargetAnnotation-methods
#' @docType methods
#' @aliases plotTargetAnnotation,annotationByFeature-method
setMethod("plotTargetAnnotation", signature(x = "annotationByFeature"),
                    function(x,precedence,col,...){
                      props=getTargetAnnotationStats(x,precedence)

                      if(precedence | ( is(x,"annotationByFeature") & !is(x,"annotationByGenicParts")) ){
                        slice.names=names(props)
                        #names(props)=paste(names(props),paste(round(props),"%"),sep=" ")
                        names(props)=paste( paste(round(props),"%"),sep=" ")

                        pie(props,cex=0.9,col=col,...)
                         legend("topright",legend=slice.names,fill=col )

                      }
                      else{
                        
                        mids=barplot(props,col=col,...)  
                        text(mids,y=props,labels=paste(paste(round(props),"%",sep="")),pos=1) 
                      }

})



# SECTION 3:
# annotate ML objects with annotations read-in and converted to GRanges objects

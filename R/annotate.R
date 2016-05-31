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
  #colClasses=c("factor","integer","integer","factor","integer","factor",
  #"integer","integer","integer","integer","factor","factor") 
  #ref         =read.table(bed.file,header=F,skip=skip,
  #colClasses = colClasses,stringsAsFactors=F) # give colClasses for fast read-in
  ref=unique(ref)
  # apply strsplit on columns 11 and 12 to extract exon 
  # start positions and exon sizes
  b.start.size=cbind(as.integer(unlist(strsplit(as.character(ref$V12),"," ))),
                     as.integer(unlist(strsplit(as.character(ref$V11),"," ))))
  
  # replicate rows occurs as many as its exon number
  rep.ref     =ref[rep(1:nrow(ref),ref[,10]),] 
  #exon.id     =unlist(sapply(ref[,10],function(x) 1:x));rep.ref$V5=exon.id
  exon.id     =unlist( mapply( function(x,y) if(x=="+"){
    return(1:y)}else{return(y:1)} ,ref[,6],ref[,10]  ) )
  rep.ref$V5=exon.id

  # calculate exon start and ends
  rep.ref$V3  =rep.ref$V2+b.start.size[,1]+b.start.size[,2] 
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
  #colClasses=c("factor","integer","integer","factor","integer","factor",
  # "integer","integer","integer","integer","factor","factor") 
  #ref         =read.table(bed.file,header=F,skip=skip,
  # colClasses = colClasses,stringsAsFactors=F) # give colClasses for fast read-in

  #remove the genes with one exon only (they won't have any introns)
  ref=ref[ref[,10]>1,]
  ids=paste(ref[,1],ref[,2],ref[,3],ref[,4],sep="")
  ref=cbind(ref,id=ids)
  ref=unique(ref)
  # apply strsplit on columns 11 and 12 to extract 
  # exon start positions and exon sizes
  b.start.size=cbind(as.integer(unlist(strsplit(as.character(ref$V12),"," ))),
                     as.integer(unlist(strsplit(as.character(ref$V11),"," ))))
  
  # replicate rows occurs as many as its exon number
  rep.ref     =ref[rep(1:nrow(ref),ref[,10]),] 
  #exon.id     =unlist(sapply(ref[,10],function(x) 1:x));rep.ref$V5=exon.id
  exon.id     =unlist( mapply( function(x,y) if(x=="+"){
    return(1:y)}else{return(y:1)} ,ref[,6],ref[,10]  ) );rep.ref$V5=exon.id
  
  # calculate exon start and ends
  rep.ref$V3  =rep.ref$V2+b.start.size[,1]+b.start.size[,2] 
  rep.ref$V2  =rep.ref$V2+b.start.size[,1]
  rep.ref     =rep.ref[,c(1:6,13)]
  
  # now try to get the exons by cbinding consecutive exons
  temp.ref    =cbind(rep.ref[1:(nrow(rep.ref)-1),],rep.ref[2:nrow(rep.ref),])
  temp.ref    =temp.ref[temp.ref[,7]==temp.ref[,14],]
  temp.ref[,2]=temp.ref[,3]
  temp.ref[,3]=temp.ref[,9]
  rep.ref     =temp.ref[,1:6]
  
  # subtract 1 from - strand exon numbers
  rep.ref[rep.ref[,6]=="-",5]=rep.ref[rep.ref[,6]=="-",5]-1 
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
  #chr= sum(grepl("chr",bed.df[,1]) )

  if(type=="exons" | type=="introns")
  {
    #does it have 12>= columns
    ex=(ncol(bed.df)>=12) 
    
    return(num.col & col1.2  & ex)
  }
  return(num.col & col1.2  )

}

#######################################
# SECTION 1: S4 functions
#######################################

# convert a data frame read-in from a bed file to a GRanges object
#  
# @param bed  a data.frame where column order and content resembles 
# a bed file with 12 columns
# @usage convert.bed.df(bed)
# @return \code{\link{GRanges}} object
#
# @note one bed track per file is only accepted, the bed files with 
# multiple tracks will cause en error
#
#  @export
#  @docType methods
#  @rdname convert.bed.df-methods
setGeneric("convert.bed.df",function(bed) standardGeneric("convert.bed.df"))

#  @aliases convert.bed.df,data.frame-method
#  @rdname convert.bed.df-methods
setMethod("convert.bed.df" ,signature(bed = "data.frame" ),
                            function(bed){
                              
  if(! check.bed.validity(bed)){stop("this is not a valid bed file")}
  if(ncol(bed)>=6){
    # check if 6th column is really strand
    if( sum( unique(bed[,6]) %in% c("+","-",".") ) != length(unique(bed[,6])) ){
      stop("Strand column of the bed file or data.frame is wrong")}
    #convert to granges
    strands=as.character(bed$V6)
    strands[strands=="."]="*"
    grange=GRanges(seqnames=as.character(bed$V1),
                  ranges=IRanges(start=bed$V2+1, end=bed$V3),
                  strand=strands, score=bed$V5,name=bed$V4)
  }
  if(ncol(bed)==5){
    grange=GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2+1, end=bed$V3),
                   strand=rep("*",nrow(bed)), score=bed$V5,name=bed$V4 )
  }
  if(ncol(bed)==4){
    grange=GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2+1, end=bed$V3),
                   strand=rep("*",nrow(bed)),name=bed$V4)
  }    
  if(ncol(bed)==3){
    grange=GRanges(seqnames=bed$V1,ranges=IRanges(start=bed$V2+1, end=bed$V3),
                   strand=rep("*",nrow(bed)))
  }                           
  return(grange)
})

# convert a data frame read-in from a bed file to a GRanges object for exons
#  
# @param bed.df  a data.frame where column order and content resembles a 
# bed file with 12 columns
# @usage convert.bed2exons(bed.df)
# @return \code{\link{GRanges}} object
#
# @note one bed track per file is only accepted, the bed files with multiple 
# tracks will cause en error
#
#  @export 
#  @docType methods
#  @rdname convert.bed2exons-methods
setGeneric("convert.bed2exons",function(bed.df) 
  standardGeneric("convert.bed2exons"))

#  @aliases convert.bed2exons,data.frame-method
#  @rdname convert.bed2exons-methods
setMethod("convert.bed2exons" ,signature(bed.df = "data.frame" ),
                            function(bed.df){
      
    if(! check.bed.validity(bed.df,"exon")){
      stop("this is not a valid bed file")
    }
    bed12.to.exons(bed.df)
})

# convert a data frame read-in from a bed file to a GRanges object for introns
#  
# @param bed.df  a data.frame where column order and content resembles a bed 
# file with 12 columns
# @usage convert.bed2introns(bed.df)
# @return \code{\link{GRanges}} object
#
# @note one bed track per file is only accepted, the bed files with multiple 
# tracks will cause en error
#
#  @export
#  @docType methods
# @rdname convert.bed2introns-methods
setGeneric("convert.bed2introns",function(bed.df) 
  standardGeneric("convert.bed2introns"))

#  @aliases convert.bed2introns,data.frame-method
#  @rdname convert.bed2introns-methods
setMethod("convert.bed2introns" ,signature(bed.df = "data.frame" ),
                            function(bed.df){
                              
  if(! check.bed.validity(bed.df,"exon")){stop("this is not a valid bed file")}                      
  bed12.to.introns(bed.df)
})


#' read a bed file and convert it to GRanges
#' 
#' The function reads a BED file and coverts it to a 
#' \code{\link[GenomicRanges]{GRanges}} object  
#'    
#' @param location  location of the file, a character string 
#' such as: "/home/user/my.bed"
#' @param remove.unsual if TRUE(default) remove the chromomesomes with unsual 
#' names, mainly random chromsomes etc
#'
#' @usage readBed(location,remove.unsual=TRUE)
#' @return \code{\link{GRanges}} object
#' @examples
#' bed.file=system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit")
#' bed.gr=readBed(location=bed.file,remove.unsual=TRUE)
#' 
#' @note one bed track per file is only accepted, the bed files with multiple 
#' tracks will cause an error
#'
#' @export
#' @docType methods
#' @rdname readBed-methods
setGeneric("readBed", function(location,remove.unsual=TRUE) 
  standardGeneric("readBed"))

#' @aliases readBed,character-method read.bed
#' @rdname readBed-methods
setMethod("readBed", signature(location = "character"),#remove.unsual="logical" ),
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

#' Read  transcript features from a BED file
#' 
#' The function returns a \code{\link[GenomicRanges]{GRangesList}} containing 
#' exon, intron, TSS(transcription start site) and promoter locations
#'
#' @param location location of the bed file with 12 or more columns
#' @param remove.unsual remove the chromomesomes with unsual names, mainly 
#'        random chromsomes etc
#' @param up.flank  up-stream from TSS to detect promoter boundaries
#' @param down.flank down-stream from TSS to detect promoter boundaries
#' @param unique.prom     get only the unique promoters, promoter boundaries 
#'        will not have a gene name if you set this option to be TRUE
#' @return a \code{\link[GenomicRanges]{GRangesList}} containing locations of 
#' exon/intron/promoter/TSS
#' @examples
#' gene.obj=readTranscriptFeatures(system.file("extdata", 
#' "refseq.hg18.bed.txt", package = "methylKit"))
#' 
#' @note  one bed track per file is only accepted, the bed files with 
#' multiple tracks will cause en error
#'
#' @export
#' @docType methods
#' @rdname readTranscriptFeatures-methods
setGeneric("readTranscriptFeatures", function(location,remove.unsual=TRUE,
                                                up.flank=1000,down.flank=1000,
                                                unique.prom=TRUE) 
  standardGeneric("readTranscriptFeatures"))

#' @aliases readTranscriptFeatures,character-method read.transcript.features
#' @rdname readTranscriptFeatures-methods
setMethod("readTranscriptFeatures", signature(location = "character"),
                    function(location,remove.unsual,up.flank ,down.flank ,
                             unique.prom){
      
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
            strand=as.character(prom.df$V6),score=rep(0,nrow(prom.df)),
            name=prom.df$V4)
    }else{
      prom.df=unique(bed[,c(1,2,3,6)])
      prom=GRanges(seqnames=as.character(prom.df$V1),
            ranges=IRanges(start=prom.df$V2, end=prom.df$V3),
            strand=as.character(prom.df$V6),score=rep(0,nrow(prom.df)),
            name=rep(".",nrow(prom.df)) )
      
    }
  
    GRangesList(exons=exons,introns=introns,promoters=prom,TSSes=tssg)
})

#' Get upstream and downstream adjacent regions to a genomic feature
#' 
#' The function returns flanking regions on either side of a genomic feature. 
#' It is useful for getting flanking regions such as CpG island shores.
#' 
#' @param grange \code{\link[GenomicRanges]{GRanges}} object for the feature
#' @param flank  number of basepairs for the flanking regions
#' @param clean  If set to TRUE, flanks overlapping with other main features 
#'               will be trimmed, and overlapping flanks will be removed
#'        this will remove multiple counts when other features overlap with flanks
#'
#' @examples
#' 
#' # read the bed file as GRanges object
#' bed.file=system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit")
#' bed.gr=readBed(location=bed.file,remove.unsual=TRUE)
#' 
#' # get flanks on the either side
#' bed.flanks=getFlanks(bed.gr,flank=2000,clean=TRUE)
#' 
#' @return \code{\link[GenomicRanges]{GRanges}} object for flanking regions
#' @export
#' @docType methods
#' @rdname getFlanks-methods
setGeneric("getFlanks", function(grange,flank=2000,clean=TRUE) 
  standardGeneric("getFlanks"))

#' @aliases getFlanks,GRanges-method
#' @rdname getFlanks-methods
setMethod("getFlanks", signature(grange= "GRanges"),
                    function(grange,flank=2000,clean=T){
          
    shores=c( IRanges::flank(grange,flank),IRanges::flank(grange,flank,FALSE) )
    if(clean){
      # merge overlapping shores remove CpG coordinates from all shores
      shores=IRanges::reduce(IRanges::setdiff(shores, grange)) 
    }
    shores
})

#' a function to read-in genomic features and their upstream and downstream 
#' adjecent regions such as CpG islands and their shores
#'
#' @param location for the bed file of the feature 
#' @param flank    number of basepairs for the flanking regions
#' @param clean    If set to TRUE, flanks overlapping with other main features 
#'                 will be trimmed
#' @param remove.unsual  remove chromsomes with unsual names random, Un and 
#' antyhing with "_" character
#' @param feature.flank.name the names for feature and flank ranges, it should 
#' be a character vector of length 2. example: c("CpGi","shores")
#' @return a \code{\link[GenomicRanges]{GenomicRangesList}} contatining one 
#' GRanges object for flanks and one for GRanges object for the main feature.
#'   
#' @examples
#'  # location of the example CpG file
#'  my.loc=system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit")
#'  cpg.obj=readFeatureFlank(location=my.loc,
#'  feature.flank.name=c("CpGi","shores"))
#'
#' @export
#' @docType methods
#' @rdname readFeatureFlank-methods
setGeneric("readFeatureFlank", function(location,remove.unsual=TRUE,
                                          flank=2000,clean=TRUE,
                                          feature.flank.name=NULL) 
  standardGeneric("readFeatureFlank") )

#' @aliases readFeatureFlank,character-method read.feature.flank
#' @rdname readFeatureFlank-methods
setMethod("readFeatureFlank", signature(location = "character"),
                    function(location,remove.unsual,flank ,
                             clean,feature.flank.name){
      feat=readBed(location,remove.unsual)
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
# annotate granges objects with annotations that read-in and converted to 
# GRanges objects
##############################################################################

#######################################
# SECTION 2: Define new classes
#######################################

# A set of objects that will hold statistics about feature and annotation overlap

#' An S4 class for  overlap of target features with a generic annotation 
#'
#' This object is desgined to hold statistics and information about genomic 
#' feature overlaps it extends \code{\link{list}} class.
#'          
#' @section Slots:\describe{
#'                  \item{\code{members}}{a matrix showing overlap of target 
#'                  features with annotation genomic features}
#'                  \item{\code{annotation}}{a named vector of percentages of 
#'                  overlap between feature and annotation}'
#'                  \item{\code{precedence}}{a named vector of percentages of 
#'                  overlap between feature and annotation}
#'                  \item{\code{num.annotation}}{a named vector of numbers of 
#'                  overlap between feature and annotation}
#'                  \item{\code{num.precedence}}{a named vector of numbers of 
#'                  overlap between feature and annotation}
#'                  \item{\code{no.of.OlapFeat}}{vector}
#'                  \item{\code{perc.of.OlapFeat}}{vector}
#' }
#' @examples
#' data(methylKit)
#' cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
#'                             package = "methylKit"),
#'                             feature.flank.name=c("CpGi","shores"))
#' 
#' # the following function returns annotationByFeature object
#' ann=annotateWithFeatureFlank(methylDiff.obj,cpg.obj$CpGi,cpg.obj$shores,
#'                                feature.name="CpGi",flank.name="Shores")
#' ann
#' 
#' @seealso see \code{\link[methylKit]{annotateWithFeatureFlank}} and 
#' \code{\link[methylKit]{annotateWithFeature}} on how to create this object.
#'          see following functions that operates on this object:
#'            \code{\link[methylKit]{getMembers}}, 
#'            \code{\link[methylKit]{getTargetAnnotationStats}},
#'            \code{\link[methylKit]{getFeatsWithTargetsStats}}, 
#'            \code{\link[methylKit]{plotTargetAnnotation}}
#' @name annotationByFeature-class
#' @aliases annotationByFeature
#' @docType class
#' @rdname annotationByFeature-class
#' @export
setClass("annotationByFeature", representation(members         ="matrix",
                                      annotation      ="numeric",
                                      precedence    ="numeric",
                                      num.annotation  ="numeric",
                                      num.precedence="numeric",
                                      no.of.OlapFeat  ="numeric",
                                      perc.of.OlapFeat="numeric"))

#' An S4 class for overlap of target features with gene annotation 
#'
#' This object is desgined to hold basic statistics and information about 
#' genomic feature overlaps. It extends 
#' \code{\link[methylKit]{annotationByFeature}} class.
#' The class contains an additional slot for data containing distance to nearest 
#' transcription start site (TSS).
#'          
#' @section Slots:\describe{
#'                  \item{\code{members}}{a matrix showing overlap of target 
#'                  features with gene annotation  features}
#'                  \item{\code{annotation}}{a named vector of percentages of 
#'                  overlap between feature and gene annotation}
#'                  \item{\code{precedence}}{a named vector of percentages of 
#'                  overlap between feature and  gene annotation}
#'                  \item{\code{num.annotation}}{a named vector of numbers of 
#'                  overlap between feature and gene annotation}
#'                  \item{\code{num.precedence}}{a named vector of numbers of 
#'                  overlap between feature and gene annotation}
#'                  \item{\code{no.of.OlapFeat}}{a named vector of numbers of 
#'                  overlap between gene annotation and the feature}
#'                  \item{\code{perc.of.OlapFeat}}{a named vector of percentages 
#'                  of overlap between gene annotation and the feature }
#'                  \item{\code{dist.to.TSS}}{a data frame showing distances to 
#'                  TSS and gene/TSS names and strand}
#' }
#' @examples
#' data(methylKit)
#' gene.obj=readTranscriptFeatures(system.file("extdata", 
#'          "refseq.hg18.bed.txt", package = "methylKit"))
#' 
#' # the following function returns an annotationByGenicParts object
#' ann=annotateWithGenicParts(methylDiff.obj, gene.obj)
#' 
#' @seealso see \code{\link[methylKit]{annotateWithGenicParts}} on 
#' how to create this object, see following functions that operates on this object  
#' \code{\link[methylKit]{getTargetAnnotationStats}}, 
#' \code{\link[methylKit]{getMembers}}, 
#' \code{\link[methylKit]{getAssociationWithTSS}},
#' \code{\link[methylKit]{getTargetAnnotationStats}}, 
#' code{\link[methylKit]{getFeatsWithTargetsStats}},
#' \code{\link[methylKit]{plotTargetAnnotation}}
#' 
#' @name annotationByGenicParts-class
#' @aliases annotationByGenicParts
#' @docType class
#' @rdname annotationByGenicParts-class
#' @export
setClass("annotationByGenicParts", representation(dist.to.TSS   ="data.frame"),
         contains="annotationByFeature")


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
  
  cat("summary of target set annotation with genic parts\n")
  cat(nrow(object@members));cat(" rows in target set\n--------------\n")
  cat("--------------\n")
  cat("percentage of target features overlapping with annotation :\n")
  print(object@annotation);cat("\n\n")
  cat("percentage of target features overlapping with annotation ",
      "(with promoter>exon>intron precedence) :\n")
  print(object@precedence) ;cat("\n\n")
  cat("percentage of annotation boundaries with feature overlap :\n")
  print(object@perc.of.OlapFeat);cat("\n\n")  
  cat("summary of distances to the nearest TSS :\n")
  print(summary(abs(object@dist.to.TSS[,2])))
})


#' @rdname show-methods
#' @aliases show,annotationByFeature-method
setMethod("show", "annotationByFeature", function(object) {
  
  cat("summary of target set annotation with feature annotation\n")
  cat(nrow(object@members));cat(" rows in target set\n--------------\n")
  cat("--------------\n")
  cat("percentage of target features overlapping with annotation :\n")
  print(object@annotation);cat("\n\n")
  cat("percentage of annotation boundaries with feature overlap :\n")
  print(object@perc.of.OlapFeat);cat("\n\n")  
})

#######################################
# SECTION 2: S3 FUNCTIONS 
# these shouldn't be exported
#######################################

annotate.gr.WithGenicParts<-function(gr,prom,exon,intron,strand=F)
{
  #require(GenomicRanges)
  if( ! strand){strand(gr)="*"}
  memb=data.frame(matrix(rep(0,length(gr)*3),ncol=3) )  
  colnames(memb)=c("prom","exon","intron")
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
  return(list(members=memb,annotation=annotation,precedence=precedence,
              num.annotation=num.annotation,
              num.precedence=num.precedence,
              numberOfOlapFeat=numberOfOlapFeat,percOfOlapFeat=percOfOlapFeat) )
  
}

distance2nearestFeature<-function(g.idh,tss)
{
  #require(GenomicRanges)
  elementMetadata(g.idh)=DataFrame(elementMetadata(g.idh),
                                   orig.row=1:length(g.idh))
  id.col    =ncol( elementMetadata(g.idh))+5 # get the row number column

  # get the id column for tss in the merged data.frame below
  tss.id.col=ncol( elementMetadata(g.idh))+5+7 
  strnd.col=ncol( elementMetadata(g.idh))+5+7-2 # get the id column for tss
 
  # a merged data.frame is returned by this function
  met.tss   = .nearest.2bed(g.idh, tss) 

  dist.col=ncol(met.tss)
  cond=met.tss$end<met.tss$start.y & met.tss$strand.y=="+"
  met.tss[cond,dist.col] = -1* met.tss[cond,dist.col]

  cond2=met.tss$end>met.tss$start.y & met.tss$strand.y=="-"
  met.tss[cond2,dist.col] = -1* met.tss[cond2,dist.col]

  res=met.tss[order(met.tss[,id.col]),c(id.col,dist.col,tss.id.col,strnd.col)]
  names(res)=c("target.row","dist.to.feature","feature.name","feature.strand")

  return(res)
}


# get the nearest features between a subject and query bed file
# get grange objects in bed zormat
# g.bed is GenomicRanges object, biatch!!!
# subject is also GenomicRanges object
.nearest.2bed<-function(g.bed,subject)
{
  chrs1=S4Vectors::levels(seqnames(g.bed))
  chrs2=S4Vectors::levels(seqnames(subject))
  chrs=chrs1[chrs1 %in% chrs2]
  res.df=NULL
  for(i in 1:length(chrs))
  {
    # find the nearest for this range
    ind  =nearest(ranges(g.bed[seqnames(g.bed)==chrs[i],]),
                  ranges(subject[seqnames(subject)==chrs[i],]))
    sbdf1 =S4Vectors::as.data.frame(g.bed[seqnames(g.bed)==chrs[i],] )
    sbdf2 =S4Vectors::as.data.frame(subject[seqnames(subject)==chrs[i],] )
    sbdf2 =sbdf2[ind,]
    names(sbdf2)=paste(names(sbdf2),".y",sep="")
    temp  =cbind(sbdf1,sbdf2)
    res.df=rbind(res.df,temp)
  }

  res.dist=rep(0,nrow(res.df))
  dist1=abs(res.df$start-res.df$end.y)+1
  dist2=abs(res.df$start.y-res.df$end)+1
  totti=res.df$width + res.df$width.y #total length
  res.dist[pmax(dist1,dist2)>totti]=pmin(dist1,dist2)[pmax(dist1,dist2) > totti] 
  res.df=cbind(res.df,dist=res.dist)
  return(res.df)
}





#######################################
# SECTION 2: S4 FUNCTIONS 
#######################################


#' Annotate given object with gene annotation
#'
#' The function annotates given genomic feature or methylKit object with gene 
#' annotation such as promoter, exon, intron & intergenic.
#' It also gets the distance to nearest TSS (transcription start site) for each 
#' genomic feature or methylation event.
#' 
#'
#' @param target  a \code{\link[methylKit]{methylDiff}} or a 
#' \code{\link[GenomicRanges]{GRanges}} object storing chromosome locations to 
#' be annotated
#' @param GRangesList.obj   A \code{\link[GenomicRanges]{GRangesList}} object 
#' containing GRanges object for promoter,exons,introns and TSSes, or simply 
#' output of \code{\link[methylKit]{readTranscriptFeatures}} function
#' @param strand If set to TRUE, annotation features and target features will be 
#' overlapped based on strand  (def:FALSE)
#' @usage annotateWithGenicParts(target,GRangesList.obj,strand=FALSE)
#' @return \code{\link[methylKit]{annotationByGenicParts}} object
#' @examples
#' data(methylKit)
#' gene.obj=readTranscriptFeatures(system.file("extdata", 
#' "refseq.hg18.bed.txt", package = "methylKit"))
#' annotateWithGenicParts(methylDiff.obj, gene.obj)
#' @seealso 
#' \code{\link[methylKit]{getMembers}}, 
#' \code{\link[methylKit]{getTargetAnnotationStats}},
#' \code{\link[methylKit]{getFeatsWithTargetsStats}}, 
#' \code{\link[methylKit]{getAssociationWithTSS}},
#' \code{\link[methylKit]{plotTargetAnnotation}}
#' @export
#' @docType methods
#' @rdname annotateWithGenicParts-methods
setGeneric("annotateWithGenicParts", 
           function(target,GRangesList.obj,strand=FALSE) 
             standardGeneric("annotateWithGenicParts") )

#' @aliases annotateWithGenicParts,GRanges,GRangesList-method annotate.WithGenicParts
#' @rdname annotateWithGenicParts-methods
setMethod("annotateWithGenicParts", signature(target= "GRanges",
                                               GRangesList.obj="GRangesList"),
                    function(target,GRangesList.obj,strand){

a.list    =annotate.gr.WithGenicParts(target,GRangesList.obj$promoters,
                                      GRangesList.obj$exons,
                                      GRangesList.obj$introns,strand=strand)
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

#' @aliases annotateWithGenicParts,methylDiff,GRangesList-method
#' @rdname annotateWithGenicParts-methods
setMethod("annotateWithGenicParts", signature(target = "methylDiff",
                                               GRangesList.obj="GRangesList"),
                    function(target,GRangesList.obj,strand){
    gr=as(target,"GRanges")
    annotateWithGenicParts(gr,GRangesList.obj,strand)
})


#' Annotate an object with two sets of genomic features
#'
#' The function annotates given genomic feature or methylKit object with two 
#' sets of annotation. 
#' It is primarily useful when annotating objects with CpG islands and shores.
#' 
#'  
#' @param target    a \code{\link[methylKit]{methylDiff}} or a 
#' \code{\link[GenomicRanges]{GRanges}} object storing chromosome locations to 
#' be annotated
#' @param feature   a \code{\link[GenomicRanges]{GRanges}} object storing 
#' chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
#' @param flank     a \code{\link[GenomicRanges]{GRanges}} object storing 
#' chromosome locations of the flanks of the feature
#' @param feature.name     string for the name of the feature
#' @param flank.name     string for the name o f the flanks
#' @param strand   If set to TRUE, annotation features and target features will 
#' be overlapped based on strand  (def:FALSE)
#' @return returns an \code{\link[methylKit]{annotationByFeature}} object
#' @examples
#' data(methylKit)
#' cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
#' package = "methylKit"),feature.flank.name=c("CpGi","shores"))
#' 
#' annotateWithFeatureFlank(methylDiff.obj,cpg.obj$CpGi,cpg.obj$shores,
#' feature.name="CpGi",flank.name="Shores")
#'
#' @seealso
#' \code{\link[methylKit]{getMembers}}, 
#' \code{\link[methylKit]{getTargetAnnotationStats}},
#' \code{\link[methylKit]{getFeatsWithTargetsStats}}, 
#' \code{\link[methylKit]{plotTargetAnnotation}}
#' 
#' @export
#' 
#' @docType methods
#' @rdname annotateWithFeatureFlank-methods
setGeneric("annotateWithFeatureFlank", function(target,feature,flank,
                                                  feature.name="feat",
                                                  flank.name="flank",
                                                  strand=FALSE) 
  standardGeneric("annotateWithFeatureFlank") )

#' @aliases annotateWithFeatureFlank,GRanges,GRanges,GRanges-method 
#'          annotate.WithFeature.Flank
#' @rdname annotateWithFeatureFlank-methods
setMethod( "annotateWithFeatureFlank", signature(target = "GRanges",
                                                   feature="GRanges",
                                                   flank="GRanges"),
                    function(target, feature, flank,feature.name,
                             flank.name,strand){
                      
    if( ! strand){strand(target)="*"}
    memb=data.frame(matrix(rep(0,length(target)*2),ncol=2) )  ;colnames(memb)=c(feature.name,flank.name)
    memb[countOverlaps(target,feature)>0,1]=1
    memb[countOverlaps(target,flank)>0,2]=1
  
    annotation=c(100*sum(memb[,1]>0)/nrow(memb) ,
                 100*sum(memb[,2]>0)/nrow(memb) ,
                 100*sum(rowSums(memb)==0)/nrow(memb) )
    names(annotation)=c(feature.name,flank.name,"other")
  
    num.annotation=c( sum(memb[,1]>0) , sum(memb[,2]>0) , 
                      sum(rowSums(memb)==0) )
    names(num.annotation)=c(feature.name,flank.name,"other")                      
    
    
    precedence=c(100*sum(memb[,1]>0)/nrow(memb) ,
                   100*sum(memb[,2]>0 & memb[,1]==0)/nrow(memb) ,
                   100*sum(rowSums(memb)==0)/nrow(memb) )
    names(precedence)=c(feature.name,flank.name,"other")
    
    num.precedence=c( sum(memb[,1]>0)  , sum(memb[,2]>0 & memb[,1]==0) , 
                      sum(rowSums(memb)==0)  )
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

#' @aliases annotateWithFeatureFlank,methylDiff,GRanges,GRanges-method
#' @rdname annotateWithFeatureFlank-methods
setMethod("annotateWithFeatureFlank", signature(target= "methylDiff",
                                                  feature="GRanges",
                                                  flank="GRanges"),
                    function(target, feature, flank,feature.name,
                             flank.name,strand){
                      gr=as(target,"GRanges")
                      annotateWithFeatureFlank(gr,feature, flank,feature.name,
                                                 flank.name,strand)
})

 
#' Annotate object with a set of genomic features
#'
#' The function annotates given genomic feature or methylKit object with a set 
#' of annotation. 
#' It is primarily useful when annotating objects with simple genomic features, 
#' such as enhancer locations.
# 
#' @param target    a \code{\link[methylKit]{methylDiff}} or a 
#' \code{\link[GenomicRanges]{GRanges}} object storing chromosome locations to 
#' be annotated
#' @param feature  a \code{\link[GenomicRanges]{GRanges}} object storing 
#' chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
#' @param strand   If set to TRUE, annotation features and target features will 
#' be overlapped based on strand  (def:FALSE)
#' @param extend   specifiying a positive value will extend the feature on 
#' both sides as much as \code{extend}
#' @param feature.name name of the annotation feature. For example: H3K4me1,
#' CpGisland etc.
#' @examples
#' data(methylKit)
#' cpg.gr=readBed(system.file("extdata", "cpgi.hg18.bed.txt", 
#' package = "methylKit"),remove.unsual=TRUE)
#' 
#' annotateWithFeature(methylDiff.obj,cpg.gr,strand=FALSE,extend=0,
#' feature.name="CpGi")
#' 
#' @return returns an \code{\link[methylKit]{annotationByFeature}} object
#' 
#' @seealso
#' \code{\link[methylKit]{getMembers}}, 
#' \code{\link[methylKit]{getTargetAnnotationStats}},
#' \code{\link[methylKit]{getFeatsWithTargetsStats}}, 
#' \code{\link[methylKit]{plotTargetAnnotation}}
#' 
#' @export
#' @docType methods
#' @rdname annotateWithFeature-methods
setGeneric("annotateWithFeature", function(target,feature,strand=FALSE,
                                            extend=0,feature.name="feat1") 
  standardGeneric("annotateWithFeature") )

#' @aliases annotateWithFeature,GRanges,GRanges-method annotate.WithFeature
#' @rdname annotateWithFeature-methods
setMethod("annotateWithFeature", signature(target = "GRanges",feature="GRanges"),
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

#' @aliases annotateWithFeature,methylDiff,GRanges-method
#' @rdname annotateWithFeature-methods
setMethod("annotateWithFeature", signature(target = "methylDiff",feature="GRanges"),
                    function(target, feature, strand,extend,feature.name){                      
  gr=as(target,"GRanges")
  annotateWithFeature(gr, feature, strand,extend,feature.name)
})

# ACCESSOR FUNCTIONS
#annotationByFeature
#annotationBygenicparts

                      
#' Get the membership slot of annotationByFeature
#'
#' Membership slot defines the overlap of target features with annotation 
#' features as a matrix.
#'
#' 
#' @param x an \code{\link[methylKit]{annotationByFeature}}  object
#' 
#' @return a matrix showing overlap of target features with annotation features. 
#' 1 for overlap, 0 for non-overlap. 
#' Each row in the matrix corresponds to a genomic feature that is annoted by 
#' one of the following functions:
#' \code{\link[methylKit]{annotateWithFeature}},
#' \code{\link[methylKit]{annotateWithFeatureFlank}},
#' \code{\link[methylKit]{annotateWithGenicParts}}
#' 
#' @examples
#' data(methylKit)
#' cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
#' package = "methylKit"),feature.flank.name=c("CpGi","shores"))
#' 
#' ann=annotateWithFeatureFlank(methylDiff.obj,cpg.obj$CpGi,cpg.obj$shores,
#' feature.name="CpGi",flank.name="Shores")
#' mat=getMembers(ann)
#' head(mat)
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


#' Get the percentage of target features overlapping with annotation from 
#' annotationByFeature
#'
#' This function retrieves percentage/number of target features overlapping 
#' with annotation
#'  
#' @param x an \code{\link[methylKit]{annotationByFeature}} or an 
#' \code{\link[methylKit]{annotationByGenicParts}} object
#' @param percentage TRUE|FALSE. If TRUE percentage of target features will be 
#' returned. If FALSE, number of target features will be returned
#' @param precedence TRUE|FALSE. If TRUE there will be a hierachy of annotation 
#' features when calculating numbers (with promoter>exon>intron precedence)
#' That means if a feature overlaps with a promoter it will be counted as 
#' promoter overlapping only, or if it is overlapping with a an exon but not 
#' a promoter, 
#' it will be counted as exon overlapping only whether or not it overlaps with 
#' an intron.
#'
#' @examples
#' data(methylKit)
#' cpg.obj=readFeatureFlank(system.file("extdata", 
#' "cpgi.hg18.bed.txt", package = "methylKit"),
#' feature.flank.name=c("CpGi","shores"))
#' 
#' ann=annotateWithFeatureFlank(methylDiff.obj,cpg.obj$CpGi,cpg.obj$shores,
#' feature.name="CpGi",flank.name="Shores")
#' getTargetAnnotationStats(ann,percentage=TRUE,precedence=TRUE)
#' 
#' @return a \code{vector} of percentages or counts showing quantity of 
#' target features overlapping with annotation
#' 
#' @export
#' @docType methods
#' @rdname getTargetAnnotationStats-methods
setGeneric("getTargetAnnotationStats", def=function(x,percentage=TRUE,
                                                    precedence=TRUE)
  standardGeneric("getTargetAnnotationStats"))

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

#' Get the percentage/count of annotation features overlapping with target features
#'
#' This function retrieves percentage/number of annotation features overlapping 
#' with targets. 
#' For example, if \code{annotationByFeature}  object is containing statistics 
#' of differentially methylated 
#' regions overlapping with gene annotation. This function will return 
#' number/percentage of introns,exons and promoters
#' overlapping with differentially methylated regions.
#'  
#' @param x an \code{\link[methylKit]{annotationByFeature}} or an 
#' \code{\link[methylKit]{annotationByGenicParts}} object
#' @param percentage TRUE|FALSE. If TRUE percentage of annotation features will 
#' be returned. If FALSE, number of annotation features will be returned
#'
#' @return a vector of percentages or counts showing quantity of annotation 
#' features overlapping with target features
#' @examples
#' data(methylKit)
#' cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt",
#'  package = "methylKit"),feature.flank.name=c("CpGi","shores"))
#' 
#' ann=annotateWithFeatureFlank(methylDiff.obj,cpg.obj$CpGi,cpg.obj$shores,
#' feature.name="CpGi",flank.name="Shores")
#' getFeatsWithTargetsStats(ann,percentage=TRUE)
#' 
#' @export
#' @docType methods
#' @rdname getFeatsWithTargetsStats-methods
setGeneric("getFeatsWithTargetsStats", def=function(x,percentage=TRUE) 
  standardGeneric("getFeatsWithTargetsStats"))


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
#' This accessor function gets the nearest TSS, its distance to target 
#' feature,strand and name of TSS/gene from annotationByGenicParts object
#' @param x an \code{\link[methylKit]{annotationByGenicParts}}  object
#' 
#' @return a \code{\link{data.frame}} containing row number of the target 
#' features,distance of target to nearest TSS, TSS/Gene name, TSS strand
#' @examples
#' data(methylKit)
#' gene.obj=readTranscriptFeatures(system.file("extdata", 
#' "refseq.hg18.bed.txt", package = "methylKit"))
#' ann=annotateWithGenicParts(methylDiff.obj, gene.obj)
#' df=getAssociationWithTSS(ann)
#' head(df)
#' 
#' @export
#' @docType methods
#' @rdname annotationByGenicParts-methods
setGeneric("getAssociationWithTSS", def=function(x) 
  standardGeneric("getAssociationWithTSS"))

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
#' This function plots a pie or bar chart for showing percentages of targets 
#' annotated by genic parts or other query features.
#' @param x an \code{\link[methylKit]{annotationByFeature}} or an 
#' \code{\link[methylKit]{annotationByGenicParts}} object
#' @param precedence TRUE|FALSE. If TRUE there will be a hierachy of annotation 
#' features when calculating numbers (with promoter>exon>intron precedence). 
#'  This option is only valid when x is a \code{annotationByGenicParts} object
#' @param col a vector of colors for piechart or the bar plot
#' @param ... graphical parameters to be passed to \code{pie} or \code{barplot} 
#' functions
#'
#' @examples
#' data(methylKit)
#' gene.obj=readTranscriptFeatures(system.file("extdata", 
#' "refseq.hg18.bed.txt", package = "methylKit"))
#' ann=annotateWithGenicParts(methylDiff.obj, gene.obj)
#' plotTargetAnnotation(ann,precedence=FALSE)
#' 
#'
#'
#' @return plots a piechart or a barplot for percentage of the target features
#' overlapping with annotation
#' 
#' @export
#' @docType methods
#' @rdname plotTargetAnnotation-methods
setGeneric("plotTargetAnnotation", def=function(x,precedence=TRUE,
            col=rainbow(length(x@annotation)),...) 
  standardGeneric("plotTargetAnnotation"))

#' @rdname plotTargetAnnotation-methods
#' @docType methods
#' @aliases plotTargetAnnotation,annotationByFeature-method
setMethod("plotTargetAnnotation", signature(x = "annotationByFeature"),
                    function(x,precedence,col,...){
    props=getTargetAnnotationStats(x,percentage=TRUE,
                                   precedence=precedence)

    if(precedence | ( is(x,"annotationByFeature") & 
                       !is(x,"annotationByGenicParts")) ){
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

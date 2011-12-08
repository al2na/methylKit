# PART THAT DEALS with annotation of CpGs and differential methylation events


# SECTION 1:
# reading annotation to Granges
# makes GRanges object from a given bed6 or bed12 file to granges object

# SECTION 1: S3 functions

# extracts exons from a bed12 file and puts them into Granges object
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
    grange.exons=GRanges(seqnames=as.character(rep.ref$V1),
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

# extracts exons from a bed12 file and puts them into Granges object
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
    grange.exons=GRanges(seqnames=as.character(rep.ref$V1),
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


# SECTION 1: S4 functions


# convert a data frame read-in from a bed file to a Granges object
setGeneric("convert.bed.df",function(bed) standardGeneric("convert.bed.df"))
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

# convert a data frame read-in from a bed file to a Granges object for exons
setGeneric("convert.bed2exons",function(bed.df) standardGeneric("convert.bed2exons"))
setMethod("convert.bed2exons" ,signature(bed.df = "data.frame" ),
                            function(bed.df){
                              
                            if(! check.bed.validity(bed.df,"exon")){stop("this is not a valid bed file")}
                            bed12.to.exons(bed.df)
})

# convert a data frame read-in from a bed file to a Granges object for introns
setGeneric("convert.bed2introns",function(bed.df) standardGeneric("convert.bed2introns"))
setMethod("convert.bed2introns" ,signature(bed.df = "data.frame" ),
                            function(bed.df){
                              
                            if(! check.bed.validity(bed.df,"exon")){stop("this is not a valid bed file")}                      
                            bed12.to.introns(bed.df)
})


# read a bed file and convert it to GRanges
# arguments:
# location: location of the file
# remove.unsual  : remove the chromomesomes with unsual names, mainly random chromsomes etc
# DETAILS: one bed track per file is only accepted, the bed files with multiple tracks will cause en error
setGeneric("read.bed", function(location,remove.unsual=T) standardGeneric("read.bed"))
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

# function reading exon intron, promoter structure from a given bed file
# ARGUMENTS      :
# location       : location of the bed file with 12 or more columns
# remove.unsual  : remove the chromomesomes with unsual names, mainly random chromsomes etc
# up.flank       : up-stream from TSS to detect promoter boundaries
# down.flank     : down-stream from TSS to detect promoter boundaries
# unique.prom    : get only the unique promoters, promoter boundaries will not have a gene name if you set this option to be TRUE
# DETAILS: one bed track per file is only accepted, the bed files with multiple tracks will cause en error
setGeneric("read.transcript.features", function(location,remove.unsual=T,up.flank=1000,down.flank=1000,unique.prom=T) standardGeneric("read.transcript.features"))
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

# a function to get upstream and downstream adjecent regions to a genomic feature such as CpG islands
# ARGUMENTS:
# GRanges object for the feature
# flank: number of basepairs for the flanking regions
# clean: If set to TRUE, flanks overlapping with other main features will be trimmed, and overlapping flanks will be removed
#        this will remove multiple counts when other features overlap with flanks
setGeneric("getFlanks", function(grange,flank=2000,clean=T) standardGeneric("getFlanks"))
setMethod("getFlanks", signature(grange= "GRanges"),
                    function(grange,flank=2000,clean=T){
          
                    shores=c( flank(grange,flank),flank(grange,flank,FALSE) )
                    if(clean){
                      shores=GenomicRanges::reduce(GenomicRanges::setdiff(shores, grange)) # ,erge overlapping shores remove CpG coordinates from all shores, das ist so cool!!
                    }
                    shores
})

# a function to read-in genomic features and their upstream and downstream adjecent regions such as CpG islands and their shores
# ARGUMENTS:
# location: for the bed file of the feature 
# flank   : number of basepairs for the flanking regions
# clean   : If set to TRUE, flanks overlapping with other main features will be trimmed
# remove.unsual : remove chromsomes with unsual names random, Un and antyhing with "_" character
#feature.flank.name: the names for feature and flank ranges, it should be a character vector of length 2. example: c("CpGi","shores")
# VALUE   :
# a GRangesList contatining one GRanges object for flanks and one for GRanges object for the main feature
setGeneric("read.feature.flank", function(location,remove.unsual=T,flank=2000,clean=T,feature.flank.name=NULL) standardGeneric("read.feature.flank") )
setMethod("read.feature.flank", signature(location = "character"),
                    function(location,remove.unsual,flank ,clean,feature.flank.name){
                    feat=read.bed(location,remove.unsual)
                    flanks=getFlanks(feat,flank=flank,clean=clean)
                    x=GRangesList(features=feat,flanks=flanks)
                    if(!is.null(feature.flank.name) & length(feature.flank.name)==2)
                    {
                      names(x)=feature.flank.name
                    }
                    x
})

# SECTION 2:
# annotate granges objects with annotations that read-in and converted to Granges objects


# A set of objects that will hold statistics about feature and annotation overlap
setClass("annotationByFeature", representation(members         ="matrix",
                                      annotation      ="numeric",
                                      hierarchical    ="numeric",
                                      num.annotation  ="numeric",
                                      num.hierarchical="numeric",
                                      no.of.OlapFeat  ="numeric",
                                      perc.of.OlapFeat="numeric"))

setClass("annotationByGenicParts", representation(dist.to.TSS   ="data.frame"),contains="annotationByFeature")


#new.obj=new("annotationByGenicParts",
#            members=matrix(c(1,2,3,4)),annotation=c(1,2,0,3,4),
#            hierarchical=c(a=1,b=2,c=0,d=3,e=4),
#            num.annotation  =c(1,2,0,3,4),
#            num.hierarchical=c(1,2,0,3,4),
#            no.of.OlapFeat  =c(1,2,0,3,4),
#            perc.of.OlapFeat=c(1,2,0,3,4),
#            dist.to.TSS     =c(1,2,0,3,4) )
                                    

setMethod("show", "annotationByGenicParts", function(object) {
  
  cat("summary of target set annotation with genic parts\n");cat(nrow(object@members));cat(" rows in target set\n--------------\n")
  cat("--------------\n")
  cat("percentage of target features overlapping with annotation :\n");print(object@annotation);cat("\n\n")
  cat("percentage of target features overlapping with annotation (with promoter>exon>intron precedence) :\n"); print(object@hierarchical) ;cat("\n\n")
  cat("percentage of annotation boundaries with feature overlap :\n");print(object@perc.of.OlapFeat);cat("\n\n")  
  cat("summary of distances to the nearest TSS :\n")
  print(summary(abs(object@dist.to.TSS[,2])))
})


# SECTION 2: S3 FUNCTIONS 
# these shouldn't be exported

annotate.gr.WithGenicParts<-function(gr,prom,exon,intron,strand=F)
{
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

  hierarchical=c(promoter  =100*sum(memb$prom>0)/nrow(memb) ,
               exon      =100*sum(memb$exon>0 & memb$prom==0)/nrow(memb) ,
               intron    =100*sum(memb$intron>0 & memb$exon==0 & memb$prom==0)/nrow(memb) ,
               intergenic=100*sum(rowSums(memb)==0)/nrow(memb) )

  num.hierarchical=c( promoter  =sum(memb$prom>0) ,
                      exon      =sum(memb$exon>0 & memb$prom==0),
                      intron    =sum(memb$intron>0 & memb$exon==0 & memb$prom==0) ,
                      intergenic=sum(rowSums(memb)==0) )  
  
  numberOfOlapFeat=c(promoter=sum(countOverlaps(prom,gr)>0),
    exon=sum(countOverlaps(exon,gr)>0),
    intron=sum(countOverlaps(intron,gr)>0) )
  percOfOlapFeat =100*numberOfOlapFeat/c(length(prom),length(exon),length(intron))
  return(list(members=memb,annotation=annotation,hierarchical=hierarchical,num.annotation=num.annotation,
              num.hierarchical=num.hierarchical,
              numberOfOlapFeat=numberOfOlapFeat,percOfOlapFeat=percOfOlapFeat) )
  
}

distance2nearestFeature<-function(g.idh,tss)
{
  
  elementMetadata(g.idh)=DataFrame(elementMetadata(g.idh),orig.row=1:length(g.idh))
  id.col    =ncol(elementMetadata(g.idh))+5 # get the row number column
  tss.id.col=ncol(elementMetadata(g.idh))+5+6 # get the id column for tss
  met.tss   = .nearest.2bed(g.idh, tss)

  dist.col=ncol(met.tss)
  met.tss[met.tss$end<met.tss$start.y & met.tss$strand.y=="+",dist.col] = -1* met.tss[met.tss$end<met.tss$start.y & met.tss$strand.y=="+",dist.col]

  met.tss[met.tss$end>met.tss$start.y & met.tss$strand.y=="-",dist.col] = -1* met.tss[met.tss$end>met.tss$start.y & met.tss$strand.y=="-",dist.col]

  res=met.tss[order(met.tss[,id.col]),c(id.col,dist.col,tss.id.col,tss.id.col-1)]
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






# function to annotate given GRanges object with promoter,exon,intron & intergenic values
# ARGUMENTS:
# GRanges.obj      : a granges object storing chromosome locations to be annotated
# GRangesList.obj  : A GRangesList object containing GRanges object for promoter,exons,introns and TSSes, or simply output of read.transcript.features function
# strand           : If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
# VALUE   :
# a annotationByGenicParts object
setGeneric("annotate.WithGenicParts", function(target,GRangesList.obj,strand=F) standardGeneric("annotate.WithGenicParts") )
setMethod("annotate.WithGenicParts", signature(target= "GRanges",GRangesList.obj="GRangesList"),
                    function(target,GRangesList.obj,strand){
                      
                      a.list    =annotate.gr.WithGenicParts(target,GRangesList.obj$promoters,GRangesList.obj$exons,GRangesList.obj$introns,strand=strand)
                      dist2TSS  =distance2nearestFeature(target,GRangesList.obj$TSSes)

                      new("annotationByGenicParts",
                                  members         =as.matrix(a.list$members),
                                  annotation      =a.list$annotation,
                                  hierarchical    =a.list$hierarchical,
                                  num.annotation  =a.list$num.annotation,
                                  num.hierarchical=a.list$num.hierarchical,
                                  no.of.OlapFeat  =a.list$numberOfOlapFeat,
                                  perc.of.OlapFeat=a.list$percOfOlapFeat,
                                  dist.to.TSS     = dist2TSS )
})

setMethod("annotate.WithGenicParts", signature(target = "methylDiff",GRangesList.obj="GRangesList"),
                    function(target,GRangesList.obj,strand){
                      gr=as(target,"GRanges")
                      annotate.WithGenicParts(gr,GRangesList.obj,strand)
})


# function to annotate given GRanges object with promoter,exon,intron & intergenic values
# ARGUMENTS:
# gr      : a granges object storing chromosome locations to be annotated
# feature : a granges object storing chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
# flank   : a granges object storing chromosome locations of the flanks of the feature
# strand  : If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
# VALUE   :
# a annotationByFeature object
setGeneric("annotate.WithFeature.Flank", function(target,feature,flank,feature.name="feat",flank.name="flank",strand=F) standardGeneric("annotate.WithFeature.Flank") )
setMethod( "annotate.WithFeature.Flank", signature(target = "GRanges",feature="GRanges",flank="GRanges"),
                    function(target, feature, flank,feature.name,flank.name,strand){
                      
                      if( ! strand){strand(target)="*"}
                      memb=data.frame(matrix(rep(0,length(gr)*2),ncol=2) )  ;colnames(memb)=c(name1,name2)
                      memb[countOverlaps(target,feature)>0,1]=1
                      memb[countOverlaps(target,flank)>0,2]=1
                    
                      annotation=c(100*sum(memb[,1]>0)/nrow(memb) ,
                                   100*sum(memb[,2]>0)/nrow(memb) ,
                                   100*sum(rowSums(memb)==0)/nrow(memb) )
                      names(annotation)=c(name1,name2,"other")

                      num.annotation=c( sum(memb[,1]>0) , sum(memb[,2]>0) , sum(rowSums(memb)==0) )
                      names(num.annotation)=c(name1,name2,"other")                      
                      
                      
                      hierarchical=c(100*sum(memb[,1]>0)/nrow(memb) ,
                                     100*sum(memb[,2]>0 & memb[,1]==0)/nrow(memb) ,
                                     100*sum(rowSums(memb)==0)/nrow(memb) )
                      names(hierarchical)=c(name1,name2,"other")
                      num.hierarchical=c( sum(memb[,1]>0)  , sum(memb[,2]>0 & memb[,1]==0) , sum(rowSums(memb)==0)  )
                      names(num.hierarchical)=c(name1,name2,"other")                      
                      
                      numberOfOlapFeat=c(sum(countOverlaps(feature,target)>0),
                                         sum(countOverlaps(flank,target)>0) )
                      names(numberOfOlapFeat)=c(name1, name2)
                      percOfOlapFeat =100*numberOfOlapFeat/c(length(feature),length(flank) )
                    
                      #return(list(members=memb,annotation=annotation,hierarchical=hierarchical,
                      #            numberOfOlapFeat=numberOfOlapFeat,percOfOlapFeat=percOfOlapFeat) )                      
                      new("annotationByFeature",
                          members         =as.matrix(memb),
                          annotation      =annotation,
                          hierarchical    =hierarchical,
                          num.annotation  =num.annotation,
                          num.hierarchical=num.hierarchical,
                          no.of.OlapFeat  =numberOfOlapFeat,
                          perc.of.OlapFeat=percOfOlapFeat)
                      
})


setMethod("annotate.WithFeature.Flank", signature(target= "methylDiff",feature="GRanges",flank="GRanges"),
                    function(target, feature, flank,feature.name,flank.name,strand){
                      gr=as(target,"GRanges")
                      annotate.WithFeature.Flank(gr,feature, flank,feature.name,flank.name,strand)
})


# function to annotate given GRanges object with promoter,exon,intron & intergenic values
# ARGUMENTS:
# gr      : a granges object storing chromosome locations to be annotated
# feature : a granges object storing chromosome locations of a feature (can be CpG islands, ChIP-seq peaks, etc)
# strand  : If set to TRUE, annotation features and target features will be overlapped based on strand  (def:FAULT)
# VALUE   :
# a annotationByFeature object
setGeneric("annotate.WithFeature", function(target,feature,strand=F,extend=0,feature.name="feat1") standardGeneric("annotate.WithFeature") )
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
                          hierarchical    =0,
                          num.annotation  =num.annotation,
                          num.hierarchical=0,
                          no.of.OlapFeat  =numberOfOlapFeat,
                          perc.of.OlapFeat=percOfOlapFeat)



})


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
#' For example, if a target feature overlaps with an exon
#' @param x a \code{annotationByFeature}  object
#' 
#' @return RETURNS a matrix showing overlap of target features with annotation features. 1 for overlap, 0 for non-overlap
#' 
#' @aliases getMembers,-methods getMembers,annotationByFeature-method
#' @export
#' docType methods
#' rdname annotationByFeature-methods                      
setGeneric("getMembers", def=function(x) standardGeneric("getMembers"))

#' @rdname annotationByFeature-methods
#' @aliases getMembers,annotationByFeature,ANY-method
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
#' @param hierarchical TRUE|FALSE. If TRUE there will be a hierachy of annotation features when calculating numbers (with promoter>exon>intron precedence)
#' That means if a feature overlaps with a promoter it will be counted as promoter overlapping only, or if it is overlapping with a an exon but not a promoter, 
#' it will be counted as exon overlapping only whether or not it overlaps with an intron.
#'
#' @return RETURNS  a vector of percentages or counts showing quantity of target features overlapping with annotation
#' 
#' @aliases getTargetAnnotation,-methods getTargetAnnotation,annotationByFeature-method
#' @export
#' docType methods
#' rdname annotationByFeature-methods
setGeneric("getTargetAnnotation", def=function(x,percentage=T,hierarchical=T) standardGeneric("getTargetAnnotation"))

#' docType methods
#' @rdname annotationByFeature-methods
#' @aliases getTargetAnnotation,annotationByFeature,ANY-method
setMethod("getTargetAnnotation", signature(x = "annotationByFeature"),
                    function(x,percentage ,hierarchical ){                      
                      if(percentage){
                        if(hierarchical){return(x@hierarchical)
                        }else{return(x@annotation)}
                      }else{
                        if(hierarchical){return(x@num.hierarchical)
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
#' 
#' @aliases getFeaturesWithTargets,-methods getFeaturesWithTargets,annotationByFeature-method
#' @export
#' docType methods
#' rdname annotationByFeature-methods
setGeneric("getFeaturesWithTargets", def=function(x,percentage=T) standardGeneric("getFeaturesWithTargets"))

#' docType methods
#' @rdname annotationByFeature-methods
#' @aliases getFeaturesWithTargets,annotationByFeature,ANY-method
setMethod("getFeaturesWithTargets", signature(x = "annotationByFeature" ),
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
#' 
#' @aliases getAssociationWithTSS,-methods getAssociationWithTSS,annotationByGenicParts-method
#' @export
#' docType methods
#' rdname annotationByGenicParts-methods
setGeneric("getAssociationWithTSS", def=function(x) standardGeneric("getAssociationWithTSS"))

#' @rdname annotationByGenicParts-methods
#' docType methods
#' @aliases getAssociationWithTSS,annotationByGenicParts,ANY-method
setMethod("getAssociationWithTSS", signature(x = "annotationByGenicParts"),
                    function(x){
                      return(x@dist.to.TSS)
                      
})

# PLOTTING FUNCTIONS

#' Get distance to nearest TSS and gene id from annotationByGenicParts
#'
#' This accessor function gets the nearest TSS, its distance to target feature,strand and name of TSS/gene from annotationByGenicParts object
#' @param x a \code{annotationByFeature}  object
#' @param hierarchical TRUE|FALSE. If TRUE there will be a hierachy of annotation features when calculating numbers (with promoter>exon>intron precedence)
#' @param col a vector of colors for piechart or the par plot
#' @param ... graphical parameters to be passed to \code{pie} or \code{barplot} functions
#'
#' @usage \code{plotTargetAnnotation(x,hierarchical=T,col,...)}
#'
#' @return plots a piechart or a barplot for percentage of the target features overlapping with annotation
#' 
#' @aliases plotTargetAnnotation,-methods plotTargetAnnotation,annotationByFeature-method
#' @export
#' docType methods
#' rdname plotTargetAnnotation-methods
setGeneric("plotTargetAnnotation", def=function(x,hierarchical=T,col=rainbow(length(x@annotation)),...) standardGeneric("plotTargetAnnotation"))

#' @rdname plotTargetAnnotation-methods
#' docType methods
#' @aliases plotTargetAnnotation,annotationByFeature,ANY-method
setMethod("plotTargetAnnotation", signature(x = "annotationByFeature"),
                    function(x,hierarchical,col,...){
                      props=getTargetAnnotation(x,hierarchical)

                      if(hierarchical){
                        slice.names=names(props)
                        #names(props)=paste(names(props),paste(round(props),"%"),sep=" ")
                        names(props)=paste( paste(round(props),"%"),sep=" ")

                        pie(props,cex=0.9,col=col,...)
                         legend("topright",legend=slice.names,fill=rainbow(length(props)) )

                      }
                      else{
                        
                        mids=barplot(props,col=col,...)  
                        text(mids,y=props,labels=paste(paste(round(props),"%",sep="")),pos=1) 
                      }

})



# SECTION 3:
# annotate ML objects with annotations read-in and converted to Granges objects

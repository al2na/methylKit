
#' @export
read <- function(location,sample.id,assembly,dbtype=NA,
                          pipeline="amp",header=TRUE,skip=0,sep="\t",
                          context="CpG",resolution="base",
                          treatment,dbdir=getwd(),mincov=10)
{
  .Deprecated("methRead")
  ## use new function, or remainder of myOldFunc
  methRead(location,sample.id,assembly,dbtype ,
           pipeline ,header ,skip ,sep ,
           context ,resolution ,
           treatment,dbdir ,mincov)
}



#' @export
read.bismark <- function(location,sample.id,assembly,
                                  save.folder=NULL,
                                  save.context=c("CpG"),read.context="CpG",
                                  nolap=FALSE,mincov=10,
                                  minqual=20,phred64=FALSE
                                  ,treatment,save.db=FALSE)
{
  .Deprecated("processBismarkAln")
  ## use new function, or remainder of myOldFunc
  processBismarkAln(location,sample.id,assembly,
                    save.folder ,
                    save.context ,read.context ,
                    nolap ,mincov,
                    minqual ,phred64 
                    ,treatment,save.db)
}




#' @export
adjust.methylC<-function(...){
  .Deprecated("adjustMethylC")
  adjustMethylC(...)
}


#export(get.methylDiff)
#' @export
get.methylDiff<-function(...){
  .Deprecated("getMethylDiff")
  getMethylDiff(...)
}

#' @export
annotate.WithFeature<-function(...){
  .Deprecated("annotateWithFeature")
  annotateWithFeature(...)
}

#' @export
annotate.WithFeature.Flank<-function(...){
  .Deprecated("annotateWithFeatureFlank")
  annotateWithFeatureFlank(...)
}


#' @export
annotate.WithGenicParts<-function(...){
  .Deprecated("annotateWithGenicParts")
  annotateWithGenicParts(...)
}


#' @export
read.bed<-function(...){
  .Deprecated("readBed")
  readBed(...)
}

#' @export
read.feature.flank<-function(...){
  .Deprecated("readFeatureFlank")
  readFeatureFlank(...)
}


#' @export
read.transcript.features<-function(...){
  .Deprecated("readTranscriptFeatures")
  readTranscriptFeatures(...)
}



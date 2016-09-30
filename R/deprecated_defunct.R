
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


#' Deprecated/Defunct functions
#' 
#' These are deprecated or defunct functions. Most of them 
#' are replaced by genomation functions.
#' 
#' @docType methods
#' @aliases annotate.WithFeature
#' @rdname genomation-deprecated
#' @export
annotate.WithFeature<-function(){
  .Deprecated("genomation::annotateWithFeature")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")}

#' @rdname  genomation-deprecated
#' @export
annotate.WithFeature.Flank<-function(){
  .Deprecated("genomation::annotateWithFeatureFlank")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")
}


#' @rdname genomation-deprecated
#' @export
annotate.WithGenicParts<-function(){
  .Deprecated("genomation::annotateWithGeneParts")
  message("Use functions in genomation package from Bioconductor")
}


#' @rdname genomation-deprecated
#' @export
read.bed<-function(){
  .Deprecated("genomation::readBed")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")}

#' @rdname genomation-deprecated
#' @export
read.feature.flank<-function(){
  .Deprecated("genomation::readFeatureFlank")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")}

#' @rdname genomation-deprecated
#' @export
read.transcript.features<-function(){
  .Deprecated("genomation::readTranscriptFeatures")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")}



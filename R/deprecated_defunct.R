#' Deprecated/Defunct functions
#' 
#' These are deprecated or defunct functions. Most of them 
#' are replaced by genomation functions. See the vignette for 
#' examples on how to use genomation functions for annotation
#' purposes.
#' 
#' @docType methods
#' @rdname methylKit-defunct
#' @name methylKit-defunct
#' @aliases annotate.WithFeature
#' @export
annotate.WithFeature<-function(){
  .Defunct("genomation::annotateWithFeature")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")}

#' @rdname  methylKit-defunct
#' @aliases annotate.WithFeature.Flank
#' @export
annotate.WithFeature.Flank<-function(){
  .Defunct("genomation::annotateWithFeatureFlank")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")
}


#' @rdname methylKit-defunct
#' @aliases annotate.WithGenicParts
#' @export
annotate.WithGenicParts<-function(){
  .Defunct("genomation::annotateWithGeneParts")
  message("Use functions in genomation package from Bioconductor")
}


#' @rdname methylKit-defunct
#' @aliases read.bed
#' @export
read.bed<-function(){
  .Defunct("genomation::readBed")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")}

#' @rdname methylKit-defunct
#' @aliases read.feature.flank
#' @export
read.feature.flank<-function(){
  .Defunct("genomation::readFeatureFlank")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")}

#' @rdname methylKit-defunct
#' @aliases read.transcript.features
#' @export
read.transcript.features<-function(){
  .Defunct("genomation::readTranscriptFeatures")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")}

#' @rdname methylKit-defunct
#' @aliases getFeatsWithTargetsStats
#' @export
getFeatsWithTargetsStats<-function(){
  .Defunct("genomation::getFeatsWithTargetsStats")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")}

#' @rdname methylKit-defunct
#' @aliases getFlanks 
#' @export
getFlanks<-function(){
  .Defunct("genomation::getFlanks")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")} 

#' @rdname methylKit-defunct
#' @aliases getMembers 
#' @export
getMembers<-function(){
  .Defunct("genomation::getMembers")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")} 

#' @rdname methylKit-defunct
#' @aliases getTargetAnnotationStats 
#' @export
getTargetAnnotationStats<-function(){
  .Defunct("genomation::getTargetAnnotationStats")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")} 

#' @rdname methylKit-defunct
#' @aliases plotTargetAnnotation 
#' @export
plotTargetAnnotation<-function(){
  .Defunct("genomation::plotTargetAnnotation")
  message("Use functions in genomation package from Bioconductor\n",
          "See vignette for examples.")} 


#' @rdname methylKit-defunct
#' @name methylKit-defunct
#' @aliases read
#' @export
read <- function()
{
  .Defunct("methRead")
}



#' @rdname methylKit-defunct
#' @name methylKit-defunct
#' @aliases read.bismark
#' @export
read.bismark <- function()
{
  .Defunct("processBismarkAln")
  
}



#' @rdname methylKit-defunct
#' @name methylKit-defunct
#' @aliases adjust.methylC
#' @export
adjust.methylC<-function(){
  .Defunct("adjustMethylC")
}

#' @rdname methylKit-defunct
#' @name methylKit-defunct
#' @aliases get.methylDiff
#' @export
get.methylDiff<-function(){
  .Defunct("getMethylDiff")
}

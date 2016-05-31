
#' @export
read <- function(location,sample.id,assembly,dbtype=NA,
                          pipeline="amp",header=T,skip=0,sep="\t",
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



#export(adjust.methylC)
#export(get.methylDiff)

#export(annotate.WithFeature)
#export(annotate.WithFeature.Flank)
#export(annotate.WithGenicParts)

#export(read.bed)
#export(read.bismark)
#export(read.feature.flank)
#export(read.transcript.features)


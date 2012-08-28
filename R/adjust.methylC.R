
#' Adjust measured 5mC levels using 5hmC levels
#' 
#' Measured 5mC levels via bisulfite sequencing might be a combination of 5hmC and 5mC levels since bisulfite sequencing can not distinguish
#' between the two. This function can adjust 5mC levels of a bisulfite sequencing experiment
#' if the user supplies corresponding 5hmC levels from the same sample.
#'
#' @param mc a \code{methylRawList} or \code{methylRaw} containing 5mC levels of a sample or set of samples
#' @param hmc a \code{methylRawList} or \code{methylRaw} containing 5hmC levels of a sample or set of samples. 
#'            If a \code{methylRawList} given the sample order should be same as "mc" \code{methylRawList} object.
#'
#' @return returns adjusted 5-methyl cytosine levels in the form of \code{methylRawList} or \code{methylRaw} object depending on the input object
#' @usage adjust.methylC(mc,hmc)
#' @examples
#' 
#' # read 5hmC and 5mC files
#' hmc.file=system.file("extdata", "test1.myCpG.txt", package = "methylKit")
#' mc.file =system.file("extdata", "test2.myCpG.txt", package = "methylKit")
#' 
#' my5hmC=read( hmc.file,sample.id="hmc",assembly="hg18")
#' my5mC =read( mc.file,sample.id="mc",assembly="hg18")
#' 
#' # adjusting the 5mC levels using 5hmC levels
#' adjusted.5mC=adjust.methylC(my5mC,my5hmC)
#' 
#' @references
#' 1. Booth, Branco, et al. (2012). Quantitative Sequencing of 5-Methylcytosine and 5-Hydroxymethylcytosine at Single-Base Resolution. Science, 934
#' 
#' 2. Yu, Hon, et al. (2012). Base-resolution analysis of 5-hydroxymethylcytosine in the Mammalian genome. Cell, 149(6), 1368-80.
#' @export
#' @docType methods
#' @rdname adjust.methylC
setGeneric("adjust.methylC", function(mc,hmc) standardGeneric("adjust.methylC") )

#' @rdname adjust.methylC
#' @aliases adjust.methylC,methylRaw,methylRaw-method
setMethod("adjust.methylC", c("methylRaw","methylRaw"),function(mc,hmc){
  
  lst=new("methylRawList",list(mc,hmc),treatment=c(1,0))
  data=getData(unite(lst))
  
  diff=(data$numCs1)-round(data$coverage1*(data$numCs2/data$coverage2))
  diff[diff<0]=0
  data$numCs1=diff
  data$numTs1=data$coverage1-data$numCs1
  colnames(data)[6:8]=c("coverage","numCs","numTs")
  new("methylRaw",data[,1:8],sample.id=mc@sample.id,  assembly=mc@assembly, 
                             context =mc@context,     resolution=mc@resolution)
  
})


#' @rdname adjust.methylC
#' @aliases adjust.methylC,methylRawList,methylRawList-method
setMethod("adjust.methylC", c("methylRawList","methylRawList"),function(mc,hmc){
  
  # check lengths equal if not give error
  if(length(mc) != length(hmc)){stop("lengths of methylRawList objects should be same\n")}
  
  my.list=list()
  for(i in 1:length(mc)){
    my.list[[i]]=adjust.methylC(mc[[i]],hmc[[i]])
  }
  new("methylRawList",my.list,treatment=mc@treatment )
  
})


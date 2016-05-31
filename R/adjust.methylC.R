
#' Adjust measured 5mC levels using 5hmC levels
#' 
#' Measured 5mC levels via bisulfite sequencing might be a combination of 5hmC 
#' and 5mC levels since bisulfite sequencing can not distinguish
#' between the two. This function can adjust 5mC levels of a bisulfite 
#' sequencing experiment
#' if the user supplies corresponding 5hmC levels from the same sample.
#'
#' @param mc a \code{methylRawList}, \code{methylRaw}, \code{methylRawDB} 
#' or \code{methylRawListDB} containing 5mC levels of a sample or set of samples
#' @param hmc a \code{methylRawList}, \code{methylRaw}, \code{methylRawDB} or
#'  \code{methylRawListDB} containing 5hmC levels of a sample or set of samples. 
#'            If a \code{methylRawList} or \code{methylRawListDB} given the 
#'            sample order should be same as "mc" \code{methylRawList} or
#'             \code{methylRawListDB} object.
#' @param chunk.size Number of rows to be taken as a chunk for processing the 
#' \code{methylRawDB} or \code{methylRawListDB} objects (default: 1e6)
#' @param save.db A Logical to decide whether the resulting object should be 
#' saved as flat file database or not, default: explained in Details sections  
#' @param ... optional Arguments used when save.db is TRUE
#'            
#'            \code{suffix}
#'                  A character string to append to the name of the output 
#'                  flat file database, 
#'                  only used if save.db is true, default actions: 
#'                  append \dQuote{_filtered} to current filename 
#'                  if database already exists or generate new file with 
#'                  filename \dQuote{sampleID_filtered}
#'                  
#'            \code{dbdir} 
#'                  The directory where flat file database(s) should 
#'                  be stored, defaults
#'                  to getwd(), working directory for newly stored databases
#'                  and to same directory for already existing database
#'                  
#'            \code{dbtype}
#'                   The type of the flat file database, currently only 
#'                   option is "tabix"
#'                   (only used for newly stored databases)
#'
#'
#' @return returns adjusted 5-methyl cytosine levels in the form of 
#' \code{methylRawList}, \code{methylRaw}, \code{methylRawDB} or 
#' \code{methylRawListDB} object depending on the input object
#' @usage adjustMethylC(mc,hmc,save.db,...,chunk.size)
#' @examples
#' 
#' # read 5hmC and 5mC files
#' hmc.file=system.file("extdata", "test1.myCpG.txt", package = "methylKit")
#' mc.file =system.file("extdata", "test2.myCpG.txt", package = "methylKit")
#' 
#' my5hmC=methRead( hmc.file,sample.id="hmc",assembly="hg18")
#' my5mC =methRead( mc.file,sample.id="mc",assembly="hg18")
#' 
#' # adjusting the 5mC levels using 5hmC levels
#' adjusted.5mC=adjustMethylC(my5mC,my5hmC)
#' 
#' 
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with 
#' \code{methylRawDB} or \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing 
#' large-sized objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work 
#' for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the
#'  \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects 
#' as \code{methylRawDB} and \code{methylRawListDB}, 
#' while being per default FALSE for \code{methylRaw} and \code{methylRawList}. 
#' If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the database
#'  allows the calculation in-memory, 
#' then you might change the value of this parameter.
#' 
#' @references
#' 1. Booth, Branco, et al. (2012). Quantitative Sequencing of 5-Methylcytosine 
#' and 5-Hydroxymethylcytosine at Single-Base Resolution. Science, 934
#' 
#' 2. Yu, Hon, et al. (2012). Base-resolution analysis of 
#' 5-hydroxymethylcytosine in the Mammalian genome. Cell, 149(6), 1368-80.
#' @export
#' @docType methods
#' @rdname adjustMethylC
setGeneric("adjustMethylC", function(mc,hmc,
                                      save.db=FALSE,...,chunk.size=1e6) 
  standardGeneric("adjustMethylC") )

#' @rdname adjustMethylC
#' @aliases adjustMethylC,methylRaw,methylRaw-method
setMethod("adjustMethylC", c("methylRaw","methylRaw"),
          function(mc,hmc,save.db=FALSE,...){
  
  lst=new("methylRawList",list(mc,hmc),treatment=c(1,0))
  data=getData(unite(lst))
  
  diff=(data$numCs1)-round(data$coverage1*(data$numCs2/data$coverage2))
  diff[diff<0]=0
  data$numCs1=diff
  data$numTs1=data$coverage1-data$numCs1
  colnames(data)[5:7]=c("coverage","numCs","numTs")
  
  if(!save.db) {
    new("methylRaw",data[,1:7],sample.id=mc@sample.id,  assembly=mc@assembly, 
                               context =mc@context,     resolution=mc@resolution)
  } else {
    
    #print(names(as.list(match.call())))
    # catch additional args 
    args <- list(...)
    #print(args)
    
    if( !( "dbdir" %in% names(args)) ){
      dbdir <- .check.dbdir(getwd())
    } else { dbdir <- .check.dbdir(args$dbdir) }
#                         if(!( "dbtype" %in% names(args) ) ){
#                           dbtype <- "tabix"
#                         } else { dbtype <- args$dbtype }
    if(!( "suffix" %in% names(args) ) ){
      suffix <- paste0("_","adjusted")
    } else { 
      suffix <- args$suffix
      suffix <- paste0("_",suffix)
    }
    
    # create methylRawDB
    obj <- makeMethylRawDB(df=data[,1:7],dbpath=dbdir,dbtype="tabix",
                           sample.id=paste0(mc@sample.id,suffix),
                           assembly=mc@assembly,context=mc@context,
                           resolution=mc@resolution)
    obj@sample.id <- mc@sample.id
    
    obj
  }  
  
})


#' @rdname adjustMethylC
#' @aliases adjustMethylC,methylRawList,methylRawList-method adjust.methylC
setMethod("adjustMethylC", c("methylRawList","methylRawList"),
          function(mc,hmc,save.db=FALSE,...){
  
  # check lengths equal if not give error
  if(length(mc) != length(hmc)){
    stop("lengths of methylRawList objects should be same\n")}
  
  my.list=list()
  for(i in 1:length(mc)){
    my.list[[i]]=adjustMethylC(mc[[i]],hmc[[i]])
  }
  
  if(!save.db) {
  new("methylRawList",my.list,treatment=mc@treatment )
  
  } else {
    args <- list(...)
    
    if( !( "dbdir" %in% names(args)) ){
      dbdir <- .check.dbdir(getwd())
    } else { dbdir <- .check.dbdir(args$dbdir) }
    
    my.list=list()
    for(i in 1:length(mc)){
      my.list[[i]]=adjustMethylC(mc[[i]],hmc[[i]],
                                  save.db=TRUE,dbdir=basename(dbdir),...)
    }
    
    new("methylRawListDB", my.list,treatment=mc@treatment)
}
  
})




#' Get methylation percentage from sorted Bismark alignments
#'
#' The function calls methylation percentage per base from sorted Bismark 
#' SAM or BAM 
#' files and reads methylation information as methylKit objects.
#'  Bismark is a popular aligner for 
#' high-throughput bisulfite sequencing experiments and it outputs its results in 
#' SAM format by default, and can be converted to BAM. 
#' Bismark SAM/BAM format contains
#' aligner specific tags which are absolutely necessary for methylation 
#' percentage calling using \code{processBismarkAln}. 
#' SAM/BAM files from other aligners will not work with this function.
#'
#' @param location location of sam or bam file(s). If multiple files are 
#'                  given this 
#'                  argument must be a list.
#' @param sample.id the id(s) of samples in the same order as file.  
#'                  If multiple sam files are given this arugment must be a list.
#' @param save.folder The folder which will be used to save methylation call files,
#'                     if set to NULL no methylation call file will be saved 
#'                     as a text file.
#'                     The files saved can be read into R in less time using 
#'                     \code{methRead} 
#'                     function in \code{methylKit} 
#' @param save.context A character vector consisting following strings: 
#' "CpG","CHG","CHH". 
#'                    The methylation percentages for these methylation contexts
#'                     will be saved to save.folder
#' @param read.context One of the 'CpG','CHG','CHH' or 'none' strings. 
#'                     Determines what type of methylation context will be 
#'                     read-in 
#'                     to the memory which can be immediately used for analysis.
#'                     If given as 'none', processBismarkAln will not return
#'                      any object,
#'                     but if a save.folder argument given it will save the 
#'                     methylation percentage call files.
#' @param assembly string that determines the genome assembly. Ex: mm9,hg18 etc.
#' This is just a string for book keeping. It can be any string. Although,
#' when using multiple files from the same assembly, this string should be 
#' consistent in each object.
#' @param nolap   if set to TRUE and the SAM/BAM file has paired-end reads,
#'                 the one 
#'                read of the overlapping paired-end read pair will be ignored 
#'                for methylation calling.
#' @param mincov  minimum read coverage to call a methylation status for a base.
#' @param minqual minimum phred quality score to call a methylation status for a base.
#' @param phred64 logical (default: FALSE) you will not need to set this TRUE, 
#'                Currently bismark gives only phred33 scale
#' @param treatment treatment vector only to be used when location and sample.id 
#'                  parameters are \code{list}s and you are trying to read-in 
#'                  multiple samples that are related to eachother in down-stream 
#'                  analysis. 
#' @param save.db A Logical to decide whether the resulting object should be saved 
#'                as flat file database or not ( default: FALSE). With 
#'                the default value, a text file containing methylation values
#'                will be saved.
#'                If TRUE, database will either be saved to location 
#'                \code{save.folder} or 
#'                if this is NULL, to a new folder in the current 
#'                working directory 
#'                named after this scheme: 
#'                "methylDB <Date> <3randomlettersornumbers>"
#'
#' @return \code{methylRaw} or \code{methylRawList} object
#'
#' @note
#' SAM files should be sorted with samtools sort or unix sort. Other sorting
#' methods can alter the order of fields(columns) in the SAM file and that will
#' result in an error when using \code{processBismarkAln()}.
#' 
#' @export
#' @docType methods
#' @rdname processBismarkAln-methods
#' @aliases processBismarkAln read.bismark
#'
#' @examples
#' 
#' # reading one bismark file:
#' my.file=system.file("extdata", "test.fastq_bismark.sorted.min.sam", 
#'                               package = "methylKit")
#' obj=processBismarkAln(my.file,"test",assembly="hg18",save.folder=NULL,
#'                  save.context="CpG",read.context="CpG")
#'  
#' # reading multiple files
#' file.list2=list(system.file("extdata", "test.fastq_bismark.sorted.min.sam", 
#' package = "methylKit"),
#'                system.file("extdata", "test.fastq_bismark.sorted.min.sam", 
#'                package = "methylKit"),
#'               system.file("extdata", "test.fastq_bismark.sorted.min.sam", 
#'               package = "methylKit"),
#'                system.file("extdata", "test.fastq_bismark.sorted.min.sam",
#'                 package = "methylKit"))
#' 
#'  objs=processBismarkAln(location=file.list2
#'              ,sample.id=list("test1","test2","ctrl1","ctrl1"),assembly="hg18",
#'              save.folder=NULL,save.context=NULL,read.context="CpG",
#'              nolap=FALSE,mincov=10,minqual=20,phred64=FALSE,
#'              treatment=c(1,1,0,0))
#' 
setGeneric("processBismarkAln", function(location,sample.id,assembly,
                                         save.folder=NULL,
                                    save.context=c("CpG"),read.context="CpG",
                                    nolap=FALSE,mincov=10,
                                    minqual=20,phred64=FALSE
                                    ,treatment,save.db=FALSE) 
  standardGeneric("processBismarkAln"))

#' @aliases processBismarkAln,character,character,character-method
#' @rdname processBismarkAln-methods
setMethod("processBismarkAln", signature(location = "character",
                                         sample.id= "character",
                                    assembly= "character"),
function(location,sample.id,assembly,save.folder,save.context
         ,read.context,
         nolap,mincov,minqual,phred64,save.db){
  
  # check if file exists
  if(! file.exists(location) ){
    stop("File supplied as the 'location' argument doesn't exist\n",
    "can not find file at: ",location,"\n")}
  
  # check output types
  if(! all( save.context %in% c("CpG","CHG","CHH")) ){
    stop("wrong 'save.context' argument supplied, only 'CpG', 'CHG' or ", 
         "'CHH' accepted")
  }
  
  # read.type
   if( !( read.context %in% c("CpG","CHG","CHH","none") ) ){
    stop("wrong 'read.context' argument supplied, only 'CpG', 'CHG', 'CHH' or ",
          "'none' accepted")
  }
      
  # if output.folder NULL create temp
  # check output folder if not create
  out.files=list("CpG"="","CHG"="","CHH"="") # list of output files            
  temp.files=FALSE
  if(is.null(save.folder)){

    for(mytype in read.context){
      out.files[[mytype]]=tempfile(pattern = paste("methylKit_temp",mytype,sep=".") )
    }
    temp.files=TRUE # set there are temp files to be deleted
    
  }else{
    try(
      dir.create(save.folder, showWarnings = TRUE, 
                 recursive = TRUE, mode = "0777")
      )

                          
    for(mytype in unique(c(save.context,read.context)) ){
      out.files[[mytype]]=paste(save.folder,"/",sample.id,"_",mytype,".txt",sep="") 
    }
    
  }
  
  # call the Rcpp function 
  methCall(read1 = location, type = "bam", nolap = nolap, minqual = minqual, 
           mincov = mincov, phred64 = phred64, CpGfile = out.files[["CpG"]], 
           CHHfile = out.files[["CHH"]], CHGfile = out.files[["CHG"]] ) 
  

  # read the result
  if(read.context != "none"){
    cat("Reading methylation percentage per base for sample:",sample.id,"\n\n")
    if(save.db) { dbtype="tabix"; 
      if(is.null(save.folder)) dbdir=getwd() else  dbdir = save.folder
      obj=methRead(location=out.files[[read.context]],
               sample.id=paste(sample.id,tolower(read.context),sep = "_"),
               assembly=assembly,dbtype=dbtype,pipeline="bismark",header=T, 
               context=read.context,dbdir = dbdir,mincov=mincov)
      obj@sample.id <- sample.id
    }
    else {
      obj=methRead(location=out.files[[read.context]],sample.id=sample.id,
               assembly=assembly,pipeline="bismark",header=T, 
               context=read.context,mincov=mincov)
    }
    if(temp.files ){dummy=lapply(out.files,unlink)}
    
    return(obj)
  }else{return("no object returned")}

})



#' @rdname processBismarkAln-methods
#' @aliases processBismarkAln,list,list,character-method
setMethod("processBismarkAln", signature(location = "list",sample.id="list",
                                         assembly="character"),
          function(location,sample.id,assembly,save.folder,save.context,read.context,
                             nolap,mincov,minqual,phred64,treatment,save.db){
            #check if the given arugments makes sense
            if(length(location) != length(sample.id)){
              stop("length of 'location'  and 'name' should be same\n")
            }
            if( (length(treatment) != length(sample.id)) & (length(treatment) !=0) ){
              stop("length of 'treatment', 'name' and 'location' should be same\n")
            }
            if(!save.db) { 
              # read each given location and record it as methylraw object
              outList=list()
              for(i in 1:length(location))
              {
                data=processBismarkAln(location[[i]],sample.id[[i]],assembly,
                                  save.folder,save.context,read.context,
                                  nolap,mincov,minqual,phred64)# read data
                outList[[i]]=data  
              }
              myobj=new("methylRawList",outList,treatment=treatment)
              
              myobj
            }
            else {
              if(is.null(save.folder)) save.folder=getwd() 
              # read each given location and record it as methylrawDB object
              outList=list()
              for(i in 1:length(location))
              {
                data=processBismarkAln(location[[i]],sample.id[[i]],assembly,
                                  save.folder,save.context,read.context,
                                  nolap,mincov,minqual,phred64,save.db)# read data
                outList[[i]]=data  
              }
              return(new("methylRawListDB",outList,treatment=treatment))
            }
          })

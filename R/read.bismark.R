

#' Function to read in pecent methylation scores from sorted Bismark SAM files
#'
#' The function calls methylation percentage per base from sorted Bismark SAM files. Bismark is a popular aligner for 
#' high-throughput bisulfite sequencing experiments and it outputs its results in SAM format by default. Bismark SAM format contains
#' aligner specific tags which are absolutely necessary for methylation percentage calling. SAM files from other aligners will not work with this function.
#'
#' @param location location of sam file(s). If multiple files are given this arugment must be a list.
#' @param sample.id the id(s) of samples in the same order as file.  If multiple sam files are given this arugment must be a list.
#' @param save.folder The folder which will be used to save methylation call files, if set to NULL no methylation call file will be saved as a text file.
#'                     The files saved can be read into R in less time using \code{read} function in \code{methylKit} 
#' @param save.context A character vector consisting following strings: "CpG","CHG","CHH". The methylation percentages for these methylation contexts will be saved to save.folder
#' @param read.context One of the 'CpG','CHG','CHH' or 'none' strings. Determines what type of methylation context will be read-in to the memory which can be immediately used for analysis.
#'                     If given as 'none', read.bismark will not return any object, but if a save.folder argument given it will save the methylation percentage call files.
#' @param assembly string that determines the genome assembly. Ex: mm9,hg18 etc.
#' @param nolap   if set to TRUE and the SAM file has paired-end reads, the one read of the overlapping paired-end read pair will be ignored for methylation calling.
#' @param mincov  minimum read coverage to call a methylation status for a base.
#' @param minqual minimum phred quality score to call a methylation status for a base.
#' @param phred64 logical ( default: FALSE) you will not need to set this TRUE, Currently bismark gives only phred33 scale
#' @param treatment treatment vector only to be used when location and sample.id parameters are \code{list}s and you are trying to read-in multiple samples that are related to eachother in down-stream analysis. 
#'
#' @return \code{methylRaw} or \code{methylRawList} object
#'
#' @usage read.bismark(location,sample.id,assembly,save.folder=NULL,save.context=c("CpG"),read.context="CpG",nolap=FALSE,mincov=10,minqual=20,phred64=FALSE,treatment)
#'
#' @export
#' @docType methods
#' @rdname read.bismark-methods
#'
#' @examples
#' 
#' # reading one bismark file:
#' my.file=system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit")
#' obj=read.bismark(my.file,"test",assembly="hg18",save.folder=NULL,save.context="CpG",read.context="CpG")
#'  
#' # reading multiple files
#' file.list2=list(system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
#'                system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
#'               system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
#'                system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"))
#' 
#'  objs=read.bismark(location=file.list2
#'              ,sample.id=list("test1","test2","ctrl1","ctrl1"),assembly="hg18",save.folder=NULL,save.context=NULL,read.context="CpG",
#'                                     nolap=FALSE,mincov=10,minqual=20,phred64=FALSE,treatment=c(1,1,0,0))
#' 
setGeneric("read.bismark", function(location,sample.id,assembly,save.folder=NULL,save.context=c("CpG"),read.context="CpG",
                                    nolap=FALSE,mincov=10,minqual=20,phred64=FALSE,treatment) standardGeneric("read.bismark"))

#' @aliases read.bismark,character,character,character-method
#' @rdname read.bismark-methods
setMethod("read.bismark", signature(location = "character",sample.id= "character",assembly= "character"),
                    function(location,sample.id,assembly,save.folder,save.context,read.context,
                             nolap,mincov,minqual,phred64){
                      
                      # check if file exists
                      #if(! file.exists(location) ){stop("File given at the 'location' argument doesn't exist")}
                      
                      # check output types
                      if(! all( save.context %in% c("CpG","CHG","CHH")) ){
                        stop("wrong 'save.context' argument supplied, only 'CpG', 'CHG' or 'CHH' accepted")
                      }
                      
                      # read.type
                       if( !( read.context %in% c("CpG","CHG","CHH","none") ) ){
                        stop("wrong 'read.context' argument supplied, only 'CpG', 'CHG', 'CHH' or 'none' accepted")
                      }
                          
                      # if output.folder NULL create temp
                      # check output folder if not create
                      out.files=list() # list of output files            
                      temp.files=FALSE
                      if(is.null(save.folder)){
                    
                        for(mytype in read.context){
                          out.files[[mytype]]=tempfile(pattern = paste("methylKit_temp",mytype,sep=".") )
                        }
                        temp.files=TRUE # set there are temp files to be deleted
                        
                      }else{
                        try(
                          dir.create(save.folder, showWarnings = TRUE, recursive = TRUE, mode = "0777")
                          )

                        out.files=list()                        
                        for(mytype in unique(c(save.context,read.context)) ){
                          out.files[[mytype]]=paste(save.folder,"/",sample.id,"_",mytype,".txt",sep="") 
                        }
                        
                      }
                      
                      # create the system command accordingly
                      my.opt.str=paste("--read1",location,"--minqual",minqual,"--mincov",mincov,"--type paired_sam")
                      if(phred64){ my.opt.str=paste(my.opt.str,"--phred64") }
                      
                      if("CpG" %in% names(out.files)){my.opt.str=paste(my.opt.str,"--CpG",out.files[["CpG"]] )}
                      if("CHG" %in% names(out.files)){my.opt.str=paste(my.opt.str,"--CHG",out.files[["CHG"]] )}
                      if("CHH" %in% names(out.files)){my.opt.str=paste(my.opt.str,"--CHH",out.files[["CHH"]] )}                     
                      if(nolap){my.opt.str=paste(my.opt.str,"--nolap" )}
                    
                      # get location of the perl script
                      ex.loc=(system.file("exec","methCall.pl", package="methylKit"))
                      if(ex.loc == ""){
                        cmd=paste("perl","/Users/ala2027/Dropbox/PAPERS/R-devel/methylkit/exec/methCall.pl",my.opt.str)
                      }else{
                        cmd=paste("perl",ex.loc,my.opt.str)
                      }
                      #cmd=paste("perl","~/Dropbox\\ Encore/Dropbox/temp/data/methCall.pl",my.opt.str)
                      
                      # then call perl to process the file
                      cat("calling for metylation percentage per base for sample:",sample.id," \n")
                      #cat(cmd)
                      status=try( system(cmd) )
                      
                      if(status != 0){stop("\nError in methylation calling...\nMake sure the file is sorted correctly and it is a legitimate Bismark SAM file\n")}
                      
                      # read the result
                      if(read.context != "none"){
                        cat("Reading methylation percentage per base for sample:",sample.id,"\n\n")
                        obj=read(location=out.files[[read.context]],sample.id=sample.id,assembly=assembly,pipeline="bismark",header=T, context=read.context)
                        if(temp.files ){dummy=lapply(out.files,unlink)}
                        
                        return(obj)
                      }else{return("no object returned")}

})



#' @rdname read.bismark-methods
#' @aliases read.bismark,list,list,character-method
setMethod("read.bismark", signature(location = "list",sample.id="list",assembly="character"),
          function(location,sample.id,assembly,save.folder,save.context,read.context,
                             nolap,mincov,minqual,phred64,treatment){
            #check if the given arugments makes sense
            if(length(location) != length(sample.id)){
              stop("length of 'location'  and 'name' should be same\n")
            }
            if( (length(treatment) != length(sample.id)) & (length(treatment) !=0) ){
              stop("length of 'treatment', 'name' and 'location' should be same\n")
            }
            
            # read each given location and record it as methylraw object
            outList=list()
            for(i in 1:length(location))
            {
              data=read.bismark(location[[i]],sample.id[[i]],assembly,
                                save.folder,save.context,read.context,
                                nolap,mincov,minqual,phred64)# read data
              outList[[i]]=data  
            }
            myobj=new("methylRawList",outList,treatment=treatment)
            
            myobj
          })

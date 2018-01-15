 
#---------------------------------------------------------------------------------------
# regular R functions to be used in S4 functions

# reads gzipped files with data.table::fread
fread.gzipped<-function(filepath,...){
  #require(R.utils)
  #require(data.table)
  
  
  
  # decompress first, fread can't read gzipped files
  if (R.utils::isGzipped(filepath)){
    
    if(.Platform$OS.type == "unix") {
      filepath=paste("gunzip -c",filepath)
    } else {
      filepath <- R.utils::gunzip(filepath,temporary = FALSE, overwrite = TRUE,
                                  remove = FALSE)
    }
    
    
  }
  
  ## Read in the file
  fread(filepath,...)
  
}
# reads a table in a fast way to a dataframe
.readTableFast<-function(filename,header=TRUE,skip=0,sep="auto")
{
  #tab5rows <- read.table(filename, header = header,skip=skip,sep=sep, 
  #                       nrows = 100,stringsAsFactors=FALSE)
  #classes  <- sapply(tab5rows, class)
  #classes[classes=="logical"]="character"
  #lines2keep <- 100
  #nL <- R.utils::countLines(filename)
  #lastrows <- read.table(filename, header = header,skip=nL-lines2keep,sep=sep, 
  #                       stringsAsFactors=FALSE)
  #lastclasses  <- sapply(lastrows, class)
  #if(any(classes!=lastclasses)) {
    #classes[which(classes!=lastclasses)] = lastclasses[which(classes!=lastclasses)]
  #}

  #return( read.table(filename, header = header,skip=skip,sep=sep,colClasses = classes
  #                   )  )
  
  fread.gzipped(filename,header=header,skip=skip,sep=sep,data.table=FALSE)
}

# reformats a data.frame to a standard methylraw data.frame
# no matter what the alignment pipeline
.structureAMPoutput<-function(data,mincov)
{  
  
  # remove data beyond mincoverage
  data = data[data[,5] >= mincov,]
  
  strand=rep("+",nrow(data))
  strand[data[,4]=="R"]="-"
  

  
  numCs=round(data[,5]*data[,6]/100)
  numTs=round(data[,5]*data[,7]/100)
  

  
  
  data.frame(chr=data[,2],start=data[,3],end=data[,3]
             ,strand=strand,coverage=data[,5],numCs=numCs,numTs=numTs)
}

# reformats a generic structure data.frame to a standard methylraw data.frame
# based on the column number assignment and if freqC is fraction or not.
.structureGeneric<-function(data, pipeline,mincov)
{
    fraction=pipeline$fraction
    chr.col=pipeline$chr.col
    start.col=pipeline$start.col
    end.col=pipeline$end.col
    coverage.col=pipeline$coverage.col
    strand.col=pipeline$strand.col
    freqC.col=pipeline$freqC.col
    
    #coerce coverage column to integer
    data[,coverage.col]=round(data[,coverage.col])
    
    # remove data beyond mincoverage
    data = data[data[,coverage.col] >= mincov,]
    
    strand=rep("+",nrow(data))
    strand[data[,strand.col]=="R" | data[,strand.col]=="-"]="-"
    adj=ifelse(fraction, 1, 100)
    numCs=round(data[,coverage.col]*data[,freqC.col]/adj)
    numTs=data[,coverage.col] - numCs
    #id=paste(data[,chr.col], data[,start.col],sep=".")
    
    data.frame(chr=data[,chr.col],start=data[,start.col],end=data[,end.col]
    ,strand=strand,coverage=data[,coverage.col],numCs=numCs,numTs=numTs)
    
}

.check.pipeline.list<-function(pipeline){
    if(!all(c("fraction", "chr.col", "start.col", "end.col", "coverage.col", 
              "strand.col", "freqC.col") %in% names(pipeline))){
        stop("Miss components for pipeline for the generic read.",
             "Try amp, or bismark, or a list in the correct format for",
             "for generic methylation file reading!")
    }
    
    values=c(pipeline$chr.col, pipeline$start.col, pipeline$coverage.col, 
             pipeline$strand.col, pipeline$freqC.col)
    if(any(duplicated(values))){
        stop("Find duplicated column number among chr.col, start.col,", 
             "coverage.col, strand.col, freqC.col!")
    }
}

# checks if dbdir in read-call for methylRawDB and methylRawListDB objects exists
.check.dbdir <- function(dir){
  
  if(dir==getwd() ){
    tabixDir <- paste("methylDB",Sys.Date(),
                      paste(sample(c(0:9, letters, LETTERS),3, replace=TRUE),
                            collapse=""))
    dir.create(tabixDir)
    dir <- paste(dir,"/",tabixDir,collapse = "",sep = "")
    message(paste("creating directory ",getwd(),"/",tabixDir,sep = ""))
  }
  else{
    tempdir <- paste(getwd(),"/",dir,sep = "")
    if(! file.exists(tempdir)){
      message(paste("creating directory ","/",dir,sep = "","..."))
      dir.create(tempdir,recursive = TRUE)
    }
    dir <- tempdir
  }
  return(dir)
}

# unfies forward and reverse strand CpGs on the forward strand if the if
# both are on the same CpG
# if that's the case their values are generally correlated
.CpG.dinuc.unify<-function(cpg)
{
  
  cpgR=cpg[cpg$strand=="-",]
  cpgF=cpg[cpg$strand=="+",]
  cpgR$start=cpgR$start-1L
  cpgR$end=cpgR$end-1L
  cpgR$strand="+"
  
  #cpgR$id=paste(cpgR$chr,cpgR$start,sep=".")
  
  cpgFR=merge(data.table(cpgF),data.table(cpgR),by=c("chr","start","end"))
  #hemi =cpgFR[abs(cpgFR$freqC.x-cpgFR$freqC.y)>=50,]
  #cpgFR=cpgFR[abs(cpgFR$freqC.x-cpgFR$freqC.y)<50,]  
  res=data.frame(
    chr     =as.character(cpgFR$chr),
    start    =as.integer(cpgFR$start),
    end      =as.integer(cpgFR$start),
    strand  =rep("+",nrow(cpgFR)),
    coverage=cpgFR$coverage.x + cpgFR$coverage.y,
    numCs   =cpgFR$numCs.x + cpgFR$numCs.y ,
    numTs   =cpgFR$numTs.x + cpgFR$numTs.y ,stringsAsFactors =FALSE
  )
  Fid=paste(cpgF$chr,cpgF$start,cpgF$end)
  Rid=paste(cpgR$chr,cpgR$start,cpgR$end)
  resid=paste(res$chr,res$start,res$end)  
  res=rbind(res, cpgF[ !  Fid  %in%  resid,],cpgR[ ! Rid  %in%  resid,] )
  #res=res[order(res$chr,res$start),]
  return(res)
}

.CpG.dinuc.unifyOld<-function(cpg)
{
  
  cpgR=cpg[cpg$strand=="-",]
  cpgF=cpg[cpg$strand=="+",]
  cpgR$start=cpgR$start-1
  cpgR$end=cpgR$end-1
  cpgR$strand="+"
  
  #cpgR$id=paste(cpgR$chr,cpgR$start,sep=".")
  
  cpgFR=merge(cpgF,cpgR,by=c("chr","start","end"))
  #hemi =cpgFR[abs(cpgFR$freqC.x-cpgFR$freqC.y)>=50,]
  #cpgFR=cpgFR[abs(cpgFR$freqC.x-cpgFR$freqC.y)<50,]  
  res=data.frame(
    chr     =as.character(cpgFR$chr),
    start    =as.integer(cpgFR$start),
    end      =as.integer(cpgFR$start),
    strand  =rep("+",nrow(cpgFR)),
    coverage=cpgFR$coverage.x + cpgFR$coverage.y,
    numCs   =cpgFR$numCs.x + cpgFR$numCs.y ,
    numTs   =cpgFR$numTs.x + cpgFR$numTs.y ,stringsAsFactors =FALSE
  )
  Fid=paste(cpgF$chr,cpgF$start,cpgF$end)
  Rid=paste(cpgR$chr,cpgR$start,cpgR$end)
  resid=paste(res$chr,res$start,res$end)  
  res=rbind(res, cpgF[ !  Fid  %in%  resid,],cpgR[ ! Rid  %in%  resid,] )
  res=res[order(res$chr,res$start),]
  return(res)
}

#' Read bismark coverage file as a methylKit object
#' 
#' Bismark aligner can output methylation information per base in
#' multiple different formats. This function reads coverage files,
#' which have chr,start,end, number of cytosines (methylated bases) 
#' and number of thymines (unmethylated bases).
#' 
#' @param location a list or vector of file paths to coverage files
#'     
#' @param mincov a numeric value for minimum coverage. Bases that have coverage
#' below this value will be removed.
#' 
#' @return data.frame
.procBismarkCoverage<-function( df,mincov=10)
{

    # remove low coverage stuff
    df=df[ (df[,5]+df[,6]) >= mincov ,]
    
    
    
    
    # make the object (arrange columns of df), put it in a list
     data.frame(chr=df[,1],start=df[,2],end=df[,3],
                strand="*",
                coverage=(df[,5]+df[,6]),numCs=df[,5],
                numTs=df[,6],
                stringsAsFactors = FALSE)
  
}
  

  


#' Read bismark cytosine report file as a methylKit object
#' 
#' Bismark aligner can output methylation information per base in
#' multiple different formats. This function reads cytosine report files,
#' which have chr,start, strand, number of cytosines (methylated bases) 
#' and number of thymines (unmethylated bases),context, trinucletide context.
#' 
#' @param location file path to coverage files
#'     
#' @param mincov a numeric value for minimum coverage. Bases that have coverage
#' below this value will be removed.
#' 
#' @return data.frame
.procBismarkCytosineReport<-function(df,mincov=10){

    # remove low coverage stuff
    df=df[ (df[,4]+df[,5]) >= mincov ,]
    
    
    
    
    # make the object (arrange columns of df), put it in a list
    data.frame(chr=df[,1],start=df[,2],end=df[,2],
                                strand=df[,3],coverage=(df[,4]+df[,5]),
                                numCs=df[,4],numTs=df[,5],
               stringsAsFactors = FALSE)
}

# end of regular functions to be used in S4 functions
#---------------------------------------------------------------------------------------



valid.methylRawObj <- function(object) {
  
    
    data=getData(object)
    check1=( (object@resolution == "base") | (object@resolution == "region") )
    check2=(ncol(data)==7)
    if(check1 & check2){
      return(TRUE)
    }
    else if (! check1 ){
        cat("resolution slot has to be either 'base' or 'region':",
              "other values not allowed")
    }
    else if(! check2){
        cat("data part of methylRaw have",ncol(data),
            "columns, expected 7 columns")
    }

}


#' An S4 class for holding raw methylation data from an alignment pipeline.
#'
#' This object stores the raw mehylation data that is read in through read 
#' function and extends \code{data.frame}.The raw methylation data is basically
#' percent methylation values and read coverage values per genomic base/region.
#'
#' @section Slots:\describe{
#'  \item{\code{sample.id}:}{string for an identifier of the sample}
#'  \item{\code{assembly}:}{string for genome assembly, ex: hg18,hg19,mm9}
#'  \item{\code{context}:}{ methylation context string, ex: CpG,CpH,CHH, etc.}
#'  \item{\code{resolution}:}{ resolution of methylation information, 'base' or 
#'  'region'}
#'                 }
#' @section Details:
#' \code{methylRaw} class extends \code{\link{data.frame}} class therefore 
#' providing novice and experienced R users with a data structure that is well 
#' known and ubiquitous in many R packages.
#' 
#' @section Subsetting:
#'  In the following code snippets, \code{x} is a \code{methylRaw}.
#'  Subsetting by \code{x[i,]} will produce a new object if subsetting is done on
#'  rows. Column subsetting is not directly allowed to prevent errors in the 
#'  downstream analysis. see ?methylKit[ .
#' 
#' @section Accessors:
#' The following functions provides access to data slots of methylRaw:
#' \code{\link[methylKit]{getData}},\code{\link[methylKit]{getAssembly}},
#' \code{\link[methylKit]{getContext}}
#' 
#' @section Coercion:
#'   \code{methylRaw} object can be coerced to 
#'   \code{\link[GenomicRanges]{GRanges}} object via \code{\link{as}} function.
#' 
#' @examples
#' 
#' # example of a raw methylation data as a text file
#' read.table(system.file("extdata", "control1.myCpG.txt", 
#'                                 package = "methylKit"),
#'                                     header=TRUE,nrows=5)
#' 
#' data(methylKit)
#' 
#' # example of a methylRaw object
#' head(methylRawList.obj[[1]])
#' str(head(methylRawList.obj[[1]]))
#' 
#' library(GenomicRanges)
#' 
#' #coercing methylRaw object to GRanges object
#' my.gr=as(methylRawList.obj[[1]],"GRanges")
#' 
#' @name methylRaw-class
#' @aliases methylRaw
#' @docType class
#' @rdname methylRaw-class
#' @export
setClass("methylRaw", contains= "data.frame",representation(
  sample.id = "character", assembly = "character",context="character",
  resolution="character"),validity=valid.methylRawObj)


#' An S4 class for holding a list of methylRaw objects.
#'
#' This class stores the list of  \code{\link[methylKit]{methylRaw}} objects.
#' Functions such as \code{lapply} can be used on this list. It extends
#'  \code{\link[base]{list}} class. This object is primarily produced
#' by \code{\link[methylKit]{read}} function.
#'
#' @section Slots:\describe{
#'                  \item{\code{treatment}}{numeric vector denoting control 
#'                  and test samples}
#'                  \item{\code{.Data}}{a list of 
#'                  \code{\link{methylRaw}} objects  } 
#'                }
#'                
#' @examples
#' data(methylKit)
#' 
#' #applying functions designed for methylRaw on methylRawList object
#' lapply(methylRawList.obj,"getAssembly")
#'
#' @name methylRawList-class
#' @aliases methylRawList
#' @docType class
#' @rdname methylRawList-class
#' @export
setClass("methylRawList", representation(treatment = "numeric"),contains = "list")

#' read file(s) to methylRaw or methylRawList objects
#'
#' The function reads a list of files or single files with methylation 
#' information for bases/region in the genome and creates a methylrawList or 
#' methylraw object. 
#' The information can be stored as flat file database by creating a 
#' methylrawlistDB or methylrawDB object. 
#' @param location file location(s), either a list of locations (each a 
#' character string) or one location string
#' @param sample.id sample.id(s)
#' @param assembly a string that defines the genome assembly such as hg18, mm9. 
#' this is just a string for book keeping. It can be any string. Although,
#' when using multiple files from the same assembly, this string should be 
#' consistent in each object.
#' @param dbtype type of the flat file database, currently only option 
#'        other than NA is "tabix". When "tabix" is given the objects are 
#'        stored in tabix files, which are compressed and indexed. 
#'        The default value is NA, in which case the objects are stored in 
#'        memory.
#' @param header if the input file has a header or not (default: TRUE)
#' @param skip number of lines to skip when reading. Can be set to 1 for bed 
#' files with track line (default: 0)
#' @param sep seperator between fields, same as \code{\link{read.table}} argument 
#' (default: "\\t")
#' @param pipeline name of the alignment pipeline, it can be either "amp", 
#' "bismark","bismarkCoverage", "bismarkCytosineReport" or a list (default:'amp'). 
#' The methylation text files generated from other pipelines can be 
#' read as generic methylation text files by supplying a named 
#' \code{\link[base]{list}} argument as "pipeline" argument.
#' The named \code{list} should containt column numbers which denotes which 
#' column of the text file corresponds to values and genomic location of the 
#' methylation events. See Details for more on possible values for this 
#' argument.
#' @param resolution designates whether methylation information is base-pair 
#' resolution or regional resolution. allowed values 'base' or 'region'. 
#' Default 'base'
#' @param treatment a vector contatining 0 and 1 denoting which samples are 
#' control which samples are test
#' @param context methylation context string, ex: CpG,CpH,CHH, etc. (default:CpG)
#' @param dbdir directory where flat file database(s) should be stored, defaults
#'       to getwd(), working directory.
#' @param mincov minimum read coverage to be included in the methylKit objects.
#'               defaults to 10. Any methylated base/region in the text files
#'               below the mincov value will be ignored.

#' @examples
#' 
#' # this is a list of example files, ships with the package
#' # for your own analysis you will just need to provide set of paths to files
#' # you will not need the "system.file(..."  part
#' file.list=list(
#'          system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
#'          system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
#'          system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#'          system.file("extdata", "control2.myCpG.txt", package = "methylKit") 
#'          )
#'
#' # read the files to a methylRawList object: myobj
#' myobj=methRead(file.list,
#'             sample.id=list("test1","test2","ctrl1","ctrl2"),
#'             assembly="hg18",treatment=c(1,1,0,0))
#'             
#' # read one file as methylRaw object
#' myobj=methRead( file.list[[1]],
#'             sample.id="test1",assembly="hg18")
#'             
#' # read a generic text file containing CpG methylation values
#' # let's first look at the content of the file
#' generic.file=system.file("extdata", "generic1.CpG.txt",package = "methylKit")
#' read.table(generic.file,header=TRUE)
#' 
#' # And this is how you can read that generic file as a methylKit object            
#'  myobj=methRead( generic.file,
#'              pipeline=list(fraction=FALSE,chr.col=1,start.col=2,end.col=2,
#'                            coverage.col=4,strand.col=3,freqC.col=5),
#'              sample.id="test1",assembly="hg18")
#'             
#' # This creates tabix files that save methylation data
#' # Without specified dbdir first creates a folder named the following 
#' # in working directory:
#' # paste("methylDB",Sys.Date(),paste(sample(c(0:9, letters, LETTERS),3,
#' # replace=TRUE),collapse=""))
#' #
#' # Then, saves tabix files from methylKit objects there
#'  myobj=methRead( file.list,
#'                sample.id=list("test1","test2","ctrl1","ctrl2"),
#'                assembly="hg18",treatment=c(1,1,0,0),
#'                dbtype="tabix") 
#' 
#' # This creates a single tabix files that saves methylation data
#' # first creates a "methylDB_objects" directory
#' # Then, saves tabix file from methylKit objects there
#'  myobj=methRead(file.list[[1]],
#'                sample.id="test1",
#'                assembly="hg18",
#'                dbtype="tabix",dbdir="methylDB_objects")
#'            
#'  # tidy up                  
#'  rm(myobj)              
#'  unlink(list.files(pattern = "methylDB",full.names = TRUE),recursive = TRUE)
#'                
#'                
#'  
#'                
#' @section Details:
#'  The output of \code{methRead} is determined by specific input arguments,as there
#'   are \code{location}, \code{sample.id}, \code{assembly} and \code{dbtype}. 
#'  The first three are obligatory, while if the last argument is given database 
#'  features are enabled. 
#'  If then \code{location} refers to an uncompressed file the function will 
#'  create a flat file database and 
#'  the associated methylRawDB object will link to this database. 
#'  If then \code{location} refers to an earlier created database file then the 
#'  object will directly link to this database, 
#'  skipping the preprocessing steps. 
#'  
#'  When \code{pipeline} argument is a list, it is exptected to provide a named 
#'  list with following names.
#'  'fraction' is a logical value, denoting if the column frequency of Cs has a 
#'  range from [0-1] or [0-100]. If true it assumes range is [0-1].
#'  'chr.col" is the number of the column that has chrosome string.   
#'  'start.col' is the number of the column that has start coordinate of the 
#'   base/region of the methylation event.
#'  'end.col'  is the number of the column that has end coordinate of the 
#'   base/region of the methylation event.
#'  'coverage.col' is the number of the column that has read coverage values. 
#'  'strand.col' is the number of the column that has strand information, the 
#'  strand information in the file has to be in the form of '+' or '-', 
#'  'freqC.col' is the number of the column that has the frequency of Cs. 
#'   See examples to see how to read a generic methylation text file.
#'   
#'   Other possible values for  \code{pipeline} argument are 'amp','bismark',
#'   'bismarkCoverage' and 'bismarkCytosineReport'. For 'amp' and 'bismark' 
#'    the function expects a tabular format shown in the webpage 
#'    (http://github.com/al2na/methylKit). 
#'    "amp" and "bismark" expect identical input and are kept for historical 
#'    reasons. 'amp' was a pipeline used in Akalin et al. 2012 Plos Genetics 
#'    paper, publicly available in googlecode. 
#'    
#'    Bismark aligner can output methylation information per
#'    base in
#'    multiple formats. With \code{pipeline='bismarkCoverage'},  
#'    the function reads bismark coverage files,
#'    which have chr,start,end, number of cytosines (methylated bases) 
#'    and number of thymines (unmethylated bases) format.
#'    If bismark coverage files are used the function will not have 
#'    the strand information,so beware of that fact.
#'    With \code{pipeline='bismarkCytosineReport'}, the function expects 
#'    cytosine report files from Bismark,
#'    which have chr,start, strand, number of cytosines (methylated bases) 
#'    , number of thymines (unmethylated bases),context and trinucletide context
#'    format. 
#'    
#'    The function can also read gzipped files. On unix systems, this is achieved
#'    by using \code{zcat filename} and feeding that into \code{data.table::fread}
#'    . On Windows, the file is first uncompressed then read into R using 
#'    \code{data.table::fread}.
#'    
#' @return returns methylRaw, methylRawList, methylRawDB, methylRawListDB object
#' 
#' @export
#' @docType methods
#' @rdname methRead-methods
setGeneric("methRead", function(location,sample.id,assembly,dbtype=NA,
                            pipeline="amp",header=TRUE,skip=0,sep="\t",
                            context="CpG",resolution="base",
                            treatment,dbdir=getwd(),mincov=10) 
  standardGeneric("methRead"))

#' @rdname methRead-methods
#' @aliases methRead,character,character,character-method
setMethod("methRead", signature(location = "character",sample.id="character",
                            assembly="character"),
          
          function(location,sample.id,assembly,dbtype,pipeline,header,skip,
                   sep,context,resolution,dbdir,mincov){ 
    if(! file.exists(location)){
      stop(location,", That file doesn't exist !!!")}
    
    if(  tools::file_ext(location)=="bgz" ) {
      if(!is.na(dbtype)) {
        if(dbtype == "tabix") {
          obj = readMethylRawDB(dbpath = location,dbtype = dbtype, 
                                sample.id = sample.id, 
                                assembly = assembly, context = context, 
                                resolution = resolution)
          return(obj)
        }
      } else {
        stop(paste("file",location,"is compressed,", "\nplease use dbtype='tabix' ",
                   "or provide uncompressed file with supported pipeline."))
      }
    } else {
  
      data<- .readTableFast(location,header=header,skip=skip,sep=sep)# read data  
      if(length(pipeline)==1 ){
        
        if(pipeline %in% c("amp","bismark") )
        {
          data<- .structureAMPoutput(data,mincov)
        }else if(pipeline == "bismarkCytosineReport"){
          data= .procBismarkCytosineReport(data,mincov)
         
        }else if(pipeline == "bismarkCoverage"){
          data= .procBismarkCoverage(data,mincov)
        }
        else{
          stop("unknown 'pipeline' argument, supported processing pipelines are: ",
               "'bismarkCytosineReport','bismarkCoverage','amp' or 'bismark' " )
        }
        
      }
      else{
        .check.pipeline.list(pipeline)
        data<- .structureGeneric(data, pipeline,mincov)
      }
      
      if(is.na(dbtype)){
        obj=new("methylRaw",data,sample.id=sample.id,assembly=assembly,
                context=context,resolution=resolution)
      }else{
        dbdir <- .check.dbdir(dir = dbdir)
        obj=makeMethylRawDB(df=data,dbpath=dbdir,dbtype=dbtype,
                            sample.id=sample.id,assembly=assembly,
                            context=context,resolution=resolution)
      }
    }
    obj 
  }
)

# @param dbtype type of the flat file database, currently only option is "tabix"
# @param dbdir directory where flat file database(s) should be stored, defaults
# @return returns a methylRawListDB object
#' @rdname methRead-methods
#' @aliases methRead,list,list,character-method read
setMethod("methRead", signature(location = "list",sample.id="list",
                            assembly="character"),
          function(location,sample.id,assembly,dbtype,pipeline,header,
                   skip,sep,context,resolution,treatment,dbdir,mincov){ 
            
  #check if the given arugments makes sense
  if(length(location) != length(sample.id)){
    stop("length of 'location'  and 'name' should be same\n")
  }
  if( (length(treatment) != length(sample.id)) & (length(treatment) !=0) ){
    stop("length of 'treatment', 'name' and 'location' should be same\n")
  }
  
  if(all(tools::file_ext(location)=="bgz")) {
    if(!is.na(dbtype)) {
      if(dbtype == "tabix") {
        
        outList=list()
        for(i in 1:length(location))
        {
          
          obj = readMethylRawDB(dbpath = location[[i]],dbtype = dbtype, 
                                sample.id = sample.id[[i]], 
                                assembly = assembly, context = context, 
                                resolution = resolution)
          outList[[i]]=obj
        }
        myobj=new("methylRawListDB",outList,treatment=treatment)
        
        return(myobj)
        
      } else {
        stop("files are compressed,\nplease use dbtype 'tabix' or provide uncompressed file with supported pipeline.")        
        }
    }
  } else if( any(tools::file_ext(location)=="bgz") ) {

    stop("one or more files are compressed,\nplease process uncompressed files beforehand and provide only either only compressed or only uncompressed files.")
    print(location)

  }  else {
  
    if(!is.na(dbtype)) dbdir <- .check.dbdir(dbdir)
    
    # read each given location and record it as methylraw object
    outList=list()
    for(i in 1:length(location))
    {
      data<- .readTableFast(location[[i]],header=header,skip=skip,
                            sep=sep)# read data
      
      if(length(pipeline)==1 )
      {
        if(pipeline %in% c("amp","bismark")){
          data<- .structureAMPoutput(data,mincov)
        }else if(pipeline == "bismarkCytosineReport"){
          data= .procBismarkCytosineReport(data,mincov)
          
        }else if(pipeline == "bismarkCoverage"){
          data= .procBismarkCoverage(data,mincov)
        } else {
          stop("unknown 'pipeline' argument, supported processing pipelines are: ",
               "'bismarkCytosineReport','bismarkCoverage','amp' or 'bismark' ", 
               "\nIf you do not have these formats, please give",
               "\na parameter list containing the format information of", 
               "\nthe data. Please refer details in the function help page")
        }
      }
      else{
        #stop("unknown 'pipeline' argument, supported alignment 
        # pipelines: amp")
        .check.pipeline.list(pipeline)
        data<- .structureGeneric(data, pipeline,mincov)
      }
      
      if(is.na(dbtype)){
        obj=new("methylRaw",data,sample.id=sample.id[[i]],assembly=assembly,
                context=context,resolution=resolution)
      }else{
        obj=makeMethylRawDB(df=data,dbpath=dbdir,dbtype=dbtype,
                            sample.id=sample.id[[i]],
                            assembly=assembly,context=context,
                            resolution=resolution)
      }
      
      outList[[i]]=obj
      
    }
    
    if(is.na(dbtype) ){ 
      myobj=new("methylRawList",outList,treatment=treatment)
    }else{
      myobj=new("methylRawListDB",outList,treatment=treatment)
    }
    
    
    
    myobj
  }
}
)



#' Filter methylRaw, methylRawDB, methylRawList and methylRawListDB object 
#' based on read coverage
#'
#' This function filters \code{methylRaw}, \code{methylRawDB}, 
#' \code{methylRawList} and \code{methylRawListDB} objects.
#' You can filter based on lower read cutoff or high read cutoff.
#'  Higher read cutoff is usefull to eliminate PCR effects
#' Lower read cutoff is usefull for doing better statistical tests.
#'
#' @param methylObj a \code{methylRaw}, \code{methylRawDB}, \code{methylRawList} 
#' or \code{methylRawListDB} object
#' @param lo.count An integer for read counts.Bases/regions having lower 
#' coverage than this count is discarded
#' @param lo.perc  A double [0-100] for percentile of read counts. Bases/regions 
#' having lower coverage than this percentile is discarded
#' @param hi.count An integer for read counts. Bases/regions having higher 
#' coverage than this is count discarded
#' @param hi.perc A double [0-100] for percentile of read counts. Bases/regions 
#' having higher coverage than this percentile is discarded
#' @param chunk.size Number of rows to be taken as a chunk for processing the 
#' \code{methylRawDB} or \code{methylRawListDB} objects, default: 1e6
#' @param save.db A Logical to decide whether the resulting object should be 
#' saved as flat file database or not, default: explained in Details sections  
#' @param ... optional Arguments used when save.db is TRUE
#'            
#'            \code{suffix}
#'                  A character string to append to the name of the output 
#'                  flat file database, 
#'                  only used if save.db is true, default actions: append 
#'                  \dQuote{_filtered} to current filename 
#'                  if database already exists or generate new file with 
#'                  filename \dQuote{sampleID_filtered}
#'                  
#'            \code{dbdir} 
#'                  The directory where flat file database(s) should be stored, 
#'                  defaults
#'                  to getwd(), working directory for newly stored databases
#'                  and to original directory for already existing database
#'                  
#'            \code{dbtype}
#'                  The type of the flat file database, currently only 
#'                  option is "tabix"
#'                  (only used for newly stored databases)
#'                  
#'                   
#' @examples
#' data(methylKit)
#' 
#' # filter out bases with covereage above 500 reads
#' filtered1=filterByCoverage(methylRawList.obj,lo.count=NULL,lo.perc=NULL,
#' hi.count=500,hi.perc=NULL)
#' 
#' # filter out bases with cread coverage above 99.9th percentile of coverage
#' # distribution
#' filtered2=filterByCoverage(methylRawList.obj,lo.count=NULL,lo.perc=NULL,
#' hi.count=NULL,hi.perc=99.9)
#' 
#' # filter out bases with covereage above 500 reads and save to database 
#' # "test1_max500.txt.bgz" 
#' # in directory "methylDB", filtered3 now becomes a \code{methylRawDB} object
#' filtered3=filterByCoverage(methylRawList.obj[[1]], lo.count=NULL,lo.perc=NULL, 
#'                            hi.count=500, hi.perc=NULL, save.db=TRUE, 
#'                            suffix="max500", dbdir="methylDB")
#'                            
#' # tidy up
#' rm(filtered3)
#' unlink("methylDB",recursive=TRUE)
#' 
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with 
#' \code{methylRawDB} or \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects
#'  which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most 
#' systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the 
#' \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects as 
#' \code{methylRawDB} and \code{methylRawListDB}, 
#' while being per default FALSE for \code{methylRaw} and \code{methylRawList}. 
#' If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the database 
#' allows the calculation in-memory, 
#' then you might change the value of this parameter.
#' 
#' 
#' 
#' 
#' @return \code{methylRaw}, \code{methylRawDB}, \code{methylRawList} or
#'  \code{methylRawListDB} object depending on input object
#' @export
#' @docType methods
#' @rdname filterByCoverage-methods
setGeneric("filterByCoverage",function(methylObj,lo.count=NULL,lo.perc=NULL,
                                       hi.count=NULL,hi.perc=NULL,
                                       chunk.size=1e6,save.db=FALSE,...) 
  standardGeneric("filterByCoverage") )

#' @aliases filterByCoverage,methylRaw-method
#' @rdname filterByCoverage-methods
setMethod("filterByCoverage", signature(methylObj="methylRaw"),
                    function(methylObj,lo.count,lo.perc,hi.count,hi.perc,
                             save.db=FALSE,...){
  if( is.null(lo.count) & is.null(lo.perc) & is.null(hi.count) & 
      is.null(hi.perc) ){return(methylObj)}
  
  data=getData(methylObj) # get the data part
  
  #figure out which cut-offs to use, maybe there is more elagent ways, 
  # quick&dirty works for now
  if(is.numeric(lo.count) ){lo.count=lo.count}
  if(is.numeric(lo.perc)){lo.count=quantile(data$coverage,lo.perc/100)}
  if(is.numeric(hi.count)){hi.count=hi.count}
  if(is.numeric(hi.perc)){hi.count=quantile(data$coverage,hi.perc/100)}
  
  if(is.numeric(lo.count)){data=data[data$coverage>=lo.count,]}
  if(is.numeric(hi.count)){data=data[data$coverage<hi.count,]}
  
  
  if(!save.db) {
  
    new("methylRaw",data,sample.id=methylObj@sample.id,
                         assembly=methylObj@assembly,
                         context=methylObj@context,
        resolution=methylObj@resolution)
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
      suffix <- paste0("_","filtered")
    } else { 
      suffix <- args$suffix
      suffix <- paste0("_",suffix)
    }
    
    # create methylRawDB
    #message(paste("creating file",paste0(methylObj@sample.id,suffix,".txt")))
    obj <- makeMethylRawDB(df=data,dbpath=dbdir,dbtype="tabix",
                           sample.id=paste0(methylObj@sample.id,suffix),
                           assembly=methylObj@assembly,
                           context=methylObj@context,
                           resolution=methylObj@resolution)
    obj@sample.id <- methylObj@sample.id
    
    #print(class(obj))
    obj
  }
  
})

#' @aliases filterByCoverage,methylRawList-method
#' @rdname filterByCoverage-methods
setMethod("filterByCoverage", signature(methylObj="methylRawList"),
                    function(methylObj,lo.count,lo.perc,hi.count,hi.perc,
                             save.db=FALSE,...){
  
  if( is.null(lo.count) & is.null(lo.perc) & is.null(hi.count) & 
      is.null(hi.perc) ){
    return(methylObj)}
  
  if(!save.db) {
    new.list=lapply(methylObj,filterByCoverage,lo.count,lo.perc,hi.count,hi.perc)
    new("methylRawList", new.list,treatment=methylObj@treatment)
  } else {
    args <- list(...)
    if( !( "dbdir" %in% names(args)) ){
      dbdir <- .check.dbdir(getwd())
    } else { dbdir <- .check.dbdir(args$dbdir) }
    new.list=lapply(methylObj,filterByCoverage,lo.count,lo.perc,hi.count,
                    hi.perc,save.db=TRUE,dbdir=basename(dbdir),...)
    new("methylRawListDB", new.list,treatment=methylObj@treatment)
  }
})


#' An S4 class for methylation events sampled in multiple experiments
#'
#' This class is designed to contain methylation information such as coverage,
#'  number of methylated bases, etc.. 
#' The methylation events contained in the class must be sampled in multiple 
#' experiments (ex: only CpG bases covered in multiple experiments are stored 
#' in the object of this class).
#' The class extends \code{data.frame} and creates an object that holds 
#' methylation information and genomic location.
#' The object belonging to this class is produced by \code{\link{unite}} function.
#'          
#' @section Slots:\describe{
#'                  \item{\code{sample.ids}:}{character vector for ids of
#'                   samples in the object}
#'
#'                  \item{\code{assembly}:}{name of the genome assembly}
#'
#'                  \item{\code{context}:}{context of methylation. 
#'                  Ex: CpG,CpH,CHH, etc}
#'
#'                  \item{\code{treatment}:}{treatment vector denoting which 
#'                  samples are test and control}
#'
#'                  \item{\code{coverage.index}:}{vector denoting which columns 
#'                  in the data correspons to coverage values}
#'
#'                  \item{\code{numCs.index}:}{vector denoting which columns in 
#'                  the data correspons to number of methylatedCs values}
#'                  \item{\code{numTs.index}:}{vector denoting which columns 
#'                  in the data correspons to number of unmethylated Cs values}
#'                  \item{\code{destranded}:}{ logical value. If \code{TRUE} 
#'                  object is destranded, if \code{FALSE} it is not.}
#'                  \item{\code{resolution}:}{ resolution of methylation 
#'                  information, allowed values: 'base' or 'region'}
#' }
#' 
#' @section Details:
#' \code{methylBase} class extends \code{\link{data.frame}} class 
#' therefore providing novice and experienced R users with a data 
#' structure that is well known and ubiquitous in many R packages.
#' 
#' 
#' @section Subsetting:
#'  In the following code snippets, \code{x} is a \code{methylBase}.
#'  Subsetting by \code{x[i,]} will produce a new object if subsetting is done on
#'  rows. Column subsetting is not directly allowed to prevent errors in the 
#'  downstream analysis. see ?methylKit[ .
#' 
#' @section Accessors:
#' The following functions provides access to data slots of methylBase:
#' \code{\link[methylKit]{getData}},\code{\link[methylKit]{getAssembly}},
#' \code{\link[methylKit]{getContext}}
#' 
#' 
#' @section Coercion:
#'   \code{methylBase} object can be coerced to 
#'   \code{\link[GenomicRanges]{GRanges}} object via \code{\link{as}} function.
#' 
#' 
#' @examples
#' data(methylKit)
#' library(GenomicRanges)
#' my.gr=as(methylBase.obj,"GRanges")
#' 
#' @name methylBase-class
#' @aliases methylBase
#' @docType class
#' @rdname methylBase-class
#' @export
setClass("methylBase",contains="data.frame",representation(
  sample.ids = "character", assembly = "character",context = "character",
  treatment="numeric",coverage.index="numeric",
                                   numCs.index="numeric",numTs.index="numeric",
  destranded="logical",resolution = "character"))


#' unite methylRawList to a single table 
#' 
#' This functions unites \code{methylRawList} and \code{methylRawListDB}
#'  objects that only bases with coverage from all samples are retained.
#' The resulting object is either of class \code{methylBase} or 
#' \code{methylBaseDB} depending on input.
#'
#' @param object a methylRawList or methylRawListDB object to be merged by 
#' common locations covered by reads
#' @param destrand if TRUE, reads covering both strands of a CpG dinucleotide 
#' will be merged, 
#'   do not set to TRUE if not only interested in CpGs (default: FALSE). 
#'   If the methylRawList object
#'   contains regions rather than bases setting destrand to TRUE will have no
#'    effect.
#' @param min.per.group an integer denoting minimum number of samples per 
#' replicate needed to cover a region/base. By default only regions/bases that 
#' are covered in all samples
#' are united as methylBase object, however by supplying an integer for 
#' this argument users can control how many samples needed to cover 
#' region/base to be united as methylBase object.
#' For example, if min.per.group set to 2 and there are 3 replicates per 
#' condition, the bases/regions that are covered in at least 2 replicates 
#' will be united and missing data for uncovered bases/regions will appear 
#' as NAs.
#' @param chunk.size Number of rows to be taken as a chunk for processing the 
#' \code{methylRawListDB} objects, default: 1e6
#' @param mc.cores number of cores to use when processing \code{methylRawListDB}
#'  objects, default: 1, but always 1 for Windows)
#' @param save.db A Logical to decide whether the resulting object should be 
#' saved as flat file database or not, default: explained in Details sections  
#' @param ... optional Arguments used when save.db is TRUE
#'            
#'            \code{suffix}
#'                  A character string to append to the name of the output 
#'                  flat file database, 
#'                  only used if save.db is true, default actions: 
#'                  The default suffix is a 13-character random string appended 
#'                  to the fixed prefix \dQuote{methylBase}, e.g. 
#'                  \dQuote{methylBase_16d3047c1a254.txt.bgz}. 
#'                  
#'            \code{dbdir} 
#'                  The directory where flat file database(s) should be stored, 
#'                  defaults
#'                  to getwd(), working directory for newly stored databases
#'                  and to same directory for already existing database
#'                  
#'            \code{dbtype}
#'                  The type of the flat file database, currently only 
#'                  option is "tabix"
#'                  (only used for newly stored databases)
#' 
#' @return a methylBase or methylBaseDB object depending on input
#' @export
#' @examples
#' 
#'  data(methylKit)
#'  ## Following 
#'  my.methylBase=unite(methylRawList.obj) 
#'  my.methylBase=unite(methylRawList.obj,destrand=TRUE)
#'  
#'  
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with 
#' \code{methylRawDB} or \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized 
#' objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work 
#' for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the 
#' \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects as 
#' \code{methylRawListDB}, 
#' while being per default FALSE for \code{methylRawList}. 
#' If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the database 
#' allows the calculation in-memory, 
#' then you might change the value of this parameter.
#' 
#'   
#' @docType methods
#' @rdname unite-methods
setGeneric("unite", function(object,destrand=FALSE,min.per.group=NULL,
                             chunk.size=1e6,mc.cores=1,save.db=FALSE,...)
  standardGeneric("unite"))

#' @rdname unite-methods
#' @aliases unite,methylRawList-method
setMethod("unite", "methylRawList",
          function(object,destrand,min.per.group,save.db=FALSE,...){
  
  
  
  
  #check if assemblies,contexts and resolutions are same type NOT IMPLEMENTED   
  if( length(unique(vapply(object,function(x) x@context,
                           FUN.VALUE="character"))) > 1)
  {
    stop("supplied methylRawList object have different methylation ",
         "contexts:not all methylation events from the same bases")
  }
  if( length(unique(vapply(object,function(x) x@assembly,
                           FUN.VALUE="character"))) > 1)
  {
    stop("supplied methylRawList object have different genome assemblies")
  }                     
  if( length(unique(vapply(object,function(x) x@resolution
                           ,FUN.VALUE="character"))) > 1)
  {
    stop("supplied methylRawList object have different methylation ",
         "resolutions:some base-pair some regional")
  } 
  
  if( (!is.null(min.per.group)) &  ( ! is.integer( min.per.group ) )  ){
    stop("min.per.group should be an integer\n",
         "try providing integers as 1L, 2L,3L etc.\n")}
  
  #merge raw methylation calls together
  df=getData(object[[1]])
  if(destrand & (object[[1]]@resolution == "base") ){df=.CpG.dinuc.unify(df)}
  df=data.table(df,key=c("chr","start","end","strand"))
  sample.ids=c(object[[1]]@sample.id)
  assemblies=c(object[[1]]@assembly)
  contexts  =c(object[[1]]@context)
  for(i in 2:length(object))
  {
    df2=getData(object[[i]])
    if(destrand & (object[[1]]@resolution == "base") ){df2=.CpG.dinuc.unify(df2)}
    #
    
    if( is.null(min.per.group) ){
      df2=data.table(df2,key=c("chr","start","end","strand"))
      # merge the dat to a data.frame
      df=merge(df,df2,by=c("chr","start","end","strand"),
               suffixes=c(as.character(i-1),as.character(i) ) ) 
      #df=df[df2, nomatch=FALSE]
    }else{
      df2=data.table(df2,key=c("chr","start","end","strand") )
      # using hacked data.table merge called merge2: temporary fix
      df=merge(df,df2,by=c("chr","start","end","strand"),
               suffixes=c(as.character(i-1),as.character(i) ) ,all=TRUE)
      #setkeyv(X,c("chr","start","end","strand"))
      #df=df[df2, nomatch=FALSE]
    }
    sample.ids=c(sample.ids,object[[i]]@sample.id)
    contexts=c(contexts,object[[i]]@context)
  }
  
  # stop if the assembly of object don't match
  if( length( unique(assemblies) ) != 1 ){
    stop("assemblies of methylrawList elements should be same\n")}
  
  
  if(  ! is.null(min.per.group) ){
    # if the the min.per.group argument is supplied, 
    #remove the rows that doesn't have enough coverage
    
    # get indices of coverage,numCs and numTs in the data frame 
    coverage.ind=seq(5,by=3,length.out=length(object))
    numCs.ind   =coverage.ind+1
    numTs.ind   =coverage.ind+2
    start.ind   =2 # will be needed to weed out NA values on chr/start/end/strand
    
    for(i in unique(object@treatment) ){
      my.ind=coverage.ind[object@treatment==i]
      ldat = !is.na(df[,my.ind,with=FALSE])
      if(  is.null(dim(ldat))  ){  # if there is only one dimension
        df=df[ldat>=min.per.group,]
      }else{
        df=df[rowSums(ldat)>=min.per.group,]
      }
    }
    # get all location columns, they are now duplicated with possible NA values
    #mat=df[,c(start.ind-1,start.ind,start.ind+1,start.ind+2)] 
    
    # get location matrix
    #locs=t(apply(mat,1,function(x) unique(x[!is.na(x)]) ) )
    #if(ncol(locs)==3){ # if the resolution is base
    #  df[,c(2:5)]=data.frame(chr=locs[,1],start=as.numeric(locs[,2]),
    # end=as.numeric(locs[,2]),strand=locs[,3])
    #}else{   # if the resolution is region
    #  df[,c(2:5)]=data.frame(chr=locs[,1],start=as.numeric(locs[,2]),
    #              end=as.numeric(locs[,3]),strand=locs[,4])
    #}
    # will be needed to weed out NA values on chr/start/end/strand
    #start.ind   =seq(10,by=7,length.out=length(object)) 
    
    #df=df[,-c(start.ind-1,start.ind,start.ind+1,start.ind+2)]
    #names(df)[2:5]=c("chr","start","end","strand")
  }
  df=as.data.frame(df)
  # get indices of coverage,numCs and numTs in the data frame 
  coverage.ind=seq(5,by=3,length.out=length(object))
  numCs.ind   =coverage.ind+1
  numTs.ind   =coverage.ind+2
  
  # change column names
  names(df)[coverage.ind]=paste(c("coverage"),1:length(object),sep="" )
  names(df)[numCs.ind]   =paste(c("numCs"),1:length(object),sep="" )
  names(df)[numTs.ind]   =paste(c("numTs"),1:length(object),sep="" )
  
  if(!save.db) {
  
    #make methylbase object and return the object
    obj=new("methylBase",(df),sample.ids=sample.ids,
            assembly=unique(assemblies),context=unique(contexts),
            treatment=object@treatment,coverage.index=coverage.ind,
            numCs.index=numCs.ind,numTs.index=numTs.ind,destranded=destrand,
            resolution=object[[1]]@resolution )
    obj
  
  } else {
  
    #print(names(as.list(match.call())))
    # catch additional args 
    args <- list(...)
    # print(args)
    
    if( !( "dbdir" %in% names(args)) ){
      dbdir <- .check.dbdir(getwd())
    } else { dbdir <- .check.dbdir(args$dbdir) }
    if(!( "dbtype" %in% names(args) ) ){
      dbtype <- "tabix"
    } else { dbtype <- args$dbtype }
    if(!( "suffix" %in% names(args) ) ){
      suffix <- NULL
    } else {
      suffix <- args$suffix
      suffix <- paste0("_",suffix)
    }

    # create methylBaseDB
    #message(paste("creating file",paste0(methylObj@sample.id,suffix,".txt")))
    obj <- makeMethylBaseDB(df=df,dbpath=dbdir,dbtype="tabix",
                            sample.ids=sample.ids,
                            assembly=unique(assemblies),
                            context=unique(contexts),
                            treatment=object@treatment,
                            coverage.index=coverage.ind,
                            numCs.index=numCs.ind,
                            numTs.index=numTs.ind,destranded=destrand,
                            resolution=object[[1]]@resolution,
                            suffix=suffix)
    obj@sample.ids <- sample.ids
    
    message(paste0("flatfile located at: ",obj@dbpath))
    
    obj
}
}
)           

#' get correlation between samples in methylBase or methylBaseDB object
#' 
#' The functions returns a matrix of correlation coefficients and/or a set
#'  of scatterplots showing the relationship between samples. The scatterplots
#'  will contain also fitted lines using \code{lm()} for linear regression
#'  and \code{lowess} for polynomial regression. 
#' 
#' @param object a methylBase or methylBaseDB object 
#' @param method a character string indicating which correlation coefficient 
#' (or covariance) is to be computed (default:"pearson", other options are 
#' "kendall" and "spearman") 
#' @param plot scatterPlot if TRUE (default:FALSE) 
#' @param nrow a numeric giving the number of lines to read in of methylBaseDB 
#' object, defaults to 2e6 
#' @return a correlation matrix object and plot scatterPlot
#' @usage getCorrelation(object,method="pearson",plot=FALSE,nrow)
#' @examples
#' 
#' data(methylKit)
#' 
#' getCorrelation(methylBase.obj,method="pearson",plot=FALSE)
#' 
#' # create methylBaseDB
#' methylBaseDB.obj <- unite(methylRawList.obj,save.db=TRUE,dbdir="methylDB")
#' 
#' getCorrelation(methylBaseDB.obj,method="pearson",plot=FALSE,nrow=10000)
#' 
#' 
#' 
#' 
#' # remove Database again
#' rm(methylBaseDB.obj)
#' unlink("methylDB",recursive=TRUE)
#'
#' @section Details: The argument 'nrow' is only evaluated if the 
#' input is a \code{methylBaseDB} object.
#' If 'nrow' is not specified \code{getCorrelation} will read the 
#' first 2M records of the given object,
#' but if you want to read all records 'nrow' has to be NULL. 
#' You should change 'nrow' if using \code{getCorrelation} with 
#' all records of the methylBaseDB object would take too long.
#' 
#' If the scatter plot is plotted, the red line in the plot is from linear 
#' regression fit and the green line is from polynomial regression fit with 
#' \code{stats::lowess}.   
#' 
#' 
#' @export
#' @docType methods
#' @rdname getCorrelation-methods
setGeneric("getCorrelation", function(object,method="pearson",
                                      plot=FALSE,nrow="numeric") 
   standardGeneric("getCorrelation"))

#' @rdname getCorrelation-methods
#' @aliases getCorrelation,methylBase-method
setMethod("getCorrelation", "methylBase",
                    function(object,method,plot){
  meth.mat = getData(object)[, object@numCs.index]/
    (getData(object)[,object@numCs.index] + 
       getData(object)[,object@numTs.index] )                                      
  names(meth.mat)=object@sample.ids
  
  print( cor(meth.mat,method=method) )

  
  panel.cor.pearson <- function(x, y, digits=2, prefix="", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,method="pearson"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }

  panel.cor.kendall <- function(x, y, digits=2, prefix="", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,method="kendall"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  
  panel.cor.spearman <- function(x, y, digits=2, prefix="", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,method="spearman"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  
  
  
  panel.my.smooth2<-function(x, y, col = par("col"), bg = NA, pch = par("pch"),
                             cex = 1, col.smooth = "darkgreen", span = 2/3, 
                             iter = 3, ...) 
  {
       par(new = TRUE)    #par(usr = c(usr[1:2], 0, 1.5) )
      smoothScatter(x, y,colramp=colorRampPalette(topo.colors(100)), bg = bg)
      ok <- is.finite(x) & is.finite(y)
      if (any(ok)) 
          lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
                col = col.smooth, ...)
          abline(lm(y[ok]~x[ok]), col="red")
  }
  
  panel.my.smooth<-function(x, y, col = par("col"), bg = NA, pch = par("pch"), 
                            cex = 0.3, col.smooth = "green", 
                            span = 2/3, iter = 3, ...) 
  {
      points(x, y, pch = 20, col = densCols(x,y,
                             colramp=colorRampPalette(topo.colors(20))), 
             bg = bg, cex = 0.1)
      ok <- is.finite(x) & is.finite(y)
      if (any(ok)){
          lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
                col = col.smooth, ...);
          abline(lm(y[ok]~x[ok]), col="red")}
  }
  panel.hist <- function(x, ...)
  {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(usr[1:2], 0, 1.5) )
      h <- hist(x, plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
  }
  
  if(plot)
  {  
   
    if(method=="spearman")
    { pairs(meth.mat, 
        lower.panel=panel.my.smooth2, 
        upper.panel=panel.cor.spearman,
        diag.panel=panel.hist,main=paste(object@context, 
                                         object@resolution ,method,"cor.") )
    }
    if(method=="kendall")
    { pairs(meth.mat, 
            lower.panel=panel.my.smooth2, 
            upper.panel=panel.cor.kendall,
            diag.panel=panel.hist,main=paste(object@context, 
                                             object@resolution ,method,"cor.") )
    }
    if(method=="pearson")
    { pairs(meth.mat, 
            lower.panel=panel.my.smooth2, 
            upper.panel=panel.cor.pearson,
            diag.panel=panel.hist,main=paste(object@context, 
                                             object@resolution ,method,"cor.") )
    }
    
    
  }
}  
 )



#' get coverage stats from methylRaw object
#' 
#' The function returns basic statistics about read coverage per base. 
#' It can also plot a histogram of read coverage values.
#' 
#' @param object a \code{methylRaw} or \code{methylRawDB} object 
#' @param plot plot a histogram of coverage if TRUE (default:FALSE) 
#' @param both.strands do stats and plot for both strands if TRUE 
#' (default:FALSE)
#' @param labels should the bars of the histrogram have labels showing 
#' the percentage of values in each bin (default:TRUE)
#' @param ... options to be passed to \code{\link[graphics]{hist}} function
#' @param chunk.size Number of rows to be taken as a chunk for processing the 
#' \code{methylRawDB} objects (default: 1e6)
#' @examples
#' data(methylKit)
#' 
#' # gets coverage stats for the first sample in methylRawList.obj object
#' getCoverageStats(methylRawList.obj[[1]],plot=TRUE,
#' both.strands=FALSE,labels=TRUE)
#' 
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with
#'  \code{methylRawDB} or \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects
#'  which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most 
#' systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the 
#' \code{chunk.size}.
#' 
#' @return a summary of coverage statistics or plot a histogram of coverage
#' @export
#' @docType methods
#' @rdname getCoverageStats-methods
setGeneric("getCoverageStats", function(object,plot=FALSE,both.strands=FALSE,
                                        labels=TRUE,...,chunk.size=1e6) 
  standardGeneric("getCoverageStats"))

#' @rdname getCoverageStats-methods
#' @aliases getCoverageStats,methylRaw-method
setMethod("getCoverageStats", "methylRaw",
                    function(object,plot,both.strands,labels,...){
                      
  if(!plot){
    qts=seq(0,0.9,0.1) # get quantiles
    qts=c(qts,0.95,0.99,0.995,0.999,1)                          
    
    if(both.strands){       
      plus.cov=object[object$strand=="+",]$coverage
      mnus.cov=object[object$strand=="-",]$coverage
      
      cat("read coverage statistics per base\n\n")
      cat("FORWARD STRAND:\n")
      cat("summary:\n")
      print( summary( plus.cov ) )
      cat("percentiles:\n")
      print(quantile( plus.cov,p=qts ))
      cat("\n\n")
      cat("REVERSE STRAND:\n")
      cat("summary:\n")
      print( summary( mnus.cov ) )
      cat("percentiles:\n")
      print(quantile( mnus.cov,p=qts ))
      cat("\n")                          
    }else{
      
      all.cov=object$coverage
      
      cat("read coverage statistics per base\n")
      cat("summary:\n")
      print( summary( all.cov ) )
      cat("percentiles:\n")
      print(quantile( all.cov,p=qts ))
      cat("\n")
    }
    
  }else{
    if(both.strands){   
      plus.cov=object[object$strand=="+",]$coverage
      mnus.cov=object[object$strand=="-",]$coverage
      
      par(mfrow=c(1,2))
      if(labels){
        a=hist(log10(plus.cov),plot=FALSE)
        my.labs=as.character(round(100*a$counts/length(plus.cov),1))
      }else{my.labs=FALSE}

      hist(log10(plus.cov),col="chartreuse4",
           xlab=paste("log10 of read coverage per",object@resolution),
           main=paste("Histogram of", object@context, 
                      "coverage: Forward strand"),
           labels=my.labs,...)
      mtext(object@sample.id, side = 3)

     if(labels){
        a=hist(log10(mnus.cov),plot=FALSE)
        my.labs=as.character(round(100*a$counts/length(mnus.cov),1))
      }else{my.labs=FALSE}
      a=hist(log10(mnus.cov),plot=FALSE)
      hist(log10(mnus.cov),col="chartreuse4",
           xlab=paste("log10 of read coverage per",object@resolution),
           main=paste("Histogram of", object@context, 
                      "coverage: Reverse strand"),
           labels=my.labs,...)
      mtext(object@sample.id, side = 3)

    }else{
      all.cov= object$coverage
     if(labels){
       a=hist(log10(all.cov),plot=FALSE)
       my.labs=as.character(round(100*a$counts/length(all.cov),1))
      }else{my.labs=FALSE}                          

      hist(log10(all.cov),col="chartreuse4",
           xlab=paste("log10 of read coverage per",object@resolution),
           main=paste("Histogram of", object@context, "coverage"),
           labels=my.labs,...)
      mtext(object@sample.id, side = 3)

    }
    
    
  }
    
                      
})

#' get Methylation stats from methylRaw or methylRawDB object
#' 
#' The function returns basic statistics about % methylation per base/region. 
#' It can also plot a histogram of % methylation values.
#' 
#' @param object a \code{methylRaw} or \code{methylRawDB} object 
#' @param plot plot a histogram of Methylation if TRUE (deafult:FALSE) 
#' @param both.strands do plots and stats for both strands seperately  if
#' TRUE (deafult:FALSE)
#' @param labels should the bars of the histrogram have labels showing the 
#' percentage of values in each bin (default:TRUE)
#' @param ... options to be passed to \code{\link[graphics]{hist}} function.
#' @param chunk.size Number of rows to be taken as a chunk for processing the
#'  \code{methylRawDB} objects (default: 1e6)
#' @examples
#' data(methylKit)
#' 
#' # gets Methylation stats for the first sample in methylRawList.obj object
#' getMethylationStats(methylRawList.obj[[1]],plot=TRUE,
#' both.strands=FALSE,labels=TRUE)
#'
#'@section Details:
#' The parameter \code{chunk.size} is only used when working with 
#' \code{methylRawDB} or \code{methylRawListDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized 
#' objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most 
#' systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the
#'  \code{chunk.size}.
#' 
#' @return a summary of Methylation statistics or plot a histogram of coverage
#' @export
#' @docType methods
#' @rdname getMethylationStats-methods
setGeneric("getMethylationStats", function(object,plot=FALSE,both.strands=FALSE,
                                           labels=TRUE,...,chunk.size=1e6) 
  standardGeneric("getMethylationStats"))

#' @rdname getMethylationStats-methods
#' @aliases getMethylationStats,methylRaw-method
setMethod("getMethylationStats", "methylRaw",
                    function(object,plot,both.strands,labels,...){
                      
  plus.met=100* object[object$strand=="+",]$numCs/object[object$strand=="+",]$coverage
  mnus.met=100* object[object$strand=="-",]$numCs/object[object$strand=="-",]$coverage
  all.met =100* object$numCs/object$coverage
  
  if(!plot){
    qts=seq(0,0.9,0.1) # get quantiles
    qts=c(qts,0.95,0.99,0.995,0.999,1)                          
    
    if(both.strands){       
      
      cat("methylation statistics per base\n\n")
      cat("FORWARD STRAND:\n")
      cat("summary:\n")
      print( summary( plus.met ) )
      cat("percentiles:\n")
      print(quantile( plus.met,p=qts ))
      cat("\n\n")
      cat("REVERSE STRAND:\n")
      cat("summary:\n")
      print( summary( mnus.met ) )
      cat("percentiles:\n")
      print(quantile( mnus.met,p=qts ))
      cat("\n")                          
    }else{
      
      
      cat("methylation statistics per base\n")
      cat("summary:\n")
      print( summary( all.met ) )
      cat("percentiles:\n")
      print(quantile( all.met,p=qts ))
      cat("\n")
    }
    
  }else{
    if(both.strands){   
      
      par(mfrow=c(1,2))
      if(labels){
        a=hist((plus.met),plot=FALSE,...)
        my.labs=as.character(round(100*a$counts/length(plus.met),1))
      }else{my.labs=FALSE}
      hist((plus.met),col="cornflowerblue",
           xlab=paste("% methylation per",object@resolution),
           main=paste("Histogram of %", object@context,"methylation: Forward strand"),
           labels=my.labs,...)
      mtext(object@sample.id, side = 3)

      if(labels){                          
        a=hist((mnus.met),plot=FALSE,...)
        my.labs=as.character(round(100*a$counts/length(mnus.met),1))
      }
      else{my.labs=FALSE}
      
      hist((mnus.met),col="cornflowerblue",
           xlab=paste("% methylation per",object@resolution),
           main=paste("Histogram of %", object@context,"methylation: Reverse strand"),
           labels=my.labs,...)
      mtext(object@sample.id, side = 3)

    }else{
      if(labels){                          
      
        a=hist((all.met),plot=FALSE,...)
        my.labs=as.character(round(100*a$counts/length(all.met),1))
      }else{my.labs=FALSE}
      hist((all.met),col="cornflowerblue",
           xlab=paste("% methylation per",object@resolution),
           main=paste("Histogram of %", object@context,"methylation"),
           labels=my.labs,...)
      mtext(object@sample.id, side = 3)

    }
    
    
  }
    
                      
})





# get distribution of difference between samples in methylBase object
# unites methylrawlist objects based on chromosomal positions of CpG dinucleotides
#setGeneric("getDifference", function(object,plot=F) standardGeneric("getCorrelation"))
#setMethod("getDifference", "methylBase",
#                    function(object,plot){
#                        meth.mat = object@data[, object@numCs.index]/(object@data[,object@numCs.index] + object@data[,object@numTs.index] )                                      
#                        names(meth.mat)=object@sample.ids                      
#                      
#                        ind=t(combn(1:4,2) )
#                        d.list=(meth.mat[,ind[1,1]]-meth.mat[,ind[1,2]])
#                        for(i in 2:nrow(ind))
#                        {
#                          d.list=cbind(d.list,meth.mat[,ind[i,1]]-meth.mat[,ind[i,2]])
#                        }
#                    }                       
#                  ) 



                        
#
#  methylBase accessor and show functions
#


#' show method for methylKit classes
#' 
#' The show method works for \code{methylRaw},\code{methylRawDB},
#' \code{methylRawList},\code{methylRawListDB},
#' \code{methylBase},\code{methylBaseDB} and \code{methylDiff} objects
#' 
#' @param object any methylKit object
#' 
#' @examples
#' data(methylKit)
#' methylDiff.obj
#' show(methylDiff.obj)
#' 
#' 
#'
#' @rdname show-methods
#' @aliases show,methylBase
setMethod("show", "methylBase", function(object) {
  
  cat("methylBase object with",nrow(object),"rows\n--------------\n")
  print(head(object))
  cat("--------------\n")
  cat("sample.ids:",object@sample.ids,"\n")
  cat("destranded",object@destranded,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("treament:", object@treatment,"\n")
  cat("resolution:", object@resolution,"\n")
})

#' @rdname show-methods
#' @aliases show,methylRaw
setMethod("show", "methylRaw", function(object) {
  
  cat("methylRaw object with",nrow(object),"rows\n--------------\n")
  print(head(object))
  cat("--------------\n")
  cat("sample.id:",object@sample.id,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("resolution:", object@resolution,"\n\n")
})

#' @rdname show-methods
#' @aliases show,methylRawList
setMethod("show", "methylRawList", function(object) {
  
  cat("methylRawList object with",length(object),"methylRaw objects\n\n")
  
  lapply(object,show)
  cat("treatment:", object@treatment,"\n")
  
})


#' get assembly of the genome
#' 
#' The function returns the genome assembly stored in any of the 
#' \code{\link{methylBase}},\code{\link{methylBaseDB}},\code{\link{methylRaw}},
#' \code{\link{methylRawDB}},\code{\link{methylDiff}} objects
#' 
#' @param x an \code{\link{methylBase}},\code{\link{methylBaseDB}},
#' \code{\link{methylRaw}},\code{\link{methylRawDB}} or \code{\link{methylDiff}} object
#' @usage getAssembly(x)
#' @examples
#' 
#' data(methylKit)
#' 
#' getAssembly(methylBase.obj)
#' getAssembly(methylDiff.obj)
#' getAssembly(methylRawList.obj[[1]])
#' 
#' 
#' @return the assembly string for the object
#' @export
#' @docType methods
#' @rdname getAssembly-methods
setGeneric("getAssembly", def=function(x) standardGeneric("getAssembly"))

#' @rdname getAssembly-methods
#' @aliases getAssembly,methylBase-method
setMethod("getAssembly", signature="methylBase", definition=function(x) {
                return(x@assembly)
        }) 

#' @rdname getAssembly-methods
#' @aliases getAssembly,methylRaw-method
setMethod("getAssembly", signature="methylRaw", definition=function(x) {
                return(x@assembly)
        })

#' get the context of methylation
#' 
#' The function returns the context of methylation. For example: "CpG","CHH" or "CHG"
#' 
#' @param x an \code{\link{methylBase}},\code{\link{methylBaseDB}},
#' \code{\link{methylRaw}},\code{\link{methylRawDB}} or an 
#' \code{\link{methylDiff}} object
#' @usage getContext(x)
#' @examples 
#' 
#' data(methylKit)
#' 
#' getContext(methylBase.obj)
#' getContext(methylDiff.obj)
#' getContext(methylRawList.obj[[1]])
#'
#' @return a string for the context methylation 
#' @export
#' @docType methods
#' @rdname getContext-methods
setGeneric("getContext", def=function(x) standardGeneric("getContext"))

#' @rdname getContext-methods
#' @aliases getContext,methylBase-method
setMethod("getContext", signature="methylBase", definition=function(x) {
                return(x@context)
        })

#' @rdname getContext-methods
#' @aliases getContext,methylRaw-method
setMethod("getContext", signature="methylRaw", definition=function(x) {
                return(x@context)
        })



#' get the data slot from the methylKit objects
#' 
#' The functions retrieves the table containing methylation information from 
#' \code{methylKit} Objects.
#' The data retrived from this function is of a \code{\link{data.frame}}. 
#' This is basically containing all relevant methylation information per 
#' genomic region or base.
#'
#' @param x an \code{\link{methylBase}},\code{\link{methylBaseDB}},
#' \code{\link{methylRaw}},\code{\link{methylRawDB}} or 
#' \code{\link{methylDiff}} object
#' @usage getData(x)
#' @examples
#' data(methylKit)
#' 
#' # following commands show first lines of returned 
#' # data.frames from getData() function
#' head(
#' getData(methylBase.obj)
#' )
#' 
#' head( getData(methylDiff.obj))
#'
#' head(getData(methylRawList.obj[[1]]))
#' 
#' 
#' @return data frame for methylation events
#' @export
#' @docType methods
#' @rdname getData-methods
setGeneric("getData", def=function(x) standardGeneric("getData"))

#' @rdname getData-methods
#' @aliases getData,methylBase-method
setMethod("getData", signature="methylBase", definition=function(x) {
                #return(as(x,"data.frame"))
                return(S3Part(x, strictS3 = TRUE))
}) 

#' @rdname getData-methods
#' @aliases getData,methylRaw-method
setMethod("getData", signature="methylRaw", definition=function(x) {
                #return(as(x,"data.frame"))
                return(S3Part(x, strictS3 = TRUE))
})

## CONVERTOR FUNCTIONS FOR methylRaw/methylRawDB and methylBase/methylBaseDB OBJECT
#convert methylRaw to GRanges
setAs("methylRaw", "GRanges", function(from)
                      {
                        from2=getData(from)
                        GRanges(seqnames=as.character(from2$chr),
                                ranges=IRanges(start=from2$start, end=from2$end),
                                       strand=from2$strand, 
                                       coverage=from2$coverage,
                                       numCs   =from2$numCs,
                                       numTs  =from2$numTs                                
                                       )

})

setAs("methylBase", "GRanges", function(from)
                      {
  from=getData(from)
  GRanges(seqnames=as.character(from$chr),
          ranges=IRanges(start=from$start, end=from$end),
                 strand=from$strand, 
                 data.frame(from[,5:ncol(from)])
                 )

})

### subset methylBase and methylRaw objects

#' selects rows from of methylKit objects
#'
#' The function returns a subset of data contained in the \code{methylKit} 
#' objects.
#' 
#' @param x an \code{\link{methylBase}},\code{\link{methylBaseDB}},
#' \code{\link{methylRaw}},\code{\link{methylRawDB}} or
#'  \code{\link{methylDiff}} object
#' @param i a numeric or logical vector. This vector corresponds to bases or 
#'          regions contained in \code{methylKit} objects.The vector is used to 
#'          subset the data.  
#' @usage select(x,i)
#' @examples
#' data(methylKit)
#' 
#' 
#' methylRawDB.obj=methRead( system.file("extdata","test1.txt.bgz",package="methylKit"),
#'                           sample.id="test1", assembly="hg18",
#'                           dbtype = "tabix",dbdir = "methylDB")
#'
#' methylBaseDB.obj=unite(methylRawList.obj,save.db=TRUE,dbdir="methylDB")
#' 
#' 
#'  # selects first hundred rows, returns a methylRaw object
#' subset1=select(methylRawList.obj[[1]],1:100)
#' subset1=select(methylRawDB.obj,1:100)
#' 
#' # selects first hundred rows, returns a methylBase object
#' subset2=select(methylBase.obj,1:100) 
#' subset2=select(methylBaseDB.obj,1:100)
#' 
#' # selects first hundred rows, returns a methylDiff object
#' subset3=select(methylDiff.obj,1:100)
#' 
#' 
#' 
#' 
#' # remove Database again
#' rm(methylBaseDB.obj)
#' rm(methylRawDB.obj)
#' unlink("methylDB",recursive=TRUE)
#' 
#' @return a \code{\link{methylBase}},\code{\link{methylRaw}} or 
#'           \code{\link{methylDiff}} object depending on the input object.
#' @export
#' @docType methods
#' @rdname select-methods
setGeneric("select", def=function(x,i) standardGeneric("select"))


#' @aliases select,methylBase-method
#' @rdname select-methods
setMethod("select", "methylBase",
          function(x, i)
          {
            if( max(i) > nrow(x)  )
              stop("subscript contains out-of-bounds indices")

            new("methylBase",getData(x)[i,],
                sample.ids = x@sample.ids, 
                assembly = x@assembly,
                context = x@context,
                treatment=x@treatment,
                coverage.index=x@coverage.index,
                numCs.index=x@numCs.index,
                numTs.index=x@numTs.index,
                destranded=x@destranded,
                resolution =x@resolution)
           }
)


#' @aliases select,methylRaw-method
#' @rdname select-methods
setMethod("select", "methylRaw",
          function(x, i)
          {
            if( max(i) > nrow(x)  )
              stop("subscript contains out-of-bounds indices")

            new("methylRaw",getData(x)[i,],
                sample.id=x@sample.id,
                assembly=x@assembly,
                context=x@context,
                resolution=x@resolution)
           }
          

)

#' extract parts of methylRaw,methylRawDB,methylBase,methylBaseDB and methylDiff data
#' 
#' The function extracts part of the data and returns a new object.
#' @name extract
#' @param x an \code{\link{methylBase}},\code{\link{methylBaseDB}},
#' \code{\link{methylRaw}},\code{\link{methylRawDB}} or 
#'          \code{\link{methylDiff}} object
#' @param i a numeric or logical vector. This vector corresponds to bases or 
#'          regions contained in \code{methylKit} objects.The vector is used to 
#'          subset the data.
#' @param j This argument can not be used for the extraction of columns.
#'          As unintentional extraction of the columns will cause an error in
#'          the downstream analysis. Using this argument will cause an error.
#'           Use \code{\link[methylKit]{getData}} to access the data part of 
#'           the objects. 
#' 
#'
#'        
#'        
#' @examples
#' data(methylKit)
#' 
#' # selects first hundred rows, returns a methylRaw object
#' subset1=methylRawList.obj[[1]][1:100] 
#' 
#' # selects first hundred rows, returns a methylBase object
#' subset2=methylBase.obj[1:100,] 
#' 
#' # selects first hundred rows, returns a methylDiff object
#' subset3=methylDiff.obj[1:100,] 
#' 
#' # This will get chromomsomes, will return a factor
#' # That means the resulting object will ceases to be a methylKit object
#' chrs=methylDiff.obj[[2]]
#' 
#' 
#' @docType methods
#' @rdname extract-methods
NULL

#' @aliases [,methylRaw,ANY,ANY,ANY-method
#' @aliases extract,methylRaw,ANY-method
#' @rdname extract-methods
setMethod("[", signature(x="methylRaw", i = "ANY", j="ANY"),  
          function(x,i,j){
            #cat(missing(i),"\n",missing(j),"\n",missing(drop))
            if(!missing(j)){
              stop(paste("subsetting on columns is not allowed for",class(x),
                         "object\nif you want to do extraction on the data part", 
                          "of the object use getData() first"),
                   call. = FALSE)
            }
            select(x,i)
          }
              )

#' @aliases [,methylBase,ANY,ANY,ANY-method
#' @aliases extract,methylBase,ANY-method
#' @rdname extract-methods
setMethod("[",signature(x="methylBase", i = "ANY", j="ANY"), 
          function(x,i,j){
            #cat(missing(i),"\n",missing(j),"\n",missing(drop))
            if(!missing(j)){
              stop(paste("subsetting on columns is not allowed for",class(x),
                         "object\nif you want to do extraction on the data part", 
                          "of the object use getData() first"),
                   call. = FALSE)
            }
            select(x,i)
          }
)

#' selects records of methylDB objects lying inside a GRanges range
#'
#' The function selects records from a \code{\link{methylBaseDB}}, 
#' \code{\link{methylRawDB}} or \code{\link{methylDiffDB}} object 
#' that lie inside the regions given by \code{ranges} of class \code{GRanges} 
#' and returns an object of class 
#' \code{\link{methylBase}}, \code{\link{methylRaw}} or \code{\link{methylDiff}} 
#' 
#' @param object an \code{\link{methylBaseDB}},\code{\link{methylRawDB}} or \code{\link{methylDiffDB}} object
#' @param ranges a GRanges object specifying the regions of interest
#' 
#' @usage selectByOverlap(object,ranges)
#' @examples
#' data(methylKit)
#' 
#' file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
#'                 system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
#'                 system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#'                 system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
#' 
#' methylRawListDB.obj=methRead(file.list,
#'                          sample.id=list("test1","test2","ctrl1","ctrl2"),
#'                          assembly="hg18",treatment=c(1,1,0,0),
#'                          dbtype = "tabix",dbdir = "methylDB")
#'
#' methylBaseDB.obj=unite(methylRawListDB.obj)
#'
#' methylDiffDB.obj = calculateDiffMeth(methylBaseDB.obj)
#' 
#' # define the windows of interest as a GRanges object, this can be any set 
#' # of genomic locations
#' library(GenomicRanges)
#' my.win=GRanges(seqnames="chr21",
#' ranges=IRanges(start=seq(from=9764513,by=10000,length.out=20),width=5000) )
#' 
#' # selects the records that lie inside the regions
#' myRaw <- selectByOverlap(methylRawListDB.obj[[1]],my.win)
#' 
#' # selects the records that lie inside the regions
#' myBase <- selectByOverlap(methylBaseDB.obj,my.win)
#' 
#' # selects the records that lie inside the regions
#' myDiff <- selectByOverlap(methylDiffDB.obj,my.win)
#' 
#' 
#' rm(methylRawListDB.obj)
#' rm(methylBaseDB.obj)
#' rm(methylDiffDB.obj)
#' unlink("methylDB",recursive=TRUE)
#' 
#' @return a \code{\link{methylBase}},\code{\link{methylRaw}} or 
#'           \code{\link{methylDiff}} object depending on the input object.
#'           
#' @author Alexander Gosdschan           
#' @export
#' @docType methods
#' @rdname selectByOverlap-methods
setGeneric("selectByOverlap", def=function(object,ranges) 
  standardGeneric("selectByOverlap"))

#' @aliases selectByOverlap,methylRaw-method
#' @rdname selectByOverlap-methods
setMethod("selectByOverlap", "methylRaw",
          function(object, ranges){
            
            if(missing(ranges) | class(ranges)!="GRanges") {
              stop("No ranges specified or given ranges object not of class ", 
                   "GRanges, please check your input!")            }
            hits <- findOverlaps(ranges,as(object,"GRanges"))@subjectHits
            
            return(object[hits])
          }
)

#' @aliases selectByOverlap,methylRawList-method
#' @rdname selectByOverlap-methods
setMethod("selectByOverlap", "methylRawList",
          function(object, ranges){
            
            if(missing(ranges) | class(ranges)!="GRanges") {
              stop("No ranges specified or given ranges object not of class ", 
                   "GRanges, please check your input!")
            }
            
            new.list <- lapply(object,selectByOverlap,ranges)
            
            new("methylRawList",new.list,treatment=object@treatment)
            
          }
)

#' @aliases selectByOverlap,methylBase-method
#' @rdname selectByOverlap-methods
setMethod("selectByOverlap", "methylBase",
          function(object, ranges){
            
            if(missing(ranges) | class(ranges)!="GRanges") {
              stop("No ranges specified or given ranges object not of class",
                   "GRanges, please check your input!")
            }
            hits <- findOverlaps(ranges,as(object,"GRanges"))@subjectHits
            
            return(object[hits])
          }
)




#' Get or Set treatment vector of methylKit object
#' 
#' The function returns or replaces the treatment vector stored in any of the 
#' following methylKit objects:
#' \code{\link{methylBase}},\code{\link{methylRawList}},\code{\link{methylBaseDB}},
#' \code{\link{methylRawListDB}},\code{\link{methylDiff}},\code{\link{methylDiffDB}}.
#'  
#' 
#' @param x a \code{methylKit} object
#' @param value a valid replacement for the treatment vector of the object
#' @usage 
#' getTreatment(x)
#' getTreatment(x) <- value
#' @examples
#' 
#' data(methylKit)
#' 
#' # The treatment vector can be printed ..
#' getTreatment(methylBase.obj)
#'  
#' # .. or replaced with a new one  
#' newObj <- methylBase.obj
#' getTreatment(newObj) <- c(1,2,3,4)
#' getTreatment(newObj)
#' 
#' 
#' @export
#' @docType methods
#' @rdname getTreatment-methods
setGeneric("getTreatment", def=function(x) standardGeneric("getTreatment"))
#' @export 
#' @rdname getTreatment-methods
setGeneric("getTreatment<-", def=function(x, value="numeric") {
  standardGeneric("getTreatment<-")})

#' @rdname getTreatment-methods
#' @aliases getTreatment,methylRawList-method
setMethod("getTreatment", signature = "methylRawList", function(x) {
  return(x@treatment)
})

#' @rdname getTreatment-methods
#' @aliases getTreatment,methylRawList-method
setReplaceMethod("getTreatment", signature = "methylRawList", function(x, value) {
  
  if(! ( length(x@treatment) == length(value) ) ){
    stop("The new treatment vector is not valid, check the length of input")
  } else {
    x@treatment <- value
    return(x)
  }
})


#' @rdname getTreatment-methods
#' @aliases getTreatment,methylBase-method
setMethod("getTreatment", signature = "methylBase", function(x) {
  return(x@treatment)
})

#' @rdname getTreatment-methods
#' @aliases getTreatment,methylBase-method
setReplaceMethod("getTreatment", signature = "methylBase", function(x, value) {
  
  if(! ( length(x@treatment) == length(value) ) ){
    stop("The new treatment vector is not valid, check the length of input")
  } else {
    x@treatment <- value
    return(x)
  }
  
})




#' Get or Set Sample-IDs of the methylKit objects
#' 
#' The function returns or replaces the sample-ids stored in any of the 
#' following \code{methylKit} objects:
#' \code{\link{methylRaw}}, \code{\link{methylRawDB}}, \code{\link{methylBase}}, 
#' \code{\link{methylBaseDB}},
#' \code{\link{methylRawList}}, \code{\link{methylRawListDB}}, 
#' \code{\link{methylDiff}}, \code{\link{methylDiffDB}}.
#' 
#' @param x an \code{\link{methylBaseDB}},\code{\link{methylRawListDB}} or 
#' \code{\link{methylDiffDB}} object
#' @param value a valid replacement vector for the sample-ids of the object 
#' @usage 
#' getSampleID(x)
#' getSampleID(x) <- value
#' @examples
#' 
#' data(methylKit)
#' 
#' #The Sample-Ids can be printed ..
#' getSampleID(methylBase.obj)
#' 
#' # .. or replaced. 
#' newObj <- methylBase.obj
#' getSampleID(newObj) <- c("sample1","sample2","sample3","sample4")
#' getSampleID(newObj)
#' 
#' 
#' @export
#' @docType methods
#' @rdname getSampleID-methods
setGeneric("getSampleID", def=function(x) standardGeneric("getSampleID"))

#' @export
#' @rdname getSampleID-methods
setGeneric("getSampleID<-", def=function(x, value) {
  standardGeneric("getSampleID<-")})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylRawList-method
setMethod("getSampleID", signature = "methylRawList", function(x) {
  names <- vapply(x,function(z) z@sample.id,FUN.VALUE = "character")
  return(names)
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylRawList-method
setReplaceMethod("getSampleID", signature = "methylRawList", function(x, value) {
  
  if(! ( length(getSampleID(x)) == length(value) ) ){
    stop("The vector of new sample ids is not valid, check the length of input")
  } else {
    treatment <- x@treatment
    x <- mapply(`getSampleID<-`, x, value)
    x <- new("methylRawList",x,treatment=treatment)
    return(x)
  }
  
})


#' @rdname getSampleID-methods
#' @aliases getSampleID,methylBase-method
setMethod("getSampleID", signature = "methylBase", function(x) {
  return(x@sample.ids)
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylBase-method
setReplaceMethod("getSampleID", signature = "methylBase", function(x, value) {
  
  if(! ( length(x@sample.ids) == length(value) ) ){
    stop("The vector of new sample ids is not valid, check the length of input")
  } else {
    x@sample.ids <- value
    return(x)
  }
  
})


#' @rdname getSampleID-methods
#' @aliases getSampleID,methylRaw-method
setMethod("getSampleID", signature = "methylRaw", function(x) {
  return(x@sample.id)
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylRaw-method
setReplaceMethod("getSampleID", signature = "methylRaw", function(x, value) {
  
  if(! ( length(value) == 1 ) ){
    stop("The vector of new sample ids is not valid, check the length of input")
  } else {
    x@sample.id <- value
    return(x)
  }
  
})


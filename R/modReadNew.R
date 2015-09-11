require(data.table)
.structureAMPoutput<-function(data)
{  
  strand=rep("+",nrow(data))
  strand[data[,4]=="R"]="-"
  numCs=round(data[,5]*data[,6]/100)
  numTs=round(data[,5]*data[,7]/100)
  
  data.frame(chr=data[,2],start=data[,3],end=data[,3]
             ,strand=strand,coverage=data[,5],numCs=numCs,numTs=numTs)
}

# reformats a generic structure data.frame to a standard methylraw data.frame
# based on the column number assignment and if freqC is fraction or not.
.structureGeneric<-function(data, pipeline)
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
    numTs   =cpgFR$numTs.x + cpgFR$numTs.y ,stringsAsFactors =F
  )
  Fid=paste(cpgF$chr,cpgF$start,cpgF$end)
  Rid=paste(cpgR$chr,cpgR$start,cpgR$end)
  resid=paste(res$chr,res$start,res$end)  
  res=rbind(res, cpgF[ !  Fid  %in%  resid,],cpgR[ ! Rid  %in%  resid,] )
  #res=res[order(res$chr,res$start),]
  return(res)
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
    cat("data part of methylRaw have",ncol(data),"columns, expected 7 columns")
  }
  
}



#' read file(s) to a methylrawList or methylraw object
#'
#' The function reads a list of files or files with methylation information for bases/region in the genome and creates a methylrawList or methylraw object
#' @param location file location(s), either a list of locations (each a character string) or one location string
#' @param sample.id sample.id(s)
#' @param assembly a string that defines the genome assembly such as hg18, mm9
#' @param header if the input file has a header or not (default: TRUE)
#' @param skip number of lines to skip when reading. Can be set to 1 for bed files with track line (default: 0)
#' @param sep seperator between fields, same as \code{\link{read.table}} argument (default: "\t")
#' @param pipeline name of the alignment pipeline, it can be either "amp" or "bismark". The methylation text files generated from other pipelines can be read as generic methylation text files by supplying a named \code{\link[base]{list}} argument as "pipeline" argument.
#' The named \code{list} should containt column numbers which denotes which column of the text file corresponds to values and genomic location of the methylation events. See Details for more.
#' @param resolution designates whether methylation information is base-pair resolution or regional resolution. allowed values 'base' or 'region'. Default 'base'
#' @param treatment a vector contatining 0 and 1 denoting which samples are control which samples are test
#' @param context methylation context string, ex: CpG,CpH,CHH, etc. (default:CpG)
#' @param dbdir directory where flat file database(s) should be stored, defaults
#'       to getwd(), working directory.
#' @param dbtype type of the flat file database, currently only option is "tabix"
#'        defaults to NULL, in which case the objects are stored in memory.
#' @examples
#' 
#' # this is a list of example files, ships with the package
#' # for your own analysis you will just need to provide set of paths to files
#' #you will not need the "system.file(..."  part
#' file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
#'                 system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
#'                 system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#'                 system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
#'
#' # read the files to a methylRawList object: myobj
#' myobj=modRead( file.list,
#'             sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",treatment=c(1,1,0,0))
#'             
#' # read one file as methylRaw object
#' myobj=modRead( file.list[[1]],
#'             sample.id="test1",assembly="hg18")
#'             
#' # read a generic text file containing CpG methylation values
#' # let's first look at the content of the file
#' generic.file=system.file("extdata", "generic1.CpG.txt", package = "methylKit")
#' read.table(generic.file,header=TRUE)
#' 
#' # And this is how you can read that generic file as a methylKit object            
#'  myobj=modRead( generic.file,pipeline=list(fraction=FALSE, chr.col=1,start.col=2,end.col=2,coverage.col=4,strand.col=3,freqC.col=5),
#'             sample.id="test1",assembly="hg18")
#'             
#' @section Details:
#'  When \code{pipeline} argument is a list, it is exptected to provide a named list with following names.
#'  'fraction' is a logical value, denoting if the column frequency of Cs has a range from [0-1] or [0-100]. If true it assumes range is [0-1].
#'  'chr.col" is the number of the column that has chrosome string.   
#'  'start.col' is the number of the column that has start coordinate of the base/region of the methylation event.
#'  'end.col'  is the number of the column that has end coordinate of the base/region of the methylation event.
#'  'coverage.col' is the number of the column that has read coverage values. 
#'  'strand.col' is the number of the column that has strand information, the strand information in the file has to be in the form of '+' or '-', 
#'  'freqC.col' is the number of the column that has the frequency of Cs. See examples to see how to read a generic methylation text file.
#'  
#' @return returns methylRaw or methylRawList
#' 
#' @export
#' @docType methods
#' @rdname modRead-methods
setGeneric("modRead", function(location,sample.id,assembly,pipeline="amp",
                               header=T,skip=0,sep="\t",
                               context="CpG",resolution="base",
                               treatment,dbdir=getwd()) standardGeneric("modRead"))



#' @rdname modRead-methods
#' @aliases modRead,character,character,character-method
setMethod("modRead", signature(location = "character",sample.id="character",assembly="character"),
          
          function(location,sample.id,assembly,pipeline,header,skip,sep,context,resolution){ 
            if(! file.exists(location)){stop(location,", That file doesn't exist !!!")}
            data<- as.data.frame( data.table::fread(location,header=header,skip=skip,sep=sep)  )  
            if(length(pipeline)==1 ){
              
              if(pipeline %in% c("amp","bismark") )
              {
                data<- .structureAMPoutput(data)
              }
              else{stop("unknown 'pipeline' argument, supported alignment pipelines: 'amp' or 'bismark' " )
              }
              
            }
            else{
              .check.pipeline.list(pipeline)
              data<- .structureGeneric(data, pipeline)
            }
            
            obj=new("methylRaw",data,sample.id=sample.id,assembly=assembly,context=context,resolution="base")
            obj         
          }
)


# reads a list of CpG methylation files and makes methylRawList object
#
# @param a list containing locations(full paths) to CpG methylation files from alignment pipeline
# @param name a list of strings that defines the experiment
# @param assembly a string that defines the genome assembly such as hg18, mm9
# @param pipeline name of the alignment pipeline, currently only supports AMP (default: AMP), or for generic read, a list object contain \code{fraction}=TRUE/FALSE, \code{chr.col}, \code{strand.col}, \code{start.col}, \code{end.col}, \code{coverage.col},\code{freqC.col}, for example: \code{list(fraction=T, chr.col=1, strand.col=2, coverage.col=3, freqC.col=4, start.col=5, end.col=5)}  
# @param header if the input files has a header or not (default: TRUE)
# @param treatment a vector contatining 0 and 1 denoting which samples are control which samples are test
# @return returns a methylRawList object
#' @rdname modRead-methods
#' @aliases modRead,list,list,character-method
setMethod("modRead", signature(location = "list",sample.id="list",assembly="character"),
          function(location,sample.id,assembly,pipeline,header,skip,sep,context,resolution,treatment){ 
            
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
              #data<- .readTableFast(location[[i]],header=header,skip=skip,sep=sep)# read data
              data<- as.data.frame( data.table::fread(location[[i]],header=header,skip=skip,sep=sep)  )  
              
              if(length(pipeline)==1 )
              {
                if(pipeline %in% c("amp","bismark")){
                  data<- .structureAMPoutput(data)
                } else {
                  stop("pipeline length is equal to 1 and is not amp or bismark. If you do not have amp or bismark format, please give a parameter list containing the format information of the data. Please refer details in the read help page")
                }
              }
              else{
                #stop("unknown 'pipeline' argument, supported alignment pipelines: amp")
                .check.pipeline.list(pipeline)
                data<- .structureGeneric(data, pipeline)
              }
              
              obj=new("methylRaw",data,sample.id=sample.id[[i]],assembly=assembly,context=context,resolution=resolution)
              outList[[i]]=obj       
            }
            myobj=new("methylRawList",outList,treatment=treatment)
            
            myobj
          })

# reads a list of CpG methylation files and makes methylRawList object
#
# @param a list containing locations(full paths) to CpG methylation files from alignment pipeline
# @param name a list of strings that defines the experiment
# @param assembly a string that defines the genome assembly such as hg18, mm9
# @param pipeline name of the alignment pipeline, currently only supports AMP (default: AMP), or for generic read, a list object contain \code{fraction}=TRUE/FALSE, \code{chr.col}, \code{strand.col}, \code{start.col}, \code{end.col}, \code{coverage.col},\code{freqC.col}, for example: \code{list(fraction=T, chr.col=1, strand.col=2, coverage.col=3, freqC.col=4, start.col=5, end.col=5)}  
# @param header if the input files has a header or not (default: TRUE)
# @param treatment a vector contatining 0 and 1 denoting which samples are control which samples are test
# @return returns a methylRawList object
#' @rdname modRead-methods
#' @aliases modRead,character,character,character-method
setMethod("modRead", signature(location = "character",sample.id="character",assembly="character"),
          
          function(location,sample.id,assembly,pipeline,header,skip,sep,context,resolution){ 
            if(! file.exists(location)){stop(location,", That file doesn't exist !!!")}
            data<- as.data.frame( data.table::fread(location,header=header,skip=skip,sep=sep)  )  
            if(length(pipeline)==1 ){
              
              if(pipeline %in% c("amp","bismark") )
              {
                data<- .structureAMPoutput(data)
              }
              else{stop("unknown 'pipeline' argument, supported alignment pipelines: 'amp' or 'bismark' " )
              }
              
            }
            else{
              .check.pipeline.list(pipeline)
              data<- .structureGeneric(data, pipeline)
            }
            
            obj=new("methylRaw",data,sample.id=sample.id,assembly=assembly,context=context,resolution="base")
            obj         
          }
)

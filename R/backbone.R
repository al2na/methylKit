 
#---------------------------------------------------------------------------------------
# regular R functions to be used in S4 functions

# reads a table in a fast way to a dataframe
.readTableFast<-function(filename,header=T,skip=0,sep="")
{
  tab5rows <- read.table(filename, header = header,skip=skip,sep=sep, nrows = 100)
  classes  <- sapply(tab5rows, class)
  return( read.table(filename, header = header,skip=skip,sep=sep, colClasses = classes)  )
}

# reformats a data.frame to a standard methylraw data.frame
# no matter what the alignment pipeline
.structureAMPoutput<-function(data)
{
  strand=rep("+",nrow(data))
  strand[data[,4]=="R"]="-"
  numCs=round(data[,5]*data[,6]/100)
  numTs=round(data[,5]*data[,7]/100)
  
  data.frame(id=data[,1],chr=data[,2],start=data[,3],end=data[,3]
             ,strand=strand,coverage=data[,5],numCs=numCs,numTs=numTs)
}

# unfies forward and reverse strand CpGs on the forward strand if the if both are on the same CpG
# if that's the case their values are generally correlated
.CpG.dinuc.unify<-function(cpg)
{

  cpgR=cpg[cpg$strand=="-",]
  cpgF=cpg[cpg$strand=="+",]
  cpgR$start=cpgR$start-1
  cpgR$end=cpgR$end-1
  cpgR$strand="+"
  
  cpgR$id=paste(cpgR$chr,cpgR$start,sep=".")

  cpgFR=merge(cpgF,cpgR,by="id")
  #hemi =cpgFR[abs(cpgFR$freqC.x-cpgFR$freqC.y)>=50,]
  #cpgFR=cpgFR[abs(cpgFR$freqC.x-cpgFR$freqC.y)<50,]  
  res=data.frame(
                 id =as.character(cpgFR$id),
                 chr     =as.character(cpgFR$chr.x),
                 start    =(cpgFR$start.x),
                 end      =(cpgFR$start.x),
                 strand  =rep("+",nrow(cpgFR)),
                 coverage=cpgFR$coverage.x + cpgFR$coverage.y,
                 numCs   =cpgFR$numCs.x + cpgFR$numCs.y ,
                 numTs   =cpgFR$numTs.x + cpgFR$numTs.y ,stringsAsFactors =F
  )
  res=rbind(res, cpgF[ !cpgF$id  %in%  res$id,],cpgR[ !cpgR$id  %in%  res$id,] )
  res=res[order(res$id),]
  return(res)
}


# end of regular functions to be used in S4 functions
#---------------------------------------------------------------------------------------






#' An S4 class for holding raw methylation data from alignment pipeline.
#'
#' This object stores the raw mehylation data that is read in through read function and extends data.frame
#'
#' @section Slots:\describe{
#'                  \item{\code{sample.id}:}{string for an identifier of the sample}
#'                  \item{\code{assembly}:}{string for genome assembly, ex: hg18,hg19,mm9}
#'                  \item{\code{context}:}{ methylation context string, ex: CpG,CpH,CHH, etc.}
#'                 }
#' @name methylRaw-class
#' @rdname methylRaw-class
#' @export
setClass("methylRaw", representation(
  sample.id = "character", assembly = "character",context="character"),contains= "data.frame")


#' An S4 class for holding a list of methylRaw objects.
#'
#' This object stores the list of raw mehylation data that is read in through read function and extends data.frame
#'
#' @section Slots:\describe{
#'                  \item{\code{treatment}:}{numeric vector denoting control and test samples}
#'                }
#' @name methylRawList-class
#' @rdname methylRawList-class
#' @export
setClass("methylRawList", representation(treatment = "numeric"),contains = "list")

#' read file(s) to a methylrawList or methylraw object
#'
#' read a list of locations or one location and create a methylrawList or methylraw object
#' @param location file location(s), either a list of locations (each a character string) or one location string
#' @param sample.id sample.id(s)
#' @param assembly a string that defines the genome assembly such as hg18, mm9
#' @param header if the input file has a header or not (default: TRUE)
#' @param pipeline name of the alignment pipeline, currently only supports AMP (default: AMP)
#' @param treatment a vector contatining 0 and 1 denoting which samples are control which samples are test
#' @param context methylation context string, ex: CpG,CpH,CHH, etc. (default:CpG)
#' @usage read(location,sample.id,assembly,pipeline="amp",header=T, context="CpG",treatment)
#' @return returns methylRaw or methylRawList
#' @export
#' @docType methods
#' @rdname read-methods
setGeneric("read", function(location,sample.id,assembly,pipeline="amp",header=T, context="CpG",treatment) standardGeneric("read"))


# read one CpG methylation file that is a result of alignment pipeline
#
# @param location full path to CpG methylation file
# @param name a string that defines the experiment
# @param assembly a string that defines the genome assembly such as hg18, mm9
# @param pipeline name of the alignment pipeline, currently only supports AMP (default: AMP)
# @param header if the input file has a header or not (default: TRUE)
#' @rdname read-methods
#' @aliases read,character,character,character-method
setMethod("read", signature(location = "character",sample.id="character",assembly="character"),
          
        function(location,sample.id,assembly,pipeline,header,context){ 
            if(! file.exists(location)){stop(location,", That file doesn't exist !!!")}
            data<- .readTableFast(location,header=header)            
            if(pipeline=="amp")
            {
              data<- .structureAMPoutput(data)
            }
            else{
                stop("unknown 'pipeline' argument, supported alignment pipelines: amp ")
            }

            obj=new("methylRaw",data,sample.id=sample.id,assembly=assembly,context=context)
            obj         
          }
)


# reads a list of CpG methylation files and makes methylRawList object
#
# @param a list containing locations(full paths) to CpG methylation files from alignment pipeline
# @param name a list of strings that defines the experiment
# @param assembly a string that defines the genome assembly such as hg18, mm9
# @param pipeline name of the alignment pipeline, currently only supports AMP (default: AMP)
# @param header if the input files has a header or not (default: TRUE)
# @param treatment a vector contatining 0 and 1 denoting which samples are control which samples are test
# @return returns a methylRawList object
#' @rdname read-methods
#' @aliases read,list,list,character-method
setMethod("read", signature(location = "list",sample.id="list",assembly="character"),
          function(location,sample.id,assembly,pipeline,header,context,treatment){ 
            
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
              data<- .readTableFast(location[[i]],header=header)# read data
              if(pipeline=="amp")
              {
                data<- .structureAMPoutput(data)
              }
              else{
                stop("unknown 'pipeline' argument, supported alignment pipelines: amp")
              }
                  
              obj=new("methylRaw",data,sample.id=sample.id[[i]],assembly=assembly,context=context)
              outList[[i]]=obj       
            }
            myobj=new("methylRawList",outList,treatment=treatment)
            
            myobj
          })


#' filter methylRaw and methylRawList object based on read coverage
#'
#' This function filters \code{methylRaw} and \code{methylRawList} objects.
#' You can filter based on lower read cutoff or high read cutoff. Higher read cutoff is usefull to eliminate PCR effects
#' Lower read cutoff is usefull for doing better statistical tests.
#'
#' @param methylObj a \code{methylRaw} or \code{methylRawList} object
#' @param lo.count An integer for read counts.Bases/regions having lower coverage than this count is discarded
#' @param lo.perc  A double [0-100] for percentile of read counts. Bases/regions having lower coverage than this percentile is discarded
#' @param hi.count An integer for read counts. Bases/regions having higher coverage than this is count discarded
#' @param hi.perc A double [0-100] for percentile of read counts. Bases/regions having higher coverage than this percentile is discarded
#' @usage filterByCoverage(methylObj,lo.count=NULL,lo.perc=NULL,hi.count=NULL,hi.perc=NULL)
#' @return \code{methylRaw} or \code{methylRawList} object depending on input object
#' @export
#' @docType methods
#' @rdname filterByCoverage-methods
setGeneric("filterByCoverage",function(methylObj,lo.count=NULL,lo.perc=NULL,hi.count=NULL,hi.perc=NULL) standardGeneric("filterByCoverage") )

#' @aliases filterByCoverage,methylRaw-method
#' @rdname filterByCoverage-methods
setMethod("filterByCoverage", signature(methylObj="methylRaw"),
                    function(methylObj,lo.count,lo.perc,hi.count,hi.perc){
                      if( is.null(lo.count) & is.null(lo.perc) & is.null(hi.count) & is.null(hi.perc) ){return(methylObj)}
                      
                      data=getData(methylObj) # get the data part
                      
                      #figure out which cut-offs to use, maybe there is more elagent ways, quick&dirty works for now
                      if(is.numeric(lo.count) ){lo.count=lo.count}
                      if(is.numeric(lo.perc)){lo.count=quantile(data$coverage,lo.perc/100)}
                      if(is.numeric(hi.count)){hi.count=hi.count}
                      if(is.numeric(hi.perc)){hi.count=quantile(data$coverage,hi.perc/100)}
                      
                      if(is.numeric(lo.count)){data=data[data$coverage>=lo.count,]}
                      if(is.numeric(hi.count)){data=data[data$coverage<hi.count,]}
                      

                      new("methylRaw",data,sample.id=methylObj@sample.id,
                                           assembly=methylObj@assembly,
                                           context=methylObj@context)

                      
})

#' @aliases filterByCoverage,methylRawList-method
#' @rdname filterByCoverage-methods
setMethod("filterByCoverage", signature(methylObj="methylRawList"),
                    function(methylObj,lo.count,lo.perc,hi.count,hi.perc){
                      new.list=lapply(methylObj,filterByCoverage,lo.count,lo.perc,hi.count,hi.perc)
                      new("methylRawList", new.list,treatment=methylObj@treatment)
})


#' An S4 class that holds base-pair resolution methylation information for multiple experiments, only bases that are covered in all experiments are held in this class
#'
#' extends data.frame and creates an object that holds methylation information and genomic location
#'          
#' @section Slots:\describe{
#'                  \item{sample.ids}{character vector for ids of samples in the object}
#'
#'                  \item{assembly}{name of the genome assembly}
#'
#'                  \item{context}{context of methylation. Ex: CpG,CpH,CHH, etc}
#'
#'                  \item{treatment}{treatment vector denoting which samples are test and control}
#'
#'                  \item{coverage.index}{vector denoting which columns in the data correspons to coverage values}
#'
#'                  \item{numCs.index}{vector denoting which columns in the data correspons to number of methylatedCs values}
#'
#'                  \item{numTs.index}{vector denoting which columns in the data correspons to number of unmethylated Cs values}
#' }
#' @name methylBase-class
#' @rdname methylBase-class
#' @export
setClass("methylBase",contains="data.frame",representation(
  sample.ids = "character", assembly = "character",context = "character",treatment="numeric",coverage.index="numeric",
                                   numCs.index="numeric",numTs.index="numeric",destranded="logical"))


#' unites methylRawList to a single table 
#' 
#' This functions unites \code{methylRawList} object that only bases with coverage from all samples are retained.
#' The resulting object is a class of \code{methylBase}
#'
#' @param .Object a methylRawList object to be merged by common locations covered by reads
#' @param destrand if TRUE, reads covering both strands of a CpG dinucleotide will be merged, 
#'   do not set to TRUE if not only interested in CpGs (default: FALSE)
#'
#' @usage unite(.Object,destrand=F)
#' @return a methylBase object
#' @aliases unite,-methods unite,methylRawList-method
#' @export
#' @docType methods
#' @rdname unite-methods
setGeneric("unite", function(.Object,destrand=F) standardGeneric("unite"))

#' @rdname unite-methods
#' @aliases unite,methylRawList-method
setMethod("unite", "methylRawList",
                    function(.Object,destrand){
  
                    
                     #merge raw methylation calls together
                     df=getData(.Object[[1]])
                     if(destrand){df=.CpG.dinuc.unify(df)}
                     sample.ids=c(.Object[[1]]@sample.id)
                     assemblies=c(.Object[[1]]@assembly)
                     contexts  =c(.Object[[1]]@context)
                     for(i in 2:length(.Object))
                     {
                       df2=getData(.Object[[i]])
                       if(destrand){df2=.CpG.dinuc.unify(df2)}
                       df=merge(df,df2[,c(1,6:8)],by="id") # merge the dat to a data.frame
                       sample.ids=c(sample.ids,.Object[[i]]@sample.id)
                       contexts=c(contexts,.Object[[i]]@context)
                     }

                     # stop if the assembly of object don't match
                     if( length( unique(assemblies) ) != 1 ){stop("assemblies of methylrawList elements should be same\n")}
          
                     # get indices of coverage,numCs and numTs in the data frame 
                     coverage.ind=seq(6,by=3,length.out=length(.Object))
                     numCs.ind   =seq(6,by=3,length.out=length(.Object))+1
                     numTs.ind   =seq(6,by=3,length.out=length(.Object))+2

                     # change column names
                     names(df)[coverage.ind]=paste(c("coverage"),1:length(.Object),sep="" )
                     names(df)[numCs.ind]   =paste(c("numCs"),1:length(.Object),sep="" )
                     names(df)[numTs.ind]   =paste(c("numTs"),1:length(.Object),sep="" )

                     #make methylbase object and return the object
                     obj=new("methylBase",as.data.frame(df),sample.ids=sample.ids,
                             assembly=unique(assemblies),context=unique(contexts),
                             treatment=.Object@treatment,coverage.index=coverage.ind,
                             numCs.index=numCs.ind,numTs.index=numTs.ind,destranded=destrand )
                     obj
                    }
          )
            

#' get correlation between samples in methylBase object
#' 
#' @param .Object a methylBase object 
#' @param plot scatterPlot if TRUE (deafult:False) 
#' @return a correlation matrix object and plot scatterPlot
#' @aliases getCorrelation,-methods getCorrelation,methylBase-method
#' @export
#' @docType methods
#' @rdname getCorrelation-methods
setGeneric("getCorrelation", function(.Object,plot=F) standardGeneric("getCorrelation"))

#' @rdname getCorrelation-methods
#' @aliases getCorrelation-method
setMethod("getCorrelation", "methylBase",
                    function(.Object,plot){
                        meth.mat = getData(.Object)[, .Object@numCs.index]/(.Object[,.Object@numCs.index] + .Object[,.Object@numTs.index] )                                      
                        names(meth.mat)=.Object@sample.ids
                        
                        print( cor(meth.mat) )
                        
                        panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
                        {
                          usr <- par("usr"); on.exit(par(usr))
                          par(usr = c(0, 1, 0, 1))
                          r <- abs(cor(x, y))
                          txt <- format(c(r, 0.123456789), digits=digits)[1]
                          txt <- paste(prefix, txt, sep="")
                          if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
                          text(0.5, 0.5, txt, cex = cex.cor * r)
                        }
                          
                        panel.my.smooth2<-function(x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = "darkgreen", span = 2/3, iter = 3, ...) 
                        {
                             par(new = TRUE)    #par(usr = c(usr[1:2], 0, 1.5) )
                            smoothScatter(x, y,colramp=colorRampPalette(topo.colors(100)), bg = bg)
                            ok <- is.finite(x) & is.finite(y)
                            if (any(ok)) 
                                lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = col.smooth, ...)
                                abline(lm(y[ok]~x[ok]), col="red")
                        }
                        
                        panel.my.smooth<-function(x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 0.3, col.smooth = "green", span = 2/3, iter = 3, ...) 
                        {
                            points(x, y, pch = 20, col = densCols(x,y,colramp=colorRampPalette(topo.colors(20))), bg = bg, cex = 0.1)
                            ok <- is.finite(x) & is.finite(y)
                            if (any(ok)){
                                lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = col.smooth, ...);
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
                          pairs(meth.mat, 
                              lower.panel=panel.my.smooth2, 
                              upper.panel=panel.cor,
                              diag.panel=panel.hist,main="CpG dinucleotide correlation")
                        }
                    }  
 )



#' get coverage stats from methylRaw object
#' 
#' @param .Object a \code{methylRaw} object 
#' @param plot plot a histogram of coverage if TRUE (deafult:FALSE) 
#' @param both.strands do stats and plot for both strands if TRUE (deafult:FALSE)
#' @return a summary of coverage statistics or plot a histogram of coverage
#' @aliases getCoverageStats,-methods getCoverageStats,methylRaw-method
#' @export
#' @docType methods
#' @rdname getCoverageStats-methods
setGeneric("getCoverageStats", function(.Object,plot=F,both.strands=F) standardGeneric("getCoverageStats"))

#' @rdname getCoverageStats-methods
#' @aliases getCoverageStats,methylRaw-method
setMethod("getCoverageStats", "methylRaw",
                    function(.Object,plot,both.strands){
                      
                      if(!plot){
                        qts=seq(0,0.9,0.1) # get quantiles
                        qts=c(qts,0.95,0.99,0.995,0.999,1)                          
                        
                        if(both.strands){       
                          plus.cov=.Object[.Object$strand=="+",]$coverage
                          mnus.cov=.Object[.Object$strand=="-",]$coverage
                          
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
                          
                          all.cov=.Object$coverage
                          
                          cat("read coverage statistics per base\n")
                          cat("summary:\n")
                          print( summary( all.cov ) )
                          cat("percentiles:\n")
                          print(quantile( all.cov,p=qts ))
                          cat("\n")
                        }
                        
                      }else{
                        if(both.strands){   
                          plus.cov=.Object[.Object$strand=="+",]$coverage
                          mnus.cov=.Object[.Object$strand=="-",]$coverage
                          
                          par(mfrow=c(1,2))
                          a=hist(log10(plus.cov),plot=F)
                          hist(log10(plus.cov),col="chartreuse4",
                               xlab="log10 of read coverage per base",
                               main="Histogram of CpG coverage: Forward strand",
                               labels=as.character(round(100*a$counts/length(plus.cov),1)))
                          
                          a=hist(log10(mnus.cov),plot=F)
                          hist(log10(mnus.cov),col="chartreuse4",
                               xlab="log10 of read coverage per base",
                               main="Histogram of CpG coverage: Reverse strand",
                               labels=as.character(round(100*a$counts/length(mnus.cov),1)))
 
                        }else{
                          all.cov= .Object$coverage
                          
                          a=hist(log10(all.cov),plot=F)
                          hist(log10(all.cov),col="chartreuse4",
                               xlab="log10 of read coverage per base",
                               main="Histogram of CpG coverage",
                               labels=as.character(round(100*a$counts/length(all.cov),1)))

                        }
                        
                        
                      }
                        
                      
})

#' get Methylation stats from methylRaw object
#' 
#' @param .Object a \code{methylRaw} object 
#' @param plot plot a histogram of Methylation if TRUE (deafult:FALSE) 
#' @param both.strands do plots and stats for both strands seperately  if TRUE (deafult:FALSE)
#' @return a summary of Methylation statistics or plot a histogram of coverage
#' @export
#' @docType methods
#' @rdname getMethylationStats-methods
setGeneric("getMethylationStats", function(.Object,plot=F,both.strands=F) standardGeneric("getMethylationStats"))

#' @rdname getMethylationStats-methods
#' @aliases getMethylationStats,methylRaw-method
setMethod("getMethylationStats", "methylRaw",
                    function(.Object,plot,both.strands){
                      
                      plus.met=100* .Object[.Object$strand=="+",]$numCs/.Object[.Object$strand=="+",]$coverage
                      mnus.met=100* .Object[.Object$strand=="-",]$numCs/.Object[.Object$strand=="-",]$coverage
                      all.met =100* .Object$numCs/.Object$coverage
                      
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
                          a=hist((plus.met),plot=F)
                          hist((plus.met),col="cornflowerblue",
                               xlab="% methylation per CpG",
                               main="Histogram of % methylation: Forward strand",
                               labels=as.character(round(100*a$counts/length(plus.met),1)))
                          
                          a=hist((mnus.met),plot=F)
                          hist((mnus.met),col="cornflowerblue",
                               xlab="% methylation per base",
                               main="Histogram of % methylation: Reverse strand",
                               labels=as.character(round(100*a$counts/length(mnus.met),1)))
 
                        }else{
                          
                          a=hist((all.met),plot=F)
                          hist((all.met),col="cornflowerblue",
                               xlab="% methylation per base",
                               main="Histogram of % methylation",
                               labels=as.character(round(100*a$counts/length(all.met),1)))

                        }
                        
                        
                      }
                        
                      
})





# get distribution of difference between samples in methylBase object
# unites methylrawlist objects based on chromosomal positions of CpG dinucleotides
#setGeneric("getDifference", function(.Object,plot=F) standardGeneric("getCorrelation"))
#setMethod("getDifference", "methylBase",
#                    function(.Object,plot){
#                        meth.mat = .Object@data[, .Object@numCs.index]/(.Object@data[,.Object@numCs.index] + .Object@data[,.Object@numTs.index] )                                      
#                        names(meth.mat)=.Object@sample.ids                      
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
#  methylBase accessor functions
#


#' get assembly of the genome
#' 
#' @param x a methylBase object 
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
#' @param x a methylBase/methylRaw/methylDiff object 
#' @return the context of methylation string
#' @export
#' @docType methods
#' @rdname getContext-methods
setGeneric("getContext", def=function(x) standardGeneric("getContext"))

#' @rdname getContext-methods
#' @aliases getContext,methylBase-method
setMethod("getContext", signature="methylBase", definition=function(x) {
                return(x@Context)
        })

#' @rdname getContext-methods
#' @aliases getContext,methylRaw-method
setMethod("getContext", signature="methylRaw", definition=function(x) {
                return(x@Context)
        })



#' gets the data slot from the methylBase object
#' 
#' The data retrived from this function is of a \code{data.frame}. This is basically containing all relevant methylation information per region
#'
#' @param x a methylBase object 
#' @return data.frame for methylation events
#' @aliases getData,-methods getData,methylBase-method
#' @export
#' @docType methods
#' @rdname getData-methods
setGeneric("getData", def=function(x) standardGeneric("getData"))

#' @rdname getData-methods
#' @aliases getData-method
setMethod("getData", signature="methylBase", definition=function(x) {
                return(as(x,"data.frame"))
}) 

#' @rdname getData-methods
#' @aliases getData,methylRaw-method
setMethod("getData", signature="methylRaw", definition=function(x) {
                return(as(x,"data.frame"))
})

## CONVERTOR FUNCTIONS FOR methylRaw and methylBase OBJECT
#convert methylRaw to GRanges
setAs("methylRaw", "GRanges", function(from)
                      {
                        from2=getData(from)
                        GRanges(seqnames=from2$chr,ranges=IRanges(start=from2$start, end=from2$end),
                                       strand=from2$strand, 
                                       id=from2$id,
                                       coverage=from2$coverage,
                                       numCs   =from2$numCs,
                                       numTs  =from2$numTs                                
                                       )

})

setAs("methylBase", "GRanges", function(from)
                      {
                        #from=getData(from1)
                        GRanges(seqnames=from$chr,ranges=IRanges(start=from$start, end=from$end),
                                       strand=from$strand, 
                                       id=from$id,
                                       data.frame(from[,6:ncol(from)])
                                       )

})


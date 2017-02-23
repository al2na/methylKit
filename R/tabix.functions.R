# tabix functions
# these functions can be used to manipulate tabix files and methylKit objects
# that depend on those files


# write to file without using scientific notation
.write.table.noSci <- function(...,scipen=999) {
  
  defsci=options()$scipen
  options(scipen = 999)
  
  write.table(...)

  options(scipen = defsci)
}

#' merge tabix files by chr, start,end, strand
#' 
#' @param tabixList list of tabix files
#' @param dir working directory
#' @param filename the output file name
#' @param mc.cores number of multiple cores. If mc.cores>1 temporary files for
#'   each chromosome will be created prior to cat, zipping and indexing the
#'   single output file
#' @param all logical parameter passed to \code{\link{merge}} function
#'   
#'   mergeTabix(tabixList,dir="~",filename="dump.meth.txt",mc.cores=1)
#'   
#' @usage mergeTabix(tabixList,dir,filename,mc.cores=1 ,all=FALSE)
#' @noRd
mergeTabix<-function(tabixList,dir,filename,mc.cores=1 ,all=FALSE){
  
  # get outfile
  
  #filename="dump.meth.txt"
  #dir="~" 
  if( class(tabixList) != "TabixFileList" ){
    tabixList <- Rsamtools::TabixFileList(tabixList)
  }
  
  # get chrs
  chrNum=table(unlist(lapply(tabixList,Rsamtools::seqnamesTabix)))
  chrs=names(chrNum)[chrNum==length(tabixList)]
  
  if(mc.cores > 1){
    # random file string
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(rndFile,filename,sep="_")
    outlist=mclapply(chrs,mergeTbxByChr,tabixList,dir,filename2,parallel=TRUE,
                     mc.cores=mc.cores,all=all)
    
    #concat subfiles

    outfile= file.path(path.expand(dir),filename) # get file name 
    if(file.exists(outfile)){
      message("overwriting ",outfile)
      unlink(outfile)
    }
  
    con=file(outfile, open = "a", blocking = TRUE) # open connection  
    for(file in list.files(path = dir, pattern = rndFile,full.names=TRUE)){
      file.append(outfile,file) # append files
    }
    close(con)
    #remove temp files
    unlink(list.files(path = dir, pattern = rndFile,full.names=TRUE))
  }else{
    
    outfile= file.path(path.expand(dir),filename) # get file name 
    if(file.exists(outfile)){
      message("overwriting ",outfile)
      unlink(outfile)
    }
    
    outlist=lapply(chrs,mergeTbxByChr,tabixList,dir,filename,
                     parallel=FALSE,all=all)
    outfile= file.path(path.expand(dir),filename) 
  }
  
  #df=fread(outfile,header=FALSE)
  
  # after iterations end
  # zip and index the text file
  makeMethTabix( outfile ,skip=0)
}


#' make a tabix file from a data.frame or data.table part of a 
#' methylKit object
#' 
#' @param df data part of a methylKit object, either data.frame or data.table
#' @param outfile name of the output file
#' 
#' @usage df2tabix(df,outfile)
#' @noRd
df2tabix<-function(df,outfile){
  
  if(file.exists(outfile)){
    message("overwriting ",outfile)
    unlink(outfile)
  }
  
  # write the file to disk
  .write.table.noSci(df,outfile,quote=FALSE,sep="\t",
              col.names=FALSE,row.names=FALSE)
  
  
  #make tabix out if the file
  makeMethTabix( outfile ,skip=0)
}


#' cat subfiles and make tabix
#' @param dir directory of subfiles to be merged
#' @param pattern a pattern in the file names
#' @param filename output file name
#' @param sort sort list of subfiles in alpha, numerical way
#' 
#' @usage catsub2tabix(dir,pattern,filename,sort=F)
#' @noRd
catsub2tabix<-function(dir,pattern,filename,sort=FALSE){
  
  outfile= file.path(path.expand(dir),filename) # get file name 
  if(file.exists(outfile)){
    message("overwriting ",outfile)
    unlink(outfile)
  }
  subfiles <- list.files(path = dir, pattern = pattern,full.names=TRUE)
  if(sort) {subfiles <- gtools::mixedsort(subfiles)}
  
  con=file(outfile, open = "a", blocking = TRUE) # open connection  
  for(file in subfiles){
    file.append(outfile,file) # append files
  }
  close(con)
  #remove temp files
  unlink(subfiles)
  
  #make tabix
  makeMethTabix( outfile ,skip=0)
  
}

#' make tabix from a flat file where first 3 columns
#' are chr,start,end in that order
#' 
#' @param filepath path to the uncompressed file
#' @param skip number of lines to skip
#' @param rm.file remove the uncompressed text file (default: yes)
#' 
#' @usage makeMethTabix(filepath,skip=0)
#' @noRd
makeMethTabix<-function(filepath,skip=0,rm.file=TRUE){
  message("compressing the file with bgzip...")
  zipped <- Rsamtools::bgzip(filepath,overwrite=TRUE)
  
  if(rm.file){file.remove(filepath)}
  
  message("making tabix index...")
  Rsamtools::indexTabix(zipped,
             seq=1, start=2, end=3,
             skip=skip, comment="#", zeroBased=FALSE)
  
}


# create short abbreviations for methylDB-slots
# used to save slots in tabix file
.encodeShortSlotNames <- function(slot.name) {
  switch(slot.name,
         sample.id="SI",
         num.records="NR", sample.ids="SI",
         assembly="AS", context="CT",
         resolution="RS", dbtype = "DT",
         treatment = "TM", coverage.index = "CI",
         numCs.index = "NC", numTs.index = "NT",
         destranded = "DS",
         NULL)
}

# decode the short abbreviations to long ones
# used to read slots from tabix file
.decodeShortSlotNames <- function(slot.name) {
  switch(slot.name,
         NR="num.records", SI="sample.ids",
         AS="assembly", CT="context",
         RS="resolution",  DT= "dbtype",
         TM = "treatment",  CI= "coverage.index",
         NC = "numCs.index",  NT= "numTs.index",
         DS = "destranded",
         NULL)
}

# format the content of slots read from tabix file
.formatShortSlotValues<- function(headerItem) {
  type <- headerItem[1]
  x <- headerItem[2]
  val <- switch(type,
                NR = as.integer(x),
                SI = unlist(strsplit(as.character(x),split = ",")),
                AS = as.character(x),
                CT = as.character(x),
                RS = as.character(x),
                DT = as.character(x),
                TM = as.integer(unlist(strsplit(as.character(x),split = ","))),
                CI = as.integer(unlist(strsplit(as.character(x),split = ","))),
                NC = as.integer(unlist(strsplit(as.character(x),split = ","))),
                NT = as.integer(unlist(strsplit(as.character(x),split = ","))),
                DS = as.logical(x),
                NULL
  )
  val <- list( val); names(val) <- type
  
  return( val)
}


#' make a tabix file from a methylKit object 
#'  
#' First writes the slots of the object into the header of the tabix file.
#' The slots are encoded like this:
#' "sample.id"/"sample.ids" --> "SI" ; "num.records" -->"NR";
#' "assembly" --> "AS"; "context" -->"CT"; 
#' "resolution" --> "RS"; "dbtype" --> "DT"; 
#' "treatment" --> "TM"; "coverage.index" --> "CI"; 
#' "numCs.index" --> "NC"; "numTs.index" --> "NT"; 
#' "destranded" --> "DS"
#' 
#' Then the data part of the object is appended to the tabix file.
#'
#' @param obj any methylKit object
#' @param filename name of the file where object is written to 
#' @param rm.txt should the uncompressed text file be removed (default: TRUE)
#' 
#' @usage 
#' @noRd
obj2tabix <- function(obj,filename,rm.txt=TRUE){
  
  # first we query each slots and ... 
  tabixHead <- sapply(slotNames(obj),
                      FUN = function(i) {
                        # convert it into a shorter form
                        ssn <- .encodeShortSlotNames(i)
                        # and then paste its values into a comma seperated string
                        if(!is.null(ssn)) paste(ssn,sep = ":",
                                                paste0(slot(obj,i),collapse = ",")
                        )
                      }
  )
  # then we remove slots we don't need in the header (as the @.Data slot) 
  tabixHead <- tabixHead[!isEmpty(tabixHead)]
  
  # then we write the creation date ..
  write(paste0("#Date:",format(Sys.time(),'%Y-%m-%d %H:%M:%S')),
        file = filename)
  # and the slots as comments 
  write(paste0("#",tabixHead),
        file = filename ,append = TRUE)
  
  
  if ("dbtype" %in% slotNames(obj)) {
    df <- headTabix(obj@dbpath,nrow = obj@num.records)
  } 
  else {
    df <- data.table::as.data.table(obj@.Data)
    write("#DT:tabix",file = filename ,append = TRUE)
    write(paste0("#NR:",nrow(df)),file = filename ,append = TRUE)
  }
  # sort the data
  df <- df[with(df,order(V1,V2,V3)),]
  # then we write the data
  write.table(x = df,
              file = filename, quote = FALSE,
              append = TRUE,col.names = FALSE,
              row.names = FALSE,sep = "\t")
  
  # and make tabix out of file
  makeMethTabix(filename,rm.file = rm.txt)
}


#' function to parse the information stored in the tabix header
#'
#' @param tbxFile tabix file
#' @noRd
readTabixHeader <- function(tbxFile) {
  
  # save header values
  h <- Rsamtools::headerTabix(tbxFile)$header
  
  # check if header is empty
  if(length(h) == 0 ) stop(paste("The Tabix File:",tbxFile,"does not include a header.\n"))
  
  # parse the slots and format the values
  headerVals <- sapply(h,
                       USE.NAMES = FALSE,
                       FUN = function(j) {
                         tmp <- unlist(strsplit(gsub("#", "", j),split = ":"))
                         val <- .formatShortSlotValues(tmp)
                         names(val) <- .decodeShortSlotNames(names(val))
                         return(val)
                       }
  )
  
  headerVals <- headerVals[!isEmpty(headerVals)]
  
  return(headerVals)
  
}


#' merge the data tables for a given chr
#' 
#' @param chr chromosome to merge 
#' @param tabixList list of tabix files
#' @param dir working directory
#' @param filename the output file name
#' @param parallel logical to determine internal processing of output
#' @param all logical parameter passed to \code{\link{merge}} function
#' @noRd
mergeTbxByChr<-function(chr,tabixList,dir,filename,parallel=FALSE,all=FALSE){
  
  #get first file on the list
  res=getTabixByChr(tbxFile = tabixList[[1]],chr = chr,return.type = "data.table")
  colnames(res)[5:7] <- paste(colnames(res)[5:7],"1",sep = ".")
  
  for(i in 2:length(tabixList)){
    
    # get tabix per Chr
    tmp=getTabixByChr(tbxFile = tabixList[[i]],chr = chr)
    colnames(tmp)[5:7] <- paste(colnames(tmp)[5:7],i,sep = ".")
    
    # merge tabix in memory
    res=merge(res,tmp,by=c("V1","V2","V3","V4"),all=all)
    
  }
  
  # order rows
  V1=V2=V3=NULL
  setorder(res, V1,V2,V3)
  
  if(!parallel){
  # write out append TRUE
  # if not parallel we can write out one file and append to it
  outfile= file.path(path.expand(dir),filename)   
  con=file(outfile, open = "a", blocking = TRUE)
  .write.table.noSci(res,con,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,
              append=TRUE)
  close(con)
  }else{
    
    # when parallel, we have to write seperate files for each chr
    # then merge them outside this function
    outfile= file.path(path.expand( dir),paste(chr,filename,sep="_")) 
    .write.table.noSci(res,outfile,quote=FALSE,sep="\t",
                col.names=FALSE,row.names=FALSE,append=TRUE)
    
  }
  return(outfile)
}



#' get data from meth tabix for a given chr
#'
# @example
# tbxFile=system.file("extdata", "ctrl1.txt", package = "methylKit")
#  getTabixByChr(chr="chr21",tbxFile)
#' @noRd
getTabixByChr<-function(tbxFile,chr="chr10",
                        return.type=c("data.table","data.frame","GRanges")){
  
  return.type <- match.arg(return.type)
  
  if( class(tbxFile) != "TabixFile" ){
    tbxFile <- TabixFile(tbxFile)
  }
  
  #res=Rsamtools:::.tabix_scan(tbxFile,space=chr,start=1L,end=500000000L) 
  param=GRanges(chr, IRanges(start=1L,end=500000000L))
  if(return.type=="data.table")
  {
    tabix2dt(scanTabix(tbxFile,param=param) )
  }else if (return.type=="data.frame"){
    tabix2df(scanTabix(tbxFile,param=param))
  }else {
    tabix2gr(scanTabix(tbxFile,param=param))
  }
}


#' get data from meth tabix for a given set of regions
#'
# @example
# granges <- as(methylRawListDB[[1]],"GRanges")
# tbxFile=methylRawListDB[[1]]@dbpath
#  getTabixByOverlap(granges=granges[1:50],tbxFile)
#' @noRd
getTabixByOverlap<-function(tbxFile,granges,return.type="data.table"){
  
  if( class(tbxFile) != "TabixFile" ){
    tbxFile <- TabixFile(tbxFile)
  }
  
  res=scanTabix(tbxFile,param=granges) 
  
  res <- list(unlist(res))
  
  if(return.type=="data.table")
  {
    tabix2dt( res )
  }else if (return.type=="data.frame"){
    tabix2df( res )
  }else {
    tabix2gr( res )
  }
  
}


#' get data from meth tabix for a given number of rows
#'
#' @example
#' tbxFile=methylRawListDB[[1]]@dbpath
#' headTabix(tbxFile)
#' @noRd
headTabix <- function(tbxFile, nrow = 10,
                    return.type = c("data.table","data.frame","GRanges") ){
  
  if(nrow < 1e6) {
    if( class(tbxFile) != "TabixFile" ){
      tbxFile <- TabixFile(tbxFile)
      open(tbxFile)
    }
    df <- getTabixByChunk( tbxFile,chunk.size=nrow,return.type)
  } 
  else {
    returnDt = if(return.type[1] == "data.table") TRUE else FALSE 
    df <- fread.gzipped(tbxFile,nrow = nrow, stringsAsFactors = TRUE, data.table = returnDt)
    
    if(return.type[1] == "GRanges"){
      return( GRanges(seqnames=as.character(df$V1),
              ranges=IRanges(start=df$V2, end=df$V3),
              strand=df$V4, df[,5:ncol(df)]) )
    } 
  }
  
  return(df)
  
}


#' get data from already opened tabixfile for a given chunkSize
#'
# @example
# tbxFile=methylRawListDB[[1]]@dbpath
#  getTabixByChunk( tbxFile,chunk.size=10)
#' @noRd
getTabixByChunk<-function(tbxFile,chunk.size=1e6,
                          return.type=c("data.table","data.frame","GRanges")){
  
  return.type <- match.arg(return.type)
  
  if( class(tbxFile) != "TabixFile" | !Rsamtools::isOpen(tbxFile, rw="read") ){
    stop("tbxFile has to be a class of TabixFile and should be open for reading ")
  }
  
  if(is.na(yieldSize(tbxFile)) | is.numeric(chunk.size)  ){
    yieldSize(tbxFile)<-chunk.size
  }
  
  res <- scanTabix(tbxFile)
  if(length(res) == 1 & length(res[[1]]) == 0)
    stop("the tabix file seems to be empty. stopping here.")
  
  if(return.type=="data.table")
  {
    tabix2dt(res)
  }else if (return.type=="data.frame"){
    tabix2df(res)
  }else {
    tabix2gr(res)
  }
}


#' convert methylKit tabix to data.table
# assuming you get a list length 1
#' @noRd
tabix2dt<-function(tabixRes){

    fread( paste0(paste(tabixRes[[1]],collapse="\n"),"\n" ),
                       stringsAsFactors=TRUE)
  
}

#' convert methylKit tabix to data.frame
# assuming you get a list length 1
#' @noRd
tabix2df<-function(tabixRes){

    fread( paste0(paste(tabixRes[[1]],collapse="\n"),"\n" ),
                       stringsAsFactors=TRUE,data.table = FALSE)
    
}

#' convert methylKit tabix to GRanges without Metacolumns
# assuming you get a list length 1
# for GRanges object with metacolumns coerce from methylDB object
#' @noRd
tabix2gr<-function(tabixRes){
  
    from <- fread(paste0(paste(tabixRes[[1]],collapse="\n"),"\n" ),
                              stringsAsFactors=TRUE, data.table = FALSE)
    
  GRanges(seqnames=as.character(from$V1),
          ranges=IRanges(start=from$V2, end=from$V3),
          strand=from$V4, from[,5:ncol(from)])

}

# applyTbxByChunk
#' Serially apply a function on chunks of tabix files
#' 
#' The function reads chunks of a tabix file and applies a function on them. 
#' The function (FUN argument) should apply on data.frames and 
#' return a data frame
#' as a result. The function is serially applied to chunks 
#' (means no parallelization). 
#' However, the function FUN itself can be a parallelized function
#' and related arguments could be passed on to the function via ... argument.
#' 
#' @param tbxFile tabix file to read. a TabixFile object
#' @param chunk.size number of rows to be taken as a chunk, default: 1e6
#' @param return.type indicates the return type for the function
#' @param FUN function to apply to chunks, it takes a data.frame and returns a 
#'            data.frame. First argument of the function should be a data frame
#' @param ... parameters to be passed to FUN. 
#' @param dir directory to create temporary files and the resulting tabix file
#' @param filename the filename for the resulting tabix file, this should not be 
#' a path, just a file name.
#'
#' @return either a path to a tabix or text file, or a data frame or data.table
#' @noRd
applyTbxByChunk<-function(tbxFile,chunk.size=1e6,dir,filename,
                          return.type=c("tabix","data.frame","data.table","text"),
                          FUN,...){
  
  return.type <- match.arg(return.type)
  FUN <- match.fun(FUN)
  
  # open tabix file with given chunk size
  if( class(tbxFile) != "TabixFile" ){
    tbxFile <- Rsamtools::TabixFile(tbxFile, yieldSize = chunk.size)

  } else {
    if(Rsamtools::isOpen(tbxFile)){close(tbxFile)}# close if already open
    Rsamtools::yieldSize(tbxFile) <-  chunk.size 
  }
  
  
  # calculate number of chunks
  recs=Rsamtools::countTabix(tbxFile)[[1]]
  chunk.num=ceiling(recs/chunk.size)
  
  open(tbxFile)
  
  if(return.type =="tabix"){
    
    # create a custom function that contains the function
    # to be applied
    myFunc<-function(chunk.num,tbxFile,dir,filename,FUN,...){
      data=getTabixByChunk(tbxFile,chunk.size=NULL,return.type="data.frame")
      res=FUN(data,...)  
      
      # for tabix
      outfile= file.path(path.expand( dir),paste(chunk.num,filename,sep="_"))
      .write.table.noSci(res,outfile,quote=FALSE,col.names=FALSE,row.names=FALSE,
                  sep="\t")
    }

    # attach a random string to the file name 
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(rndFile,filename,sep="_")
    
    # apply function to chunks
    res=lapply(1:chunk.num,myFunc,tbxFile,dir,filename2,FUN,...)
    
    # collect & cat temp files,then make tabix
    path <- catsub2tabix(dir,pattern=filename2,filename,sort = TRUE)

    return(gsub(".tbi","",path))
    
  } else if(return.type =="text"){
    
    # create a custom function that contains the function
    # to be applied
    myFunc2<-function(chunk.num,tbxFile,dir,filename,FUN,...){
      data=getTabixByChunk(tbxFile,chunk.size=NULL,return.type="data.frame")
      res=FUN(data,...)  
      
      # for text
      outfile= file.path(path.expand( dir),paste(chunk.num,filename,sep="_"))
      .write.table.noSci(res,outfile,quote=FALSE,col.names=FALSE,row.names=FALSE,
                  sep="\t")
    }
    
    # attach a random string to the file name 
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(rndFile,filename,sep="_")
    
    # apply function to chunks
    res=lapply(1:chunk.num,myFunc2,tbxFile,dir,filename2,FUN,...)
    
    
    outfile= file.path(path.expand(dir),filename) # get file name 
    if(file.exists(outfile)){
      message("overwriting ",outfile)
      unlink(outfile)
    }
    con=file(outfile, open = "a", blocking = TRUE) # open connection  
    for(file in gtools::mixedsort(
      list.files(path = dir, pattern = filename2,full.names=TRUE))){
      file.append(outfile,file) # append files
    }
    close(con)
    #remove temp files
    unlink(list.files(path = dir, pattern = filename2,full.names=TRUE))
    
    return(outfile)
  
  }else if(return.type=="data.frame"){
    
    # create a custom function that contains the function
    # to be applied
    myFunc3<-function(chunk.num,tbxFile,FUN,...){
      data=getTabixByChunk(tbxFile,chunk.size=NULL,return.type="data.frame")
      FUN(data,...)    
    }
    
    res=lapply(1:chunk.num,myFunc3,tbxFile,FUN,...)
    
    # collect and return
    data.frame(rbindlist(res))
  }else{
    
    myFunc4<-function(chunk.num,tbxFile,FUN,...){
      data=getTabixByChunk(tbxFile,chunk.size=NULL,return.type="data.table")
      FUN(data,...)    
    }
    
    res=lapply(1:chunk.num,myFunc4,tbxFile,FUN,...)
  
    
    # collect and return
    rbindlist(res)
  }
  
}


# applyTbxByCHr
#' Apply a function on tabix files chromosome by chromosome 
#' 
#' The function reads a tabix file chromosome by chromosome and applies 
#' a function on each. 
#' The function (FUN argument) should apply on data.frames and return a 
#' data frame
#' as a result. The function is parallel applied to each chromosome 
#' and related arguments could be passed on to the function via ... argument.
#' 
#' @param tbxFile tabix file to read. a TabixFile object
#' @param chrs chromosome names. Based on chromosome names the chunks of 
#' tabix file
#'        will be read into the memory. If missing use all chromosome names 
#'        in tabix file.
#' @param return.type indicates the return type for the function
#' @param FUN function to apply to chrs, it takes a data.frame and returns a 
#'            data.frame. First argument of the function should be a data frame
#' @param ... parameters to be passed to FUN. 
#' @param dir directory to create temporary files and the resulting tabix file
#' @param filename the filename for the resulting tabix file, this should not be 
#' a path, just a file name.
#' @param mc.cores number of cores to use in parallel (works only on UNIX based OS)
#' 
#' @return either a path to a tabix or text file, or a data frame or data.table
#' @noRd
applyTbxByChr<-function(tbxFile,chrs,dir,filename,
                        return.type=c("tabix","data.frame","data.table"),
                        FUN,...,mc.cores=1){
  
  return.type <- match.arg(return.type)
  FUN <- match.fun(FUN)
  if(Sys.info()['sysname']=="Windows") {mc.cores = 1}
  if(missing(chrs)) { chrs = Rsamtools::seqnamesTabix(tbxFile)}
  if(return.type =="tabix"){
   
    # create a custom function that contains the function
    # to be applied
    myFunc<-function(chr,tbxFile,dir,filename,FUN,...){
      data=getTabixByChr(chr = chr,tbxFile,return.type="data.frame")
      res=FUN(data,...)  
      
      # for tabix
      outfile= file.path(path.expand( dir),paste(chr,filename,sep="_"))
      
      .write.table.noSci(res,outfile,quote=FALSE,col.names=FALSE,row.names=FALSE,
                  sep="\t")
    }
    
    # attach a random string to the file name 
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(rndFile,filename,sep="_")
    
    res=mclapply(chrs,myFunc,tbxFile,dir,filename2,FUN,...,mc.cores = mc.cores)
    
    # collect & cat temp files,then make tabix

    path <- catsub2tabix(dir,filename2,filename)
    
    return(gsub(".tbi","",path))

    
  }else if(return.type=="data.frame"){
    
    # create a custom function that contains the function
    # to be applied
    myFunc2<-function(chr,tbxFile,FUN,...){
      data=getTabixByChr(chr = chr,tbxFile,return.type="data.frame")
      res=FUN(data,...)  
    }
    
    res=mclapply(chrs,myFunc2,tbxFile,FUN,...,mc.cores = mc.cores)
    
    # collect and return
    data.frame(rbindlist(res))
  }else{
    
    myFunc3<-function(chr,tbxFile,FUN,...){
      data=getTabixByChr(chr = chr,tbxFile,return.type="data.table")
      res=FUN(data,...)  
    }
    
    res=mclapply(chrs,myFunc3,tbxFile,FUN,...,mc.cores = mc.cores)
    
    # collect and return
    rbindlist(res)
  }
}

# applyTbxByOverlap
#' Serially apply a function on regions of tabix files
#' 
#' The function reads regions of a tabix file and applies a function on them. 
#' The function (FUN argument) should apply on data.frames and return a data frame
#' as a result. The function is serially applied to chunks (means no parallelization). 
#' However, the function FUN itself can be a parallelized function
#' and related arguments could be passed on to the function via ... argument.
#' 
#' @param tbxFile tabix file to read. a TabixFile object
#' @param ranges a GRanges object specifying the regions
#' @param chunk.size number of rows to be taken as a chunk, default: 1e6
#' @param return.type indicates the return type for the function
#' @param FUN function to apply to chunks, it takes a data.frame and returns a 
#'            data.frame. First argument of the function should be a data frame
#' @param ... parameters to be passed to FUN. 
#' @param dir directory to create temporary files and the resulting tabix file
#' @param filename the filename for the resulting tabix file, this should not be 
#' a path, just a file name.
#'
#' @return either a path to a tabix or text file, or a data frame or data.table
#' @noRd
applyTbxByOverlap<-function(tbxFile,ranges,chunk.size=1e6,dir,filename,
                          return.type=c("tabix","data.frame","data.table"),
                          FUN,...){
  
  return.type <- match.arg(return.type)
  FUN <- match.fun(FUN)
  
  # open tabix file with given chunk size
  if( class(tbxFile) != "TabixFile" ){
    tbxFile <- Rsamtools::TabixFile(tbxFile)
    
  } 
#   else {
#     if(Rsamtools::isOpen(tbxFile)){close(tbxFile)}# close if already open
#     Rsamtools::yieldSize(tbxFile) <-  chunk.size 
#   }
  
  
  # calculate number of chunks
  total.width = sum(as.numeric(width(ranges)))
  chunk.num=ceiling((total.width/width(ranges)[1])/chunk.size)
  groups <- ceiling(seq(length(ranges))/ceiling(length(ranges)/chunk.num))
  
  if( (chunk.num > length(ranges)) ){
    chunk.num <- length(ranges)
     groups <- seq(length(ranges))
  }
  #print(paste("chunks:",chunk.num))
  region.split <- split(ranges,groups)
  if(return.type =="tabix"){
    
    # create a custom function that contains the function
    # to be applied
    myFunc<-function(chunk.num,region.split,tbxFile,dir,filename,FUN,...){
      data <- try(expr = data <- getTabixByOverlap(
        tbxFile,granges = region.split[[chunk.num]],
        return.type="data.frame"),silent = TRUE)
      
      if( class(data)== "try-error") {
        
#         warning( paste("No records found in range between",
        # min(IRanges::end(region.split[[chunk.num]])),
#                        "and",max(IRanges::end(region.split[[chunk.num]])),
#                        "for Chromosome",
        # unique(as.character(region.split[[chunk.num]]@seqnames))))
        
      } else {
      
      res=FUN(data,...)  
      
      # for tabix
      outfile= file.path(path.expand( dir),paste(chunk.num,filename,sep="_"))
      .write.table.noSci(res,outfile,quote=FALSE,col.names=FALSE,row.names=FALSE,
                  sep="\t")
    }
    }
    
    # attach a random string to the file name 
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(rndFile,filename,sep="_")
    
    # apply function to chunks
    res=lapply(1:chunk.num,myFunc,region.split,tbxFile,dir,filename2,FUN,...)
    
    # collect & cat temp files,then make tabix
    path <- catsub2tabix(dir,pattern=filename2,filename,sort = TRUE)
    
    return(gsub(".tbi","",path))
    
  } else if(return.type=="data.frame"){
    
    # create a custom function that contains the function
    # to be applied
    myFunc2<-function(chunk.num,region.split,tbxFile,FUN,...){
      data <- try(expr = data <- getTabixByOverlap(
        tbxFile,granges = region.split[[chunk.num]],
        return.type="data.frame"),silent = TRUE)
      
      if( !(class(data)== "try-error") ) {
      res=FUN(data,...) 
      }
    }
    
    res=lapply(1:chunk.num,myFunc2,region.split,tbxFile,FUN,...)
    
    # collect and return
    data.frame(rbindlist(res))
  }else{
    
    myFunc3<-function(chunk.num,region.split,tbxFile,FUN,...){
      data <- try(expr = data <- getTabixByOverlap(
                                    tbxFile,
                                    granges = region.split[[chunk.num]],
                                    return.type="data.table"),silent = TRUE)
      
      if( !(class(data)[1] == "try-error") ) { ## class of data.table is both "data.table" and "data.frame
      res=FUN(data,...)  
      }
    }
    
    res=lapply(1:chunk.num,myFunc3,region.split,tbxFile,FUN,...)
    
    
    # collect and return
    rbindlist(res)
  }
  
}

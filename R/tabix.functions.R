# tabix functions
# these functions can be used to manipulate tabix files and methylKit objects
# that depend on those files


#' merge tabix files by chr, start,end, strand
#' 
#' @param tabixList list of tabix files created by Rsamtools::TabixFileList
#' @param dir working directory
#' @param filename the output file name
#' @param mc.cores number of multiple cores. If mc.cores>1 temporary files for each chromsome
#'        will be created prior to cat, zipping and indexing the single output file
#'        
#'        mergeTabix(tabixList,dir="~",filename="dump.meth.txt",mc.cores=1) 
mergeTabix<-function(tabixList,dir,filename,mc.cores=1 ){
  
  # get outfile
  
  #filename="dump.meth.txt"
  #dir="~" 
  
  
  # get chrs
  chrNum=table(unlist(lapply(tabixList,seqnamesTabix)))
  chrs=names(chrNum)[chrNum==length(tabixList)]
  
  if(mc.cores > 1){
    # random file string
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(rndFile,filename,sep="_")
    outlist=mclapply(chrs,mergeTbxByChr,tabixList,dir,filename2,parallel=TRUE,mc.cores=mc.cores)
    
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
    unlink(list.files(path = dir, pattern = rndFile,full.names=T))
  }else{
    
    outlist=mclapply(chrs,mergeTbxByChr,tabixList,dir,filename,parallel=FALSE,mc.cores=1)
    outfile= file.path(path.expand(dir),filename) 
  }
  
  #df=data.table::fread(outfile,header=FALSE)
  
  # after iterations end
  # zip and index the text file
  makeMethTabix( outfile ,skip=0)
}


#' make a tabix file from a data.frame or data.table part of a 
#' methylKit object
#' 
#' 
df2tabix<-function(df,outfile){
  
  # write the file to disk
  write.table(df,outfile,quote=FALSE,sep="\t",
              col.names=FALSE,row.names=FALSE)
  
  
  #make tabix out if the file
  makeMethTabix( outfile ,skip=0)
  
}


#' cat subfiles and make tabix
#' @param dir directory of subfiles to be merged
#' @param pattern a pattern in the file names
#' @param filename output file name
catsub2tabix<-function(dir,pattern,filename){
  
  outfile= file.path(path.expand(dir),filename) # get file name 
  if(file.exists(outfile)){
    message("overwriting ",outfile)
    unlink(outfile)
  }
  
  con=file(outfile, open = "a", blocking = TRUE) # open connection  
  for(file in list.files(path = dir, pattern = pattern,full.names=TRUE)){
    file.append(outfile,file) # append files
  }
  close(con)
  #remove temp files
  unlink(list.files(path = dir, pattern = pattern,full.names=T))
  
  #make tabix
  makeMethTabix( outfile ,skip=0)
  
}

# make tabix from a flat file where first 3 columns
# are chr,start,end in that order
makeMethTabix<-function(filepath,skip=0){
  message("compressing the file with bgzip...")
  zipped <- Rsamtools::bgzip(filepath,overwrite=TRUE)
  
  message("making tabix index...")
  Rsamtools::indexTabix(zipped,
             seq=1, start=2, end=3,
             skip=skip, comment="#", zeroBased=FALSE)
  
}


#' merge the data tables for a given chr
#' 
mergeTbxByChr<-function(chr,tabixList,dir,filename,parallel=FALSE){
  
  #get first file on the list
  res=getTabixByChr(chr,tabixList[[1]])
  for(i in 2:length(tabixList)){
    
    # get tabix per Chr
    tmp=getTabixByChr(chr,tabixList[[i]])
    
    # merge tabix in memory
    res=merge(res,tmp,by=c("V1","V2","V3"))
    
  }
  
  # order rows
  setorder(res, V1,V2,V3)
  
  if(!parallel){
  # write out append TRUE
  # if not parallel we can write out one file and append to it
  outfile= file.path(path.expand(dir),filename)   
  con=file(outfile, open = "a", blocking = TRUE)
  write.table(res,con,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
  close(con)
  }else{
    
    # when parallel, we have to write seperate files for each chr
    # then merge them outside this function
    outfile= file.path(path.expand( dir),paste(chr,filename,sep="_")) 
    write.table(res,outfile,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=TRUE)
    
  }
  return(outfile)
}



#' get data from meth tabix for a given chr
#'
#' @example
#' tbxFile="data/GSM1180321_methylC-Seq.mESC_TET2_KO.biorep1.techrep1.mkit.txt.bgz"
#'  getTabixByChr(chr="chr10",tbxFile)
getTabixByChr<-function(tbxFile,chr="chr10",return.type="data.table"){
  
  if( class(tbxFile) != "TabixFile" ){
    tbxFile <- Rsamtools::TabixFile(tbxFile)
  }
  
  res=Rsamtools:::.tabix_scan(tbxFile,space=chr,start=1L,end=500000000L) 
  
  if(return.type=="data.table")
  {
    tabix2dt(res)
  }else{
    tabix2df(res)
  }
}


#' get data from meth tabix for a given set of regions
#'
#' @example
#' tbxFile="data/GSM1180321_methylC-Seq.mESC_TET2_KO.biorep1.techrep1.mkit.txt.bgz"
#'  getTabixByChr(chr="chr10",tbxFile)
getTabixByOverlap<-function(tbxFile,granges,return.type="data.table"){
  
  if( class(tbxFile) != "TabixFile" ){
    tbxFile <- Rsamtools::TabixFile(tbxFile)
  }
  
  res=Rsamtools::scanTabix(tbxFile,param=granges) 
  
  if(return.type=="data.table")
  {
    tabix2dt(res)
  }else{
    tabix2df(res)
  }
}

#' get data from meth tabix for a given set of regions
#'
#' @example
#' tbxFile="data/GSM1180321_methylC-Seq.mESC_TET2_KO.biorep1.techrep1.mkit.txt.bgz"
#'  getTabixByChr(chr="chr10",tbxFile)
headTabix<-function(tbxFile,nrow=10,return.type="data.table"){
  
  if( class(tbxFile) != "TabixFile" ){
    tbxFile <- Rsamtools::TabixFile(tbxFile)
    open(tbxFile)
  }
  
  getTabixByChunk( tbxFile,chunk.size=nrow,return.type)
}


#' get data from meth tabix for a given chunkSize
#'
#' @example
#' tbxFile="data/GSM1180321_methylC-Seq.mESC_TET2_KO.biorep1.techrep1.mkit.txt.bgz"
#'  getTabixByChunk( tbxFile,chunk.size=10)
getTabixByChunk<-function(tbxFile,chunk.size=1e6,return.type="data.table"){
  
  if( class(tbxFile) != "TabixFile" | !Rsamtools::isOpen(tbxFile, rw="read") ){
    stop("tbxFile has to be a class of TabixFile and should be open for reading ")
  }
  
  if(is.na(Rsamtools::yieldSize(tbxFile)) | is.numeric(chunk.size)  ){
    Rsamtools::yieldSize(tbxFile)<-chunk.size
  }
  
  if(return.type=="data.table")
  {
    tabix2dt(Rsamtools::scanTabix(tbxFile) )
  }else{
    tabix2df(Rsamtools::scanTabix(tbxFile) )
  }
}


# convert methylKit tabix to data.table
# assuming you get a list length 1
tabix2dt<-function(tabixRes){
  dt <- data.table::fread( paste(tabixRes[[1]],collapse="\n" ) )
  setnames(x = dt,old = names(dt), new = c("chr","start","end","strand","coverage","numCs","numTs"))
  dt
}

# convert methylKit tabix to data.frame
# assuming you get a list length 1
tabix2df<-function(tabixRes){
  df <- data.table::fread( paste(tabixRes[[1]],collapse="\n" ),data.table=FALSE )
  setnames(x = df,old = names(df), new = c("chr","start","end","strand","coverage","numCs","numTs"))
  df
}

# convert methylKit tabix to data.frame
# assuming you get a list length 1
tabix2gr<-function(tabixRes){
}




# applyTbxByChunk
#' serially apply a function on chunks of tabix files
#' 
#' The function reads chunks and applies a function. 
#' The function (FUN argument) should apply on data.frames and return a data frame
#' as a result. The function is serially applied to chunks (means no paralleization)
#' However, the function FUN itseld can be parallelized function
#' and related arguments could be passed on to the function via ... argument.
#' 
#' @param tbxFile tabix file to read. a TabixFile object
#' @param chunk.size number of rows to be taken as a chunk, default:
#' @param return.type indicates the return type for the function
#' @param FUN function to apply to chunks, it takes a data.frame and returns a 
#'            data.frame. First argument of the function should be a data frame
#' @param ... parameters to be passed to FUN. 
#' @param dir directory to create temporary files and the resulting tabix file
#' @param filename the filename for the resulting tabix file, this should not be 
#' a path, just a file name.
#'
#' @return either a path to a tabix or text file, or a data frame or data.table
applyTbxByChunk<-function(tbxFile,chunk.size=1e6,dir,filename,
                          return.type=c("tabix","data.frame","data.table","text"),
                          FUN,...){
  
  return.type <- match.arg(return.type)
  FUN <- match.fun(FUN)
  
  # calculate number of chunks
  recs=countTabix(tbxFile$path)[[1]]
  chunk.num=ceiling(recs/chunk.size)
  
  
  if(return.type =="tabix"){
    
    # create a custom function that contains the function
    # to be applied
    myFunc<-function(chunk.num,tbxFile,dir,filename,FUN,...){
      data=getTabixByChunk(tbxFile,chunk.size=NULL,return.type="data.frame")
      res=FUN(data,...)  
      
      # for tabix
      outfile= file.path(path.expand( dir),paste(chunk.num,filename,sep="_"))
      write.table(res,outfile,quote=FALSE,col.names=FALSE,row.names=FALSE,
                  sep="\t")
    }
    
    
    # open tabix file with given chunk size
    if(isOpen(tbxFile)){close(tbxFile)}# close if already open
    open(tbxFile,yieldSize=chunk.size)
    
    # attach a random string to the file name 
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(rndFile,filename,sep="_")
    
    # apply function to chunks
    res=lapply(1:chunk.num,myFunc,dir,filename2,FUN,...)
    
    # collect & cat temp files,then make tabix
    catsub2tabix(dir,pattern=filename2,filename)
    
    return(file.path(path.expand( dir),filename) )
    
  }else if(return.type=="data.frame"){
    
    # create a custom function that contains the function
    # to be applied
    myFunc<-function(chunk.num,tbxFile,FUN,...){
      data=getTabixByChunk(tbxFile,chunk.size=NULL,return.type="data.frame")
      FUN(data,...)    
    }
    
    res=lapply(chunk.num,myFunc,tbxFile,FUN,...)
    
    # collect and return
    do.call("rbind",res)
  }else{
    
    myFunc<-function(chunk.num,tbxFile,FUN,...){
      data=getTabixByChunk(tbxFile,chunk.size=NULL,return.type="data.frame")
      FUN(data,...)    
    }
    
    res=lapply(chunk.num,myFunc,tbxFile,FUN,...)
  
    
    # collect and return
    data.table(do.call("rbind",res))
  }
  
  
  
  
  
  
}


#' apply a function on tabix files chromosome by chromosome 
#' 
#' @param tbxFile
#' @param chrs chromosome names. Based on chromosome names the chunks of tabix file
#'        will be read into the memory.
#' @param return.type indicates the return type for the function
#' @param FUN function to apply to chunks, it takes a data.frame and returns a 
#'            data.frame. First argument of the function should be a data frame
#' @param ... parameters to be passed to FUN. 
#' @param dir directory to create temporary files and the resulting tabix file
#' @param filename the filename for the resulting tabix file, this should not be 
#' a path, just a file name.

#' 
applyTbxByChr<-function(tbxFile,chrs,dir,filename,return.type=c("tabix","data.frame","data.table"),FUN,...,mc.cores){
  
  return.type <- match.arg(return.type)
  FUN <- match.fun(FUN)
  if(return.type =="tabix"){
   
    # create a custom function that contains the function
    # to be applied
    myFunc<-function(chr,tbxFile,dir,filename,FUN,...){
      data=getTabixByChr(chr,tbxFile,return.type="data.frame")
      res=FUN(data,...)  
      
      # for tabix
      outfile= file.path(path.expand( dir),paste(chr,filename,sep="_"))
      write.table(res,outfile,quote=FALSE,col.names=FALSE,row.names=FALSE,
                  sep="\t")
    }
    
    # attach a random string to the file name 
    rndFile=paste(sample(c(0:9, letters, LETTERS),9, replace=TRUE),collapse="")
    filename2=paste(rndFile,filename,sep="_")
    
    res=mclapply(chrs,myFunc,dir,filename2,...)
    
    # collect & cat temp files,then make tabix
    catsub2tabix(dir,pattern,filename2)
    
  }else if(return.type=="data.frame"){
    
    # create a custom function that contains the function
    # to be applied
    myFunc<-function(chr,tbxFile,dir,filename,...){
      data=getTabixByChr(chr,tbxFile,return.type="data.frame")
      res=FUN(data,...)  
    }
    
    res=mclapply(chrs,myFunc,tbxFile,dir,filename2,...)
    
    # collect & cat temp files,then make tabix
    do.call("rbind",res)
  }else{
    
    myFunc<-function(chr,tbxFile,dir,filename,...){
      data=getTabixByChr(chr,tbxFile,return.type="data.frame")
      res=FUN(data,...)  
    }
    
    res=mclapply(chrs,myFunc,tbxFile,dir,filename2,...)
    
    # collect & cat temp files,then make tabix
    data.table(do.call("rbind",res))
  }
  
 
  
 
  
}



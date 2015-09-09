
valid.methylRawDB <- function(object) {
  
  
  #data=getData(object,nrow=5)
  check1=( (object@resolution == "base") | (object@resolution == "region") )
  check2=file.exists(object@dbpath)
  if(check1 & check2){
    return(TRUE)
  }
  else if (! check1 ){
    message("resolution slot has to be either 'base' or 'region':",
            "other values not allowed")
    FALSE
  }
  else if(! check2){
    message("The DB file can not be found, check the value of 'dbpath'")
    FALSE
  }
  
}


# methylRawDB
# dbpath -> filepath to flat file database, in this case tabix file
# num.records -> number of lines in tabix file
# dbtype currently only "tabix" allowed
setClass("methylRawDB", slots=list(dbpath="character",num.records="numeric",
  sample.id = "character", assembly = "character",context="character",
  resolution="character",dbtype="character"),validity=valid.methylRawDB)





# PRIVATE function:
# makes a methylRawDB object from a df
# it is called from modRead function or whenever this functionality is needed
makeMethylRawDB<-function(df,dbpath,dbtype,
                          sample.id, assembly ,context,
                          resolution){
  filepath=paste0(dbpath,"/",sample.id,".txt")
  df2tabix(df,filepath)
  num.records=Rsamtools::countTabix(paste0(filepath,".bgz"))[[1]] ## 
  
  new("methylRawDB",dbpath=paste0(filepath,".bgz"),num.records=num.records,
  sample.id = sample.id, assembly = assembly,context=context,
  resolution=resolution,dbtype=dbtype)
}

setMethod("show", "methylRawDB", function(object) {
  
  cat("methylRawDB object with",object@dbtype,"rows\n--------------\n")
  #print(head(object))
  print(object@dbpath)
  cat("--------------\n")
  cat("sample.ids:",object@sample.ids,"\n")
  cat("destranded",object@destranded,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("treament:", object@treatment,"\n")
  cat("resolution:", object@resolution,"\n")
})


selectByOverlap<-function(object, ranges,inmemory=TRUE){
  
}

# methylRawListDB

valid.methylRawListDB <- function(object) {
  
  
  # if all elements are methyl
  if(!all(sapply(object,class)=="methylRawDB")){
    FALSE
  }
  else if ( length(object) != length(object@treatment) ){
    FALSE
  }
  else{
    TRUE
  }
  
}

setClass("methylRawListDB", slots=list(treatment = "vector"),contains = "list",
         validity=valid.methylRawListDB)


# methylBaseDB


# accesors



# subset each class


# get/set treatment vector

# get/set sample names
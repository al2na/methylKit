## batch control 



#' reconstruct methylBase object based on a new methylation percentage matrix
#' 
#' The function reconstructs a new methylBase object from another methylBase object
#' and percent methylation matrix. Basically, it uses the coverage columns in the given
#' methylBase object and deduces new number of methylated Cs and unmethylated Cs based
#' on given percent methylation matrix. It is ideally used to reconstruct methylBase objects
#' after batch correction on percent methylation values. The percent methylation matrix
#' rows must match methylBase object rows in order and in addition column order (the order of samples)
#' should match in both in given methylBase and  percent methylation matrix.
#' 
#' @param methMat percent methylation matrix, row order and order of the samples same as the methylBase object
#' @param mBase \code{\link{methylBase}} object to be reconstructed 
#' 
#' @return new \code{\link{methylBase}} object where methylation percentage matches
#'         input \code{methMat} and coverages matches input \code{mBase}
#' 
#' @author Altuna Akalin
#' 
#' 
#' @examples 
#' 
#' @export
#' @docType methods
#' @rdname reconstruct-methods
reconstruct<-function(methMat,mBase){
  
  # check if indeed methMat is percent methylation matrix
  if(max(methMat)<=1){
    warn("\nmake sure 'methMat' is percent methylation matrix (values between 0-100) \n")
  }

  # check if indeed methMat is percent methylation matrix  
  if(nrow(methMat) != nrow(mBase) | ncol(methMat) != length(mBase@numCs.index) ){
    stop("\nmethMat dimensions do not match number of samples\n",
         "and number of bases in methylBase object\n")
  }
  
  df=getData(mBase)
  mat=df[,mBase@numCs.index]+getData(mBase)[,mBase@numTs.index]
  
  # get new unmethylated and methylated counts
  numCs=round(methMat*mat/100)
  numTs=round((100-methMat)*mat/100)
  
  df[,mBase@numCs.index]=numCs
  df[,mBase@numTs.index]=numTs
  
  new("methylBase",df,sample.ids=sample.ids,
      assembly=mBase@assembly,context=mBase@context,
      treatment=mBase@treatment,coverage.index=mBase@coverage.index,
      numCs.index=mBase@numCs.index,numTs.index=mBase@numTs.index,
      destranded=mBase@destranded,resolution=mBase@resolution )
  
}

#' associates principal components with sample annotations
#' 
#' This function associates principal components
#' sampleAnnotation=data.frame(batch=c("a","a","b","b"))
assocComp<-function(mBase,sampleAnnotation){
  scale=TRUE
  center=TRUE
  mat=percMethylation(mBase)
  pr=prcomp(mat,scale.=scale,center=center)
  vars=100*pr$sdev**2/sum(pr$sdev**2)

  for(i in ncol(sampleAnnotation)){
    if(is.factor(sampleAnnotation[,i]) | is.character(sampleAnnotation[,i]) ){
      kruskal.test(split(pr$rotation[,1],sampleAnnotation[,i]) )
      
    }
  }
  
  list(pcs=pr$rotation,vars=vars,association=assoc)
}

#' remove principal components from a methylBase object
#' 
#' This function can remove a given principal componet from a given 
#' methylBase object. First, it calculates principal components from
#' percent methylation matrix and removes the given component(s), reconstructs
#' the methylation matrix then reconstructs number of methylated and unmethylated Cs per
#' position based on the reconstructed percent methylation matrix.
#' 
#' @param mBase \code{\link{methylBase}} object with no NA values, that means
#'               all bases should be covered in all samples.
#' @param comp vector of component numbers to be removed
#' 
#' @examples
#' 
#' data(methylKit)
#' 
#' # remove 1st principal component
#' newObj=removeComp(methylBase.obj,comp=1)
#' 
#' # remove 3rd and 4th  principal components
#' newObj=removeComp(methylBase.obj,comp=c(3,4))
#' 
#' @export
#' @docType methods
#' @rdname reconstruct-methods 
removeComp<-function(mBase,comp=NULL){
  if(is.na(comp) | is.null(comp)){
    stop("no component to remove\n")
  }
  
  if(any(comp > ncol() )){
    stop("'comp' elements can only take values between 1 and number of samples\n")
  }
  
  scale=TRUE
  center=TRUE
  mat=percMethylation(mBase)  
  mat=scale(mat,scale=scale,center=center)
  centers <- attr(mat,'scaled:center') 
  scales <- attr(mat,'scaled:scale') 
  pr=prcomp(mat,scale.=FALSE,center=FALSE)
  
  pr$rotation[,comp]=0
  res=pr$x %*% t(pr$rotation)
  
  res=(scale(res,center=(-centers/scales),scale=1/scales))
  attr(res,"scaled:center")<-NULL 
  attr(res,"scaled:scale")<-NULL 
  res[res>100]=100
  res[res<0]=0
  res
}
  
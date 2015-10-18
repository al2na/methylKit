## batch control 



#' Reconstruct methylBase object based on a new methylation percentage matrix
#' 
#' The function reconstructs a new methylBase object from an input methylBase object
#' and percent methylation matrix. Basically, it uses the read coverages in the input
#' methylBase object and deduces new number of methylated Cs and unmethylated Cs based
#' on the input percent methylation matrix. It is ideally to be used to 
#' reconstruct methylBase objects
#' after batch correction on percent methylation values. The percent methylation 
#' matrix rows must match methylBase object rows in order ,and in addition column 
#' order (the order of samples) in input methylBase must match the order in 
#' percent methylation matrix.
#' 
#' @param methMat percent methylation matrix, row order and order of the samples
#'  same as the methylBase object
#' @param mBase \code{\link{methylBase}} object to be reconstructed 
#' 
#' @return new \code{\link{methylBase}} object where methylation percentage matches
#'         input \code{methMat} and coverages matches input \code{mBase}
#' 
#' @author Altuna Akalin
#' 
#' @note Batch effect correction (if any batch effect exists) is a tricky issue. 
#' We provide some simple ways to deal with it 
#' (see \code{\link{assocComp}} and \code{\link{removeComp}} ), 
#' But if you can find other ways to correct for batch effects and want to create 
#' a methylBase object with the corrected percent methylation values, you can use this function.
#' 
#' @examples 
#' data(methylKit)
#' 
#' # get percent methylation
#' mat=percMethylation(methylBase.obj)
#' 
#' # do some changes in the matrix
#' # this is just a toy example
#' # ideally you want to correct the matrix
#' # for batch effects
#' mat[mat==100]=80
#' 
#' # reconstruct the methylBase from the corrected matrix
#' newobj=reconstruct(mat,methylBase.obj)
#' 
#' @export
#' @docType methods
#' @rdname reconstruct-methods
setGeneric("reconstruct", function(methMat,mBase) standardGeneric("reconstruct"))

#' @rdname reconstruct-methods
#' @aliases reconstruct,methylBase-method
setMethod("reconstruct",signature(mBase="methylBase"), function(methMat,mBase){
  
  # check if indeed methMat is percent methylation matrix
  if(max(methMat)<=1){
    warning("\nmake sure 'methMat' is percent methylation matrix (values between 0-100) \n")
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
  
  new("methylBase",df,sample.ids=mBase@sample.ids,
      assembly=mBase@assembly,context=mBase@context,
      treatment=mBase@treatment,coverage.index=mBase@coverage.index,
      numCs.index=mBase@numCs.index,numTs.index=mBase@numTs.index,
      destranded=mBase@destranded,resolution=mBase@resolution )
  
}
)

#' Associate principal components with sample annotations
#' 
#' This function associates principal components with sample annotations
#' such as age, gender, batch_id. Can be used to detect which batch effects
#' are associated with the variation in the methylation values.
#' 
#' @param mBase \code{\link{methylBase}} object with no NA values in the data part.
#' @param sampleAnnotation a data frame where columns are different annotations and 
#'                        rows are the samples, in the same order as in the methylBase object.
#' 
#' @return a named list of principal component matrix (named 'pcs'),% variation explained
#'         by principal compopents (named 'vars') and a p-value matrix showing association
#'         p-values between sample annotations and principal components (named 'association').
#'         
#' 
#' @author Altuna Akalin
#'
#' @examples
#' data(methylKit)
#' sampleAnnotation=data.frame(batch_id=c("a","a","b","b"),age=c(19,34,23,40))
#' as=assocComp(mBase=methylBase.obj,sampleAnnotation)
#' 
#' 
#' 
#' @export
#' @docType methods
#' @rdname assocComp-methods
setGeneric("assocComp", function(mBase,sampleAnnotation) standardGeneric("assocComp"))

#' @rdname assocComp-methods
#' @aliases assocComp,methylBase-method
setMethod("assocComp","methylBase", function(mBase,sampleAnnotation){
  scale=TRUE
  center=TRUE
  mat=percMethylation(mBase) # get matrix
  pr=prcomp(mat,scale.=scale,center=center) # get PCA
  vars=100*pr$sdev**2/sum(pr$sdev**2) # calc variation explained

  # get association p-values using different tests
  res=list()
  for(i in 1:ncol(sampleAnnotation)){
    
    # for factors do kruskal.wallis or wilcox test
    if(is.factor(sampleAnnotation[,i]) | is.character(sampleAnnotation[,i]) | is.logical(sampleAnnotation[,i])){
      annot=as.factor(sampleAnnotation[,i])
      res[[names(sampleAnnotation)[i]]]=apply(pr$rotation,2,function(x){ # cat(x)
                                      if(length(unique(annot))>2 ){
                                        kruskal.test(split(x,annot)  )$p.value
                                        }else{
                                        wilcox.test(split(x,annot)[[1]],split(x,annot)[[2]]  )$p.value 
                                        }
            })
      
    }else{# for factors do cor.test
      annot=  sampleAnnotation[,i]
      res[[names(sampleAnnotation)[i]]]=apply(pr$rotation,2,function(x){ # cat(x)

          cor.test(x,annot)$p.value
  
      })
    }
  }
  
  list(pcs=pr$rotation,vars=vars,association=do.call("rbind",res))
}
)

#' Remove principal components from a methylBase object
#' 
#' This function can remove a given principal componet from a given 
#' methylBase object. First, it calculates principal components from
#' percent methylation matrix and removes the given component(s), reconstructs
#' the methylation matrix then reconstructs number of methylated and unmethylated Cs per
#' position based on the reconstructed percent methylation matrix, and finally returns
#' a new \code{\link{methylBase}} object.
#' 
#' @param mBase \code{\link{methylBase}} object with no NA values, that means
#'               all bases should be covered in all samples.
#' @param comp vector of component numbers to be removed
#' 
#' @return new \code{\link{methylBase}} object
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
#' @rdname removeComp-methods
setGeneric("removeComp", function(mBase,comp=NULL) standardGeneric("removeComp"))

#' @rdname removeComp-methods
#' @aliases removeComp,methylBase-method
setMethod("removeComp","methylBase", function(mBase,comp){ 
  if(is.na(comp) | is.null(comp)){
    stop("no component to remove\n")
  }
  
  if(any(comp > length(mBase@sample.ids) )){
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
  reconstruct(res,mBase)
}
)
#' Pools replicates within groups to a single sample per group
#'
#' The function sums up coverage, numCs and numTs values within each group so one representative sample for each group will be created in a new methylBase object
#' @param obj  \code{methylBase} object with two groups or more and each group should have multiple samples
#' @param sample.ids  a character vector of new sample.ids ex:c("test","control"), should follow the same order as unique treatment vector,
#'                    and should be equal to the length of the unique treatment vector
#' @return a  \code{methylBase} object
#'
#' @usage pool(obj,sample.ids)
#' @author  Altuna Akalin
#' @export
#' @examples
#' library(methylKit)
#' data(methylKit)
#' 
#' # methylBase.obj has two groups, each group has two samples,
#' # the following function will pool the samples in each group
#' # so that each group will be represented by one pooled sample
#' pooled.methylBase=pool(methylBase.obj,sample.ids=c("test","control"))
#'
#' @docType methods
#' @rdname pool-methods
setGeneric("pool", function(obj,sample.ids) standardGeneric("pool"))

#' @rdname pool-methods
#' @aliases pool,methylBase-method
setMethod("pool", "methylBase",
                    function(obj,sample.ids){
                      df=getData(obj)

                      treat=unique(obj@treatment)
                      res=df[,1:5]
                      for(i in 1:length(treat) ){

                         # get indices
                         setCs=obj@numCs.index[obj@treatment==treat[i]]
                         setTs=setCs+1
                         set.cov=setCs-1

                         if(length(setCs)>1){
                           Cs=rowSums(df[,setCs],na.rm=TRUE)
                           Ts=rowSums(df[,setTs],na.rm=TRUE)
                           covs=rowSums(df[,set.cov],na.rm=TRUE)

                         }else{
                           Cs  =df[,setCs]
                           Ts  =df[,setTs]
                           covs=df[,set.cov]
                         }
                         res=cbind(res,covs,Cs,Ts) # bind new data
                         names(res)[(3:5)+3*i]=paste( c("coverage" ,"numCs" ,"numTs"),i,sep="") # change names of columns

                      }
                      coverage.ind=3*(1+(1:length(treat)))
                      obj1=new("methylBase",as.data.frame(res),sample.ids=sample.ids,
                             assembly=obj@assembly,context=obj@context,
                             treatment=treat,coverage.index=coverage.ind,
                             numCs.index=coverage.ind+1,numTs.index=coverage.ind+2,destranded=obj@destranded,resolution=obj@resolution )
                      obj1

                      
                      
})

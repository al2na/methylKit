# The function to reorganize the methylRawList and methylBase objects
# 
# 

#' reorganize methylRawList and methylBase objects by creating new objects from subset of samples
#' 
#' Create a new  \code{methylRawList} or \code{methylBase} object by selecting a subset of samples from the input object, which is 
#' a \code{methylRawList} or \code{methylBase} object. You can use the function to partition a large methylRawList or methylBase object
#' to smaller object based on sample ids or when you want to reorder samples and treatmet vector.
#'
#' @param methylObj a \code{methylRawList} or \code{methylBase} object
#' @param sample.ids a vector for sample.ids to be subset. Order is important and the order should be similar to treatment. sample.ids should be
#'        a subset or reordered version of sample ids in the input object.
#' @param treatment  treatment vector, should be same length as sample.ids vector
#'
#' @return RETURNS a \code{methylRawList} or \code{methylBase} object depending on the input object
#' @usage reorganize(methylObj,sample.ids,treatment)
#' @examples
#' # this is a list of example files, ships with the package
#' file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
#'                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
#'
#'
#' # read the files to a methylRawList object: myobj
#' myobj=read( file.list,
#'           sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))
#' meth=unite(myobj,destrand=TRUE)
#'
#' myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
#' meth2 =reorganize(meth,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
#'
#' @export
#' @docType methods
#' @rdname reorganize-methods
setGeneric("reorganize", function(methylObj,sample.ids,treatment) standardGeneric("reorganize") )


#' @rdname reorganize-methods
#' @aliases reorganize,methylBase-method
setMethod("reorganize", signature(methylObj="methylBase"),
                    function(methylObj,sample.ids,treatment){
                      
                        #sample.ids length and treatment length should be equal
                        if(length(sample.ids) != length(treatment) ){
                          stop("length of sample.ids should be equal to treatment")
                        }
                        
                        if( ! all(sample.ids %in% methylObj@sample.ids) ){
                          stop("provided sample.ids is not a subset of the sample ids of the object")
                        }
          
                        temp.id = methylObj@sample.ids # get the subset of ids
                        col.ord = order(match(temp.id,sample.ids))[1:length(sample.ids)] # get the column order in the original matrix
                        
                        
                        ind.mat=rbind(methylObj@coverage.index[col.ord],  # make a matrix indices for easy access 
                                      methylObj@numCs.index[col.ord],
                                      methylObj@numTs.index[col.ord])
                        dat    =getData(methylObj) # get data
                        
                        #newdat =cbind(dat[,1:5],dat[ind.mat[,1]],dat[ind.mat[,2]],dat[ind.mat[,3]]) # reorder columns using ind.mat
                        newdat =dat[,1:5]
                        for(i in 1:ncol(ind.mat))
                        {
                          newdat=cbind(newdat,dat[,ind.mat[,i]])
                        }
                                                
                        # get indices of coverage,numCs and numTs in the data frame 
                        coverage.ind=seq(6,by=3,length.out=length(sample.ids))
                        numCs.ind   =coverage.ind+1
                        numTs.ind   =coverage.ind+2
                        
                        # update column names
                        #colnames(newdat)[6:ncol(newdat)]=paste(c("coverage","numCs","numTs"),rep(1:length(sample.ids),each=length(sample.ids)),sep="")                        

                        names(newdat)[coverage.ind]=paste(c("coverage"),1:length(sample.ids),sep="" )
                        names(newdat)[numCs.ind]   =paste(c("numCs"),1:length(sample.ids),sep="" )
                        names(newdat)[numTs.ind]   =paste(c("numTs"),1:length(sample.ids),sep="" )
                        

                        new("methylBase",as.data.frame(newdat),sample.ids=sample.ids,
                             assembly=methylObj@assembly,context=methylObj@context,
                             treatment=treatment,coverage.index=coverage.ind,
                             numCs.index=numCs.ind,numTs.index=numTs.ind,
                             destranded=methylObj@destranded, resolution=methylObj@resolution )
 
})


#' @rdname reorganize-methods
#' @aliases reorganize,methylRawList-method
setMethod("reorganize", signature(methylObj="methylRawList"),
                    function(methylObj,sample.ids,treatment){
                      
                        #sample.ids length and treatment length should be equal
                        if(length(sample.ids) != length(treatment) ){
                          stop("length of sample.ids should be equal to treatment")
                        }
                        
                        orig.ids=sapply(methylObj,function(x) x@sample.id) # get ids from the list of methylRaw 
                        if( ! all(sample.ids %in% orig.ids) ){
                          stop("provided sample.ids is not a subset of the sample ids of the object")
                        }
          
                        col.ord=order(match(orig.ids,sample.ids))[1:length(sample.ids)] # get the column order in the original matrix
          
                        outList=list()    
                        for(i in 1:length(sample.ids)){
                          #ind=which( orig.ids==sample.ids[i] )
                          outList[[i]]=methylObj[[ col.ord[i]  ]]
                          #outList[[i]]=methylObj[[ ind  ]]

                        }
          
                        new("methylRawList",outList,treatment=treatment)
          
})
          
          
          

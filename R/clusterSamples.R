# methylKit R package
# sampleCluster.R function

#---------------------------------------------------------------------------------------
# regular R functions to be used in S4 functions

# correlation function return dist object
# x matrix, each row is a sample.
# method only support "pearson"
# abs logical. If use absolute value for correlation to distance.
# diag logical. Inlcude diagnol value or not
# upper output upper panel or lower panel.
.dist.cor=function(x, method="pearson", abs=TRUE, diag=FALSE, upper=FALSE)
{
  if (!is.na(pmatch(method, "pearson"))) 
    method <- "pearson"
  METHODS <- c("pearson")
  method <- pmatch(method, METHODS)
  if (is.na(method)) 
    stop("invalid correlation method")
  if (method == -1) 
    stop("ambiguous correlation method")
  
  xcor = cor(t(x), method=METHODS[method])
  if(abs)
    xcor = 1-abs(xcor)
  else
    xcor = 1-xcor
  if(upper)
    d <- xcor[upper.tri(xcor,diag=diag)]
  else
    d <- xcor[lower.tri(xcor,diag=diag)]
  attr(d, "Size") <- nrow(x)
  attr(d, "Labels") <- dimnames(x)[[1L]]
  attr(d, "Diag") <- diag
  attr(d, "Upper") <- upper
  attr(d, "method") <- METHODS[method]
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  return(d)
} 

# cluster function on matrix and return hierarchical plot
# x matrix each column is a sample
# dist.method method to get the distance between samples
# hclust.method the agglomeration method to be used
# plot if TRUE, plot the hierarchical clustering
.cluster=function(x, dist.method="correlation", hclust.method="ward", plot=TRUE){
  DIST.METHODS <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", 
        "binary", "minkowski")
  dist.method <- pmatch(dist.method, DIST.METHODS)

  HCLUST.METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
        "median", "centroid")
  hclust.method <- pmatch(hclust.method, HCLUST.METHODS)
  if (is.na(hclust.method)) 
    stop("invalid clustering method")
  if (hclust.method == -1) 
    stop("ambiguous clustering method")

  if(DIST.METHODS[dist.method] == "correlation")
    d = .dist.cor(t(x))
  else
    d=dist(scale(t(x)), method=DIST.METHODS[dist.method]);
  
  hc=hclust(d, HCLUST.METHODS[hclust.method]);
  
  if(plot){
    plclust(hc,hang=-1, main=paste("CpG dinucleotide methylation clustering\nDistance: ",
                                   DIST.METHODS[dist.method],sep=""), xlab = "Samples");
  }
  return(hc)
  }

# Principal Components Analysis on methylBase object
# x matrix each column is a sample
# cor a logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. (The correlation matrix can only be used if there are no constant variables.)
.pcaPlot = function(x, cor=TRUE){
  x.pr = princomp(x, cor=cor)
  loads = loadings(x.pr)
  plot(loads[,1:2], main = "CpG dinucleotide methylation PCA Analysis")
  text(loads[,1], loads[,2],adj=c(-0.4,0.3))
  return(summary(x.pr))
}

# end of regular functions to be used in S4 functions
#---------------------------------------------------------------------------------------

#' Hierarchical cluster analysis on samples in methylBase object
#' 
#' @param .Object a \code{methylBase} object
#' @param dist the distance measure to be used. This must be one of "\code{correlation}", "\code{euclidean}", "\code{maximum}", "\code{manhattan}", "\code{canberra}", "\code{binary}" or "\code{minkowski}". Any unambiguous substring can be given. (default:"\code{correlation}")
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "\code{ward}", "\code{single}", "\code{complete}", "\code{average}", "\code{mcquitty}", "\code{median}" or "\code{centroid}". (default:"\code{ward}")
#' @param plot clustering plot if TRUE (default:TRUE) 
#' @return a \code{tree} object produced by hclust and plot hierarchical clustering
#' @aliases clusterSamples,-methods clusterSamples, methylBase-method
#' @export
#' @docType methods
#' @rdname clusterSamples-methods
setGeneric("clusterSamples", function(.Object, dist="correlation", method="ward", plot=TRUE) standardGeneric("clusterSamples"))

#' @rdname clusterSamples-methods
#' @aliases clusterSamples,ANY-method
setMethod("clusterSamples", "methylBase",
                    function(.Object, dist="correlation", method="ward", plot=TRUE){
                        meth.mat = getData(.Object)[, .Object@numCs.index]/(.Object[,.Object@numCs.index] + .Object[,.Object@numTs.index] )                                      
                        names(meth.mat)=.Object@sample.ids
                        
                        .cluster(meth.mat, dist.method=dist, hclust.method=method, plot=TRUE)
                        
                        }
)

#' Principal Components Analysis on samples in methylBase object
#' 
#' @param .Object a \code{methylBase} object
#' @param cor a logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. (The correlation matrix can only be used if there are no constant variables.)
#' @return The form of the value returned by \code{PCASamples} is the table of importance of components
#' @aliases PCASamples,-methods PCASamples, methylBase-method
#' @export
#' @docType methods
#' @rdname PCASamples-methods
setGeneric("PCASamples", function(.Object, cor=TRUE) standardGeneric("PCASamples"))

#' @rdname PCASamples-methods
#' @aliases PCASamples,ANY-method
setMethod("PCASamples", "methylBase",
                    function(.Object, cor=TRUE){
                        meth.mat = getData(.Object)[, .Object@numCs.index]/(.Object[,.Object@numCs.index] + .Object[,.Object@numTs.index] )                                      
                        names(meth.mat)=.Object@sample.ids
                        
                        .pcaPlot(meth.mat, cor=TRUE)
                        
                        }
)

#numC=grep("numCs", names(methidh))
#numT=grep("numTs", names(methidh))
#meth.mat=methidh[,numC]/(methidh[,numC]+methidh[,numT])


## ACESSOR FUNCTIONS FOR methylDiff OBJECT
# a function for getting data part of methylDiff                      
setMethod(f="getData", signature="methylDiff", definition=function(x) {
                return(as(x,"data.frame"))
        }) 

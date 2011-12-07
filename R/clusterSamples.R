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
.cluster=function(x, dist.method="correlation", hclust.method="ward", plot=TRUE, treatment=treatment,sample.ids=sample.ids){
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
    #plclust(hc,hang=-1, main=paste("CpG dinucleotide methylation clustering\nDistance method: ",
    #                               DIST.METHODS[dist.method],sep=""), xlab = "Samples");
    # plot
    treatment=treatment
    sample.ids=sample.ids

    my.cols=rainbow(length(unique(treatment)))

    col.list=as.list(my.cols[treatment+1])
    names(col.list)=sample.ids

    colLab <- function(n,col.list)
      {
      if(is.leaf(n))
        {
          a <- attributes(n)

          attr(n, "nodePar") <- c(a$nodePar, list(lab.col =
          col.list[[a$label]], lab.cex=1,
          col=col.list[[a$label]], cex=1, pch=16 ))
        }
      n
      }
    
    dend = as.dendrogram(hc)
    dend_colored <- dendrapply(dend, colLab,col.list)
    
    plot(dend_colored, main = "CpG dinucleotide methylation clustering", 
         sub = paste("Distance method: \"", DIST.METHODS[dist.method],
         "\"; Clustering method: \"", HCLUST.METHODS[hclust.method],"\"",sep=""), 
         xlab = "Samples", ylab = "Height");
    # end of plot
    }
  return(hc)
  }

# Adjust the range of a vector "x" by "i" proportion of the range of x
# 
# @param x a vector of value
# @param i the proprotion of the range of x to adjust the size.
.adjlim=function(x, i)
  {
  if(length(x)>1) {
    xr=range(x); 
    xlim=c(); 
    xlim[1]=xr[1] - abs(xr[1]) * i; 
    xlim[2]=xr[2] + abs(xr[2] * i); 
    return(xlim)} 
  else 
    print("length vector x should be more than 1")
  }
# Principal Components Analysis on methylBase object
# x matrix each column is a sample
# cor a logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. (The correlation matrix can only be used if there are no constant variables.)
.pcaPlot = function(x, cor=TRUE, screeplot=FALSE, adj.lim=c(0.001,0.1), treatment=treatment,sample.ids=sample.ids){
  x.pr = princomp(x, cor=cor)
  if (screeplot){
    i=5;screeplot(x.pr, type="barplot", main="CpG dinucleotide methylation PCA Screeplot", col = rainbow(i)[i])
  }
  else{
    loads = loadings(x.pr)

    treatment=treatment
    sample.ids=sample.ids
    my.cols=rainbow(length(unique(treatment)))
    
    plot(loads[,1:2], main = "CpG dinucleotide methylation PCA Analysis",col=my.cols[treatment+1],
         xlim=.adjlim(loads[,1],adj.lim[1]), ylim=.adjlim(loads[,2], adj.lim[2]))
    text(loads[,1], loads[,2],labels=sample.ids,adj=c(-0.4,0.3), col=my.cols[treatment+1])
  }
  return(summary(x.pr))
}

# end of regular functions to be used in S4 functions
#---------------------------------------------------------------------------------------

#' CpG Dinucleotide Methylation Hierarchical Cluster Analysis
#' 
#' @param .Object a \code{methylBase} object
#' @param dist the distance measure to be used. This must be one of "\code{correlation}", "\code{euclidean}", "\code{maximum}", "\code{manhattan}", "\code{canberra}", "\code{binary}" or "\code{minkowski}". Any unambiguous substring can be given. (default:"\code{correlation}")
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "\code{ward}", "\code{single}", "\code{complete}", "\code{average}", "\code{mcquitty}", "\code{median}" or "\code{centroid}". (default:"\code{ward}")
#' @param plot a logical value indicating whether to plot hierarchical clustering. (default:TRUE) 
#'
#' @return a \code{tree} object of a hierarchical cluster analysis using a set of dissimilarities for the n objects being clustered.
#'
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
                        
                        .cluster(meth.mat, dist.method=dist, hclust.method=method, plot=plot, treatment=.Object@treatment,sample.ids=.Object@sample.ids)
                        
                        }
)

#' CpG Dinucleotide Methylation Principal Components Analysis
#' 
#' @param .Object a \code{methylBase} object
#' @param cor a logical value indicating whether the calculation should use the correlation matrix or the covariance matrix. (default: TRUE)
#' @param screeplot a logical value indicating whether to plot the variances against the number of the principal component. (default: FALSE)
#' @param adj.lim a vector indicating the propotional adjustment of xlim (adj.lim[1]) and ylim (adj.lim[2]). (default: c(0.0004,0.1))
#'
#' @return The form of the value returned by \code{PCASamples} is the summary of principal component analysis by \code{princomp}.
#'
#' @aliases PCASamples,-methods PCASamples, methylBase-method
#' @export
#' @docType methods
#' @rdname PCASamples-methods
setGeneric("PCASamples", function(.Object, cor=TRUE, screeplot=FALSE, adj.lim=c(0.0004,0.1)) standardGeneric("PCASamples"))

#' @rdname PCASamples-methods
#' @aliases PCASamples,ANY-method
setMethod("PCASamples", "methylBase",
                    function(.Object, cor=TRUE, screeplot=FALSE, adj.lim=c(0.0004,0.1)){
                        meth.mat = getData(.Object)[, .Object@numCs.index]/(.Object[,.Object@numCs.index] + .Object[,.Object@numTs.index] )                                      
                        names(meth.mat)=.Object@sample.ids
                        
                        .pcaPlot(meth.mat, cor=cor, screeplot=screeplot, adj.lim=adj.lim, treatment=.Object@treatment,sample.ids=.Object@sample.ids)
                        
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

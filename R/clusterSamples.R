# methylKit R package
# sampleCluster.R function

#---------------------------------------------------------------------------------------
# regular R functions to be used in S4 functions




rowSds <- function(x, center=NULL, ...) {
  n <- !is.na(x);
  n <- rowSums(n);
  n[n <= 1] <- NA;

  if (is.null(center)) {
    center <- rowMeans(x, ...);
  }

  x <- x - center;
  x <- x*x;
  x <- rowSums(x, ...);
  x <- x/(n-1);

 sqrt(x);
}


colSds <- function(x, ...) {
  x <- t(x);
  rowSds(x, ...);
}




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
.cluster=function(x, dist.method="correlation", hclust.method="ward", plot=TRUE,
                  treatment=treatment,sample.ids=sample.ids,context){
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

    my.cols=rainbow(length(unique(treatment)), start=1, end=0.6)

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
    
    plot(dend_colored, main = paste(context, "methylation clustering"), 
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
.pcaPlot = function(x,comp1=1,comp2=2, screeplot=FALSE, adj.lim=c(0.001,0.1), treatment=treatment,sample.ids=sample.ids,context,scale=TRUE,center=TRUE,obj.return=FALSE){
  #x.pr = princomp(x, cor=cor)
  

  x.pr = prcomp((x),scale.=scale,center=center)

  if (screeplot){
    i=5;screeplot(x.pr, type="barplot", main=paste(context,"methylation PCA Screeplot"), col = rainbow(i)[i])
  }
  else{
    #loads = loadings(x.pr)
    loads = x.pr$rotation
    treatment=treatment
    sample.ids=sample.ids
    my.cols=rainbow(length(unique(treatment)), start=1, end=0.6)

    
    plot(loads[,comp1],loads[,comp2], main = paste(context,"methylation PCA Analysis"),
         col=my.cols[treatment+1],
         xlim=.adjlim(loads[,comp1],adj.lim[1]), 
         ylim=.adjlim(loads[,comp2], adj.lim[2]),
         xlab=paste("loadings for PC",comp1,sep=""), 
         ylab=paste("loadings for PC",comp2,sep=""))
    
    text(loads[,comp1], loads[,comp2],labels=sample.ids,adj=c(-0.4,0.3), 
         col=my.cols[treatment+1])
  }
  if(obj.return){  return((x.pr))}

}


# Principal Components Analysis on methylBase object on transposed data
# x matrix each column is a sample
.pcaPlotT = function(x,comp1=1,comp2=2,screeplot=FALSE, adj.lim=c(0.001,0.1),
                     treatment=treatment,sample.ids=sample.ids,context,
                     scale=TRUE,center=TRUE,obj.return=FALSE){

  x.pr = prcomp(t(x),scale.=scale,center=center)
 
  if (screeplot){
    i=5;screeplot(x.pr, type="barplot", 
                  main=paste(context,"methylation PCA Screeplot"), 
                  col = rainbow(i)[i])
  }
  else{
    #loads = loadings(x.pr)
    #loads = x.pr$rotation
    treatment=treatment
    sample.ids=sample.ids
    my.cols=rainbow(length(unique(treatment)), start=1, end=0.6)
    pc1=x.pr$x[,comp1]
    pc2=x.pr$x[,comp2]
    
    plot(pc1,pc2, main = paste(context,"methylation PCA Analysis"),
         col=my.cols[treatment+1],
         xlim=.adjlim(pc1,adj.lim[1]), ylim=.adjlim(pc2, adj.lim[2]),
         xlab=paste("PC",comp1,sep=""), ylab=paste("PC",comp2,sep=""))
    text(pc1, pc2,labels=sample.ids,adj=c(-0.4,0.3), col=my.cols[treatment+1])
  }
  if(obj.return){  return((x.pr))}
}


# end of regular functions to be used in S4 functions
#---------------------------------------------------------------------------------------

#' Hierarchical Clustering using methylation data
#'  
#' The function clusters samples using \code{\link[stats]{hclust}} function and 
#' various distance metrics derived from percent methylation per base or per 
#' region for each sample.
#' 
#' 
#' @param .Object a \code{methylBase} or \code{methylBaseDB} object
#' @param dist the distance measure to be used. This must be one of 
#'        "\code{correlation}", "\code{euclidean}", "\code{maximum}", 
#'        "\code{manhattan}", "\code{canberra}", "\code{binary}" or "\code{minkowski}". 
#'        Any unambiguous abbreviation can be given. (default:"\code{correlation}")
#' @param method the agglomeration method to be used. This should be 
#'        (an unambiguous abbreviation of) one of "\code{ward}", "\code{single}", 
#'        "\code{complete}", "\code{average}", "\code{mcquitty}", "\code{median}" 
#'        or "\code{centroid}". (default:"\code{ward}")
#' @param sd.filter  If \code{TRUE}, the bases/regions with low variation will be 
#'        discarded prior to clustering (default:TRUE)
#' @param sd.threshold A numeric value. If \code{filterByQuantile} is \code{TRUE}, 
#'        features whose standard deviations is less than the quantile denoted by 
#'        \code{sd.threshold} will be removed.
#'         If \code{filterByQuantile} is \code{FALSE}, then features whose 
#'         standard deviations is less than the value of \code{sd.threshold} 
#'         will be removed.(default:0.5)
#' @param filterByQuantile A logical determining if \code{sd.threshold} is to 
#'        be interpreted as a quantile of all Standard Deviation values from 
#'        bases/regions (the default), or as an absolute value
#' @param plot a logical value indicating whether to plot hierarchical 
#'        clustering. (default:TRUE) 
#' @param chunk.size Number of rows to be taken as a chunk for processing the \code{methylBaseDB} objects, default: 1e6
#'        
#' @usage clusterSamples(.Object, dist="correlation", method="ward",
#'                        sd.filter=TRUE,sd.threshold=0.5,
#'                        filterByQuantile=TRUE, plot=TRUE)
#' @examples
#' data(methylKit)
#' 
#' clusterSamples(methylBase.obj, dist="correlation", method="ward", plot=TRUE)
#' 
#' 
#' 
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with \code{methylBaseDB}  objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the \code{chunk.size}.
#' 
#' @return a \code{tree} object of a hierarchical cluster analysis using a set 
#'          of dissimilarities for the n objects being clustered.
#'
#' @export
#' @docType methods
#' @rdname clusterSamples-methods
setGeneric("clusterSamples", function(.Object, dist="correlation", method="ward",
                                      sd.filter=TRUE,sd.threshold=0.5,
                                      filterByQuantile=TRUE, plot=TRUE,chunk.size=1e6) 
                                              standardGeneric("clusterSamples"))

#' @rdname clusterSamples-methods
#' @aliases clusterSamples,methylBase-method
setMethod("clusterSamples", "methylBase",
  function(.Object, dist, method ,sd.filter, sd.threshold, 
                   filterByQuantile, plot)
  {
    mat      =getData(.Object)
    # remove rows containing NA values, they might be introduced at unite step
    mat      =mat[ rowSums(is.na(mat))==0, ] 
    
    meth.mat = mat[, .Object@numCs.index]/
      (mat[,.Object@numCs.index] + mat[,.Object@numTs.index] )                                      
    names(meth.mat)=.Object@sample.ids
    
    # if Std. Dev. filter is on remove rows with low variation
    if(sd.filter){
      if(filterByQuantile){
        sds=rowSds(as.matrix(meth.mat))
        cutoff=quantile(sds,sd.threshold)
        meth.mat=meth.mat[sds>cutoff,]
      }else{
        meth.mat=meth.mat[rowSds(as.matrix(meth.mat))>sd.threshold,]
      }
    }
    
    .cluster(meth.mat, dist.method=dist, hclust.method=method, 
             plot=plot, treatment=.Object@treatment,
             sample.ids=.Object@sample.ids,
             context=.Object@context)
    
  }
)

#' Principal Components Analysis of Methylation data
#' 
#' The function does a PCA analysis using \code{\link[stats]{prcomp}} function 
#' using percent methylation matrix as an input.
#' 
#' @param .Object a \code{methylBase} or \code{methylBaseDB} object
#' @param screeplot a logical value indicating whether to plot the variances 
#'        against the number of the principal component. (default: FALSE)
#' @param adj.lim a vector indicating the propotional adjustment of 
#'        xlim (adj.lim[1]) and 
#'        ylim (adj.lim[2]). This is primarily used for adjusting the visibility
#'        of sample labels on the on the PCA plot. (default: c(0.0004,0.1))
#' @param scale logical indicating if \code{prcomp} should scale the data to 
#'        have unit variance or not (default: TRUE)
#' @param center logical indicating if \code{prcomp} should center the data 
#'        or not (default: TRUE)
#' @param comp vector of integers with 2 elements specifying which components 
#'        to be plotted.
#' @param transpose if TRUE (default) percent methylation matrix will be 
#'        transposed, this is equivalent to doing PCA on variables that are 
#'        regions/bases. The resulting plot will location of samples in the new 
#'        coordinate system if FALSE the variables for the matrix will be samples 
#'        and the resulting plot whill show how each sample (variable) 
#'        contributes to the principle component.the samples that are highly 
#'        correlated should have similar contributions to the principal components.       
#' @param sd.filter  If \code{TRUE}, the bases/regions with low variation will 
#'        be discarded prior to PCA (default:TRUE)
#' @param sd.threshold A numeric value. If \code{filterByQuantile} is \code{TRUE}, 
#'        the value should be between 0 and 1 and the features whose standard 
#'        deviations is less than the quantile denoted by \code{sd.threshold} 
#'        will be removed. If \code{filterByQuantile} is \code{FALSE}, 
#'        then features whose standard deviations is less than the value 
#'        of \code{sd.threshold} will be removed.(default:0.5)
#' @param filterByQuantile A logical determining if \code{sd.threshold} is to be 
#'        interpreted as a quantile of all standard deviation values from 
#'        bases/regions (the default), or as an absolute value
#' @param obj.return if the result of \code{prcomp} function should be returned 
#'        or not. (Default:FALSE)
#' @param chunk.size Number of rows to be taken as a chunk for processing the \code{methylRawListDB} objects, default: 1e6
#' 
#' @usage PCASamples(.Object, screeplot=FALSE, adj.lim=c(0.0004,0.1), scale=TRUE,
#' center=TRUE,comp=c(1,2),transpose=TRUE,sd.filter=TRUE,
#'            sd.threshold=0.5,filterByQuantile=TRUE,obj.return=FALSE)
#' 
#' @examples
#' data(methylKit) 
#' 
#' # do PCA with filtering rows with low variation, filter rows with standard 
#' # deviation lower than the 50th percentile of Standard deviation distribution
#' PCASamples(methylBase.obj,screeplot=FALSE, adj.lim=c(0.0004,0.1),
#'            scale=TRUE,center=TRUE,comp=c(1,2),transpose=TRUE,sd.filter=TRUE,
#'            sd.threshold=0.5,filterByQuantile=TRUE,obj.return=FALSE)
#' 
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with \code{methylBaseDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the \code{chunk.size}.
#' 
#' @return The form of the value returned by \code{PCASamples} is the summary 
#'         of principal component analysis by \code{prcomp}.
#' @note cor option is not in use anymore, since \code{prcomp} is used for PCA 
#'        analysis instead of \code{princomp}
#'  
#'
#' @export
#' @docType methods
#' @rdname PCASamples-methods
setGeneric("PCASamples", function(.Object, screeplot=FALSE, adj.lim=c(0.0004,0.1),
                                  scale=TRUE,center=TRUE,comp=c(1,2),transpose=TRUE,
                                  sd.filter=TRUE,sd.threshold=0.5,
                                  filterByQuantile=TRUE,obj.return=FALSE,chunk.size=1e6) 
          standardGeneric("PCASamples"))

#' @rdname PCASamples-methods
#' @aliases PCASamples,methylBase-method
setMethod("PCASamples", "methylBase",
  function(.Object, screeplot, adj.lim,scale,center,comp,
                             transpose,sd.filter, sd.threshold, 
                             filterByQuantile,obj.return)
  {
    
    mat      = getData(.Object)
    # remove rows containing NA values, they might be introduced at unite step
    mat      = mat[ rowSums(is.na(mat))==0, ] 
    meth.mat = mat[, .Object@numCs.index]/
      (mat[,.Object@numCs.index] + mat[,.Object@numTs.index] )                                      
    names(meth.mat)=.Object@sample.ids
    
    # if Std. Dev. filter is on remove rows with low variation
    if(sd.filter){
      if(filterByQuantile){
        sds=rowSds(as.matrix(meth.mat))
        cutoff=quantile(sds,sd.threshold)
        meth.mat=meth.mat[sds>cutoff,]
      }else{
        meth.mat=meth.mat[rowSds(as.matrix(meth.mat))>sd.threshold,]
      }
    }
    
    if(transpose){
      .pcaPlotT(meth.mat,comp1=comp[1],comp2=comp[2],screeplot=screeplot, 
                adj.lim=adj.lim, 
                treatment=.Object@treatment,sample.ids=.Object@sample.ids,
                context=.Object@context
                ,scale=scale,center=center,obj.return=obj.return)
      
    }else{
      .pcaPlot(meth.mat,comp1=comp[1],comp2=comp[2],screeplot=screeplot, 
               adj.lim=adj.lim, 
               treatment=.Object@treatment,sample.ids=.Object@sample.ids,
               context=.Object@context,
               scale=scale,center=center,  obj.return=obj.return)
    }
    
  }      
)

#numC=grep("numCs", names(methidh))
#numT=grep("numTs", names(methidh))
#meth.mat=methidh[,numC]/(methidh[,numC]+methidh[,numT])



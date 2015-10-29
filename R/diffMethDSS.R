#' calculate Differential Methylation with DSS
#' 
#' This function provides an interface to the DSS package by Hao Wo at Emory University.
#' 
#' @param meth  a methylBase object
#' @param adjust different methods to correct the p-values for multiple testing. 
#'              Default is "SLIM" from methylKit. For "qvalue" please see \code{\link[qvalue]{qvalue}} 
#'              and for all other methods see \code{\link[stats]{p.adjust}}.
#' @param mc.cores integer denoting how many cores should be used for parallel
#'              differential methylation calculations (can only be used in
#'              machines with multiple cores).

#' @usage calculateDiffMethDSS(meth, adjust="SLIM", mc.cores=1)
#'                
#' @examples
#' 
#' data(methylKit)
#' 
#' dssDiffay <- calculateDiffMethDSS(methylBase.obj, adjust="SLIM", mc.cores=1)
#' 
#' @return a methylDiff object
#' @section Details:
#' test
#' 
#' @export
#' @docType methods
#' @rdname calculateDiffMethDSS-methods



setGeneric("calculateDiffMethDSS", function(meth, adjust=c("SLIM","holm","hochberg","hommel","bonferroni","BH","BY","fdr","none","qvalue"),
                                            mc.cores=1) standardGeneric("calculateDiffMethDSS"))

setMethod("calculateDiffMethDSS", "methylBase", 
          function(meth, adjust=c("SLIM","holm","hochberg","hommel","bonferroni","BH","BY","fdr","none","qvalue"),
                   mc.cores=1){
            
            #require(DSS)
            #library(DSS)
            #library(methylKit)
            
            #cat('Test', '\n')
            #cat('Removing NA values from methylBase object...', '\n')
            #meth=select(meth, (1:dim(meth)[[1]])[ apply(meth, 1, function (x) !any(is.na(x)) ) ] )
            
            #cat('Sorting methylBade object to enable us to match strand data later ...', '\n')
            #meth=new("methylDiff",getData(meth)[ with(getData(meth), order(chr, start)), ]  ,sample.ids=meth@sample.ids,assembly=meth@assembly,context=meth@context,
            #treatment=meth@treatment,destranded=meth@destranded,resolution=meth@resolution)
  
            sample_indices=1:length(meth@sample.ids)
            #BS1=methylBaseToBSeq(meth, use_samples=sample_indices[meth@treatment==0])
            #BS2=methylBaseToBSeq(meth, use_samples=sample_indices[meth@treatment!=0])
            
            if(length(unique(meth@treatment)) > 2){
                    warning("altering 'treatment' condition.\ncalculateDiffMethDSS() works with two groups only ")
                    meth@treatment[meth@treatment!=0]=1
            }
            
            raw_output=callDML2(meth)
            
            #ids=paste(raw_output[ , 'chr'] , raw_output[, 'pos'], sep='.')
            #raw_output$id=ids
            #raw_output=raw_output[ with(raw_output, order(chr, pos) ), ]
            raw_output$strand=getData(meth)$strand
            
            output=raw_output[, c('chr', 'pos', 'pos', 'strand', 'pval', 'fdr', 'diff')]
            output$diff=100*(output$diff)
            colnames(output)=c('chr', 'start', 'end', 'strand', 'pvalue', 'qvalue', 'meth.diff')
            
            # adjust pvalues
            output$qvalue=p.adjusted(output$pvalue,method=adjust)
            
            #if(slim) {
            #  cat('Trying SLIM...', '\n')
            #  slimObj=SLIMfunc(output[ , 'pvalue'])
            #  output$qvalue=QValuesfun(output[ , 'pvalue'], slimObj$pi0_Est)
            #}
            
            obj=new("methylDiff",output,sample.ids=meth@sample.ids,assembly=meth@assembly,context=meth@context,
                    treatment=meth@treatment,destranded=meth@destranded,resolution=meth@resolution)
            obj
        }
)

# wrapper function for SLIM, p.adjust, qvalue-package
p.adjusted <- function(pvals,method=c("SLIM","holm","hochberg","hommel","bonferroni","BH","BY","fdr","none","qvalue"),
                       n=length(pvals),fdr.level=NULL,pfdr=FALSE,STA=.1,Divi=10,Pz=0.05,B=100,Bplot=FALSE){
  
  method <- match.arg(method)
  
  qvals=switch(method,
               # SLIM function
               SLIM={QValuesfun(pvals,
                                SLIMfunc(pvals,STA=STA,Divi=Divi,Pz=Pz,B=B,Bplot=Bplot)$pi0_Est)
               },
               # r-base/p-adjust functions
               holm={p.adjust(pvals,method=method,n)},
               hochberg={p.adjust(pvals,method=method,n)},
               hommel={p.adjust(pvals,method=method,n)},
               bonferroni={p.adjust(pvals,method=method,n)},
               BH={p.adjust(pvals,method=method,n)},
               BY={p.adjust(pvals,method=method,n)},
               fdr={p.adjust(pvals,method=method,n)},
               none={p.adjust(pvals,method=method,n)},
               # r-bioconductor/qvalue-package function
               qvalue={qvalue(pvals,fdr.level,pfdr)$qvalues}
  )
}

#methylBaseToBSeq<-function(meth, use_samples=NULL, use_sites=NULL, remove.na=TRUE ) {
#  
#  #require(methylKit)
#  #require(DSS)
#  require(bsseq)
#  
#  if (is.null(use_samples)) use_samples=1:length(meth@sample.ids)
#  
#  if (is.null(use_sites)) use_sites=1:(dim(meth)[[1]])
#
#  cur_data=getData(select(meth, use_sites))
#    
#  if (remove.na) cur_data=na.omit(cur_data)
#  
#  M=as.matrix(cur_data [ , as.character(   lapply( use_samples, function (x) paste0('numCs',x)  ) )]  )
#  Cov=as.matrix(cur_data [ , as.character(   lapply( use_samples, function (x) paste0('coverage',x)  ) )]  )
#  
#  pos=cur_data$start
#  
#  chr=cur_data$chr
#  
#  BS1<-BSseq( M=M, Cov=Cov, pos=pos, chr=chr)
#  BS1
#}


#####################################################################
## CODE FROM THIS POINT IS COPIED FROM DSS PACKAGE
## BY HAO WO AT EMORY UNIVERSITY
######################################################################

######################################################################
## a list of functions for DML/DMR detections from Bisulfite seq data
######################################################################
#require(bsseq)

######################################
## wrapper function to calling DML. 
## This is without smoothing.
######################################
callDML <- function(BS1, BS2, equal.disp=FALSE, threshold=0) {
  ## first merge data from two conditions according to chr and pos
  cat('Using internal DSS code...', '\n')
  n1 <- getBSseq(BS1,"Cov"); x1 <- getBSseq(BS1,"M")
  n2 <-getBSseq(BS2,"Cov"); x2 <- getBSseq(BS2,"M")
  gr1 <- getBSseq(BS1,"gr");  gr2 <- getBSseq(BS2,"gr")
  alldata <- mergeData.counts(n1, x1, gr1, n2, x2, gr2)
  
  
  ## estimate priors from counts.
  if(equal.disp) {
    ## assume equal dispersion in two conditions. Combine counts from two conditions and estimate dispersions.
    ## Should keep only those sites didn't show much differences??
    ## I'll ignore that part for now. 
    ix.X <- grep("X", colnames(alldata))
    x1 <- alldata[,ix.X]
    ix.N <- grep("N", colnames(alldata))
    n1 <- alldata[,ix.N]
    prior1 <- est.prior.logN(x1, n1)
    prior2 <- 0
  } else { # different prior for two conditions
    n1 <- getBSseq(BS1,"Cov"); x1 <- getBSseq(BS1,"M")
    prior1 <- est.prior.logN(x1, n1)
    n2 <- getBSseq(BS2,"Cov"); x2 <- getBSseq(BS2,"M")
    prior2 <- est.prior.logN(x2, n2)
  }
  
  ## grab data 
  cc <- colnames(alldata)
  ix.X1 <- grep("X.*cond1", cc);
  ix.N1 <- grep("N.*cond1", cc)
  ix.X2 <- grep("X.*cond2", cc);
  ix.N2 <- grep("N.*cond2", cc)
  ncol1 <- length(ix.X1); ncol2 <- length(ix.X2)
  x1 <- as.matrix(alldata[,ix.X1]); n1 <- as.matrix(alldata[,ix.N1])
  x2 <- as.matrix(alldata[,ix.X2]); n2 <- as.matrix(alldata[,ix.N2])
  
  ## compute means. Spatial correlations are ignored at this time
  ## the means need to be of the same dimension as X and N
  estprob1 <- compute.mean(x1, n1);   estprob2 <- compute.mean(x2, n2)
  
  ## perform Wald test 
  wald <- waldTest.DML(x1, n1, estprob1, x2, n2, estprob2,
                       prior1, prior2, threshold, equal.disp=equal.disp)
  
  ## combine with chr/pos and output
  result <- data.frame(chr=alldata$chr, pos=alldata$pos, wald)
  
  ## sort result according to chr and pos
  #ix <- sortPos(alldata$chr, alldata$pos)
  
  #return(result[ix,])
  return(result[ ,])
}


# @param mbase methylBase object
callDML2 <- function(mbase, equal.disp=FALSE, threshold=0) {
  ## first merge data from two conditions according to chr and pos
  cat('Using internal DSS code...', '\n')
  n1 = as.matrix(getData(mbase)[mbase@coverage.index[mbase@treatment==0]])
  x1 = as.matrix(getData(mbase)[mbase@numCs.index[mbase@treatment==0]])
  
  n2 = as.matrix(getData(mbase)[mbase@coverage.index[mbase@treatment==1]])
  x2 = as.matrix(getData(mbase)[mbase@numCs.index[mbase@treatment==1]])

  gr1 <-as(mbase,"GRanges")[,0] ;  gr2 <- gr1
  
  # this could be improved, we can avoid the whole merge step
  # as mbase is already merged, and takes time on big data sets
  # we can just use cbind() 
  alldata <- mergeData.counts(n1, x1, gr1, n2, x2, gr2)
  
  
  ## estimate priors from counts.
  if(equal.disp) {
    ## assume equal dispersion in two conditions. Combine counts from two conditions and estimate dispersions.
    ## Should keep only those sites didn't show much differences??
    ## I'll ignore that part for now. 
    ix.X <- grep("X", colnames(alldata))
    x1 <- alldata[,ix.X]
    ix.N <- grep("N", colnames(alldata))
    n1 <- alldata[,ix.N]
    prior1 <- est.prior.logN(x1, n1)
    prior2 <- 0
  } else { # different prior for two conditions
    n1 = as.matrix(getData(mbase)[mbase@coverage.index[mbase@treatment==0]])
    x1 = as.matrix(getData(mbase)[mbase@numCs.index[mbase@treatment==0]])
    prior1 <- est.prior.logN(x1, n1)
   
    n2 = as.matrix(getData(mbase)[mbase@coverage.index[mbase@treatment==1]])
    x2 = as.matrix(getData(mbase)[mbase@numCs.index[mbase@treatment==1]])
    prior2 <- est.prior.logN(x2, n2)
  }
  
  ## grab data 
  cc <- colnames(alldata)
  ix.X1 <- grep("X.*cond1", cc);
  ix.N1 <- grep("N.*cond1", cc)
  ix.X2 <- grep("X.*cond2", cc);
  ix.N2 <- grep("N.*cond2", cc)
  ncol1 <- length(ix.X1); ncol2 <- length(ix.X2)
  x1 <- as.matrix(alldata[,ix.X1]); n1 <- as.matrix(alldata[,ix.N1])
  x2 <- as.matrix(alldata[,ix.X2]); n2 <- as.matrix(alldata[,ix.N2])
  
  ## compute means. Spatial correlations are ignored at this time
  ## the means need to be of the same dimension as X and N
  estprob1 <- compute.mean(x1, n1);   estprob2 <- compute.mean(x2, n2)
  
  ## perform Wald test 
  wald <- waldTest.DML(x1, n1, estprob1, x2, n2, estprob2,
                       prior1, prior2, threshold, equal.disp=equal.disp)
  
  ## combine with chr/pos and output
  result <- data.frame(chr=alldata$chr, pos=alldata$pos, wald)
  
  ## sort result according to chr and pos
  #ix <- sortPos(alldata$chr, alldata$pos)
  
  #return(result[ix,])
  return(result)
}



###############################################################################
## Perform Wald tests for calling DML.
###############################################################################
waldTest.DML <- function(x1,n1,estprob1, x2,n2, estprob2, prior1, prior2,
                         threshold, equal.disp=equal.disp) {
  
  if(equal.disp)
    prior2 <- prior1
  
  ## estimated shrunk dispersions
  ## - this part is slow. Should be computed parallely
  if(equal.disp) { ## equal dispersion. Combine two groups and shrink
    x <- cbind(x1, x2); n <- cbind(n1, n2)
    #    estprob <- cbind(matrix(rep(estprob1,ncol(x1)), ncol=ncol(x1)),
    #                     matrix(rep(estprob1,ncol(x2)), ncol=ncol(x2)))
    estprob <- cbind(estprob1, estprob2)
    shrk.phi1 <- shrk.phi2 <- dispersion.shrinkage(x, n, prior1, estprob)
    
  } else { ## shrink two groups separately 
    shrk.phi1 <- dispersion.shrinkage(x1, n1, prior1, estprob1)
    shrk.phi2 <- dispersion.shrinkage(x2, n2, prior2, estprob2)
  }
  
  ## Wald test
  wald <- compute.waldStat(estprob1[,1], estprob2[,1], n1, n2, shrk.phi1, shrk.phi2)
  
  ## obtain posterior probability that the differnce of two means are greater than a threshold
  if( threshold>0 ) {dat
    p1 <- pnorm(wald$diff-threshold, sd=wald$diff.se) ## Pr(delta.mu > threshold)
    p2 <- pnorm(wald$diff+threshold, sd=wald$diff.se, lower.tail=FALSE) ## Pr(-delta.mu < -threshold)
    pp.diff <- p1 + p2
    wald <- data.frame(wald, postprob.overThreshold=pp.diff)
  }
  
  return(wald)
}

###############################################
## compute Wald test statistics
## Need to deal with missing data
## Currently work for two conditions.
## Condition means are assumed to be the same among replicates.
###############################################
compute.waldStat <- function(estprob1, estprob2, n1, n2, phi1, phi2) {
  dif <- estprob1 - estprob2
  n1m <- rowSums(n1);    n2m <- rowSums(n2)
  var1 <- rowSums(n1*estprob1*(1-estprob1)*(1+(n1-1)*phi1)) / (n1m)^2
  var2 <- rowSums(n2*estprob2*(1-estprob2)*(1+(n2-1)*phi2)) / (n2m)^2
  ##vv <- var1/ncol1+var2/ncol2
  vv <- var1 + var2
  ## bound vv a little bit??
  vv[vv<1e-5] <- 1e-5
  se <- sqrt(vv)
  stat <- dif/se
  pval <- 2 * (1 - pnorm(abs(stat))) ## p-value for hypothesis testing
  fdr <- p.adjust(pval, method="fdr")
  
  data.frame(mu1=estprob1, mu2=estprob2, diff=dif, diff.se=se, stat=stat, pval=pval, fdr=fdr)
}

#######################################################
## some utility functions for BS-seq data
#######################################################

###### make an object of BSseq given count data from several replicates
## The input is a list of  data frames with columns: chr, pos, Ntotal, Nmethy.
makeBSseqData <- function(dat, sampleNames) {
  n0 <- length(dat)
  
  if(missing(sampleNames))
    sampleNames <- paste("sample", 1:n0, sep="")
  
  alldat <- dat[[1]]
  colnames(alldat)[3:4] <- c("N1", "X1")
  
  
  if(n0 > 1) { ## multiple replicates, merge data
    for(i in 2:n0) {
      thisdat <- dat[[i]]
      colnames(thisdat)[3:4] <- paste(colnames(thisdat)[3:4], i, sep="")
      alldat <- merge(alldat, thisdat, all=TRUE)
    }
  }
  
  ## make BSseq object
  ix.X <- grep("X", colnames(alldat))
  ix.N <- grep("N", colnames(alldat))
  alldat[is.na(alldat)] <- 0
  result <- BSseq(chr=alldat$chr, pos=alldat$pos, M=as.matrix(alldat[,ix.X, drop=FALSE]),
                  Cov=as.matrix(alldat[,ix.N, drop=FALSE]))
  
  result
}


########################################################################
## A function to estimate prior parameters, assume log-normal prior.
## It takes X and N, and only use the sites with big coverages,
## then return the mean and sd of prior distribution.
## Potentially, this should take mean if smoothing is allowed.
########################################################################
est.prior.logN <- function(X, N) {
  ## keep sites with large coverage and no missing data
  ix=rowMeans(N>10)==1 & rowSums(N==0)==0
  X=X[ix,]; N=N[ix,]
  ## compute sample mean/var
  p=X/N
  mm=rowMeans(p)
  mm[mm==0]=1e-5
  mm[mm==1]=1-1e-5
  vv=rowVars(p)
  phi=vv/mm/(1-mm)
  ## exclude those with vv==0. Those are sites with unobservable phis.
  ## But this will over estimate the prior.
  ## What will be the consequences????
  phi=phi[vv>0]
  lphi=log(phi[phi>0])
  prior.mean=median(lphi, na.rm=TRUE)
  prior.sd=IQR(lphi, na.rm=TRUE) /1.39 
  
  ## It seems this over-estimates the truth. Need to use the tricks in
  ## my biostat paper to remove the over-estimation. To be done later.
  c(prior.mean, prior.sd)
}



########################################
## Merge two sets of counts
########################################
mergeData.counts <- function(n1, x1, gr1, n2, x2, gr2) {
  allchr1 <- seqnames(gr1); allpos1 <- start(gr1)
  allchr2 <- seqnames(gr2); allpos2 <- start(gr2)
  
  colnames(x1) <- paste(paste("X", 1:ncol(x1), sep=""), "cond1", sep=".")
  colnames(n1) <- paste(paste("N", 1:ncol(n1), sep=""), "cond1", sep=".")
  colnames(x2) <- paste(paste("X", 1:ncol(x2), sep=""), "cond2", sep=".")
  colnames(n2) <- paste(paste("N", 1:ncol(n2), sep=""), "cond2", sep=".")
  
  
  
  dat1 <- data.frame(chr=as.character(seqnames(gr1)), pos=start(gr1), x1, n1)
  dat2 <- data.frame(chr=as.character(seqnames(gr2)), pos=start(gr2), x2, n2)
  alldat <- merge(dat1, dat2, all=FALSE) ## only keep the ones with data in both samples
  alldat[is.na(alldat)] <- 0
  
  return(alldat)
}

########################################
## the rowVars function
########################################
rowVars <- function (x, center = NULL, ...) {
  n <- !is.na(x)
  n <- rowSums(n)
  n[n <= 1] <- NA
  if (is.null(center)) {
    center <- rowMeans(x, ...)
  }
  x <- x - center
  x <- x * x
  x <- rowSums(x, ...)
  x <- x/(n - 1)
  x
}


######################################################################################
## function to compute means. If it's based on each CG site, just take an average.
## If include smoothing, use BSmoooth.
## This function is to be expanded.
######################################################################################
compute.mean <- function(X, N) {
  ## shrink the mean estimates a little bit.
  p <- (rowSums(X)+0.5)/(rowSums(N)+1)
  res <- matrix(rep(p, ncol(X)), ncol = ncol(X))
  return(res)
}


########################################################################
## Dispersion shrinkage based on log-normal penalized likelihood.
## Takes X, N, estimated mean and prior.
##
## The shrinakge is done in log scale. So data will be shrink to the
## logarithmic means.
########################################################################
dispersion.shrinkage <- function(X, N, prior, estprob) {
  ## penalized likelihood function
  plik.logN <- function(size, X,mu,m0,tau,phi) {
    -(sum(dbb(size, X, mu, exp(phi))) + dnorm(phi, mean=m0, sd=tau, log=TRUE))
  }
  ## for CG sites with no coverage, use prior 
  shrk.phi=exp(rep(prior[1],nrow(N)))
  
  ## deal with estprob, make it a matrix if not.
  if(!is.matrix(estprob))
    estprob <- as.matrix(estprob)
  
  ## skip those without coverage, or no replicates.
  ix <- rowSums(N>0) > 0
  X2 <- X[ix, ,drop=FALSE]; N2 <- N[ix,,drop=FALSE]; estprob2 <- estprob[ix,,drop=FALSE]
  shrk.phi2 <- rep(0, nrow(X2))
  for(i in 1:nrow(X2)) {
    ## I can keep the 0's with calculation. They don't make any difference.
    shrk.one=optimize(f=plik.logN, size=N2[i,], X=X2[i,], mu=estprob2[i,], m0=prior[1], tau=prior[2],
                      interval=c(-5, log(0.99)),tol=1e-4)
    shrk.phi2[i]=exp(shrk.one$minimum)
  }
  shrk.phi[ix] <- shrk.phi2
  
  return(shrk.phi)
}

#########################################################
## beta-binomial (BB) density function.
## The BB distribution is parametrized by mean and dispersion.
#########################################################
dbb <- function (size, x, mu, phi, log=TRUE)  {
  ## first convert mu/phi to alpha/beta
  tmp=1/phi-1
  alpha=mu*tmp
  beta=tmp - alpha
  v=lchoose(size,x)-lbeta(beta, alpha)+lbeta(size-x + beta,x+alpha)
  if(!log)
    return(exp(v))
  else return(v)
}

###########################################
## sort according to chr and pos.
## Return the index
###########################################
sortPos <- function(chr, pos) {
  do.call(order, list(chr, pos))
}
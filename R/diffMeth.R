## PART THAT DEALS WITH DIFFERENTIAL METHYLATION CALCULATIONS

##############################################################################
## S3 functions to be used in S4 stuff
##############################################################################

# SLIM
######################################
#####SLIM pi0 Estimation function
#####Copyright by Tsai Lab of UGA, US, and Hong-Qiang Wang, IIM, CAS, China
#####Reference: SLIM: A Sliding Linear Model for Estimating the Proportion of True Null Hypotheses in Datasets With Dependence Structures
#####usage:
#####inputs: 
#####rawp:p-values, required
#####STA:lambda1
#####Divi: the number of segments
#####Pz: maximum of p-values of alternative tests
#####B: the number of quantile points
#####Bplot: logic, if true, drawing lambda-gamma plot
#####outputs: 
#####pi0_Est: estimated value of pi0
#####selQuantile: the value of alpha determined 

#####################################


SLIMfunc<-function(rawp,STA=.1,Divi=10,Pz=0.05,B=100,Bplot=FALSE)
{


  ####################
  m <- length(rawp) 
 
  ########################
  alpha_mtx=NULL;#
  pi0s_est_COM=NULL;
  SzCluster_mtx=NULL;
  P_pi1_mtx=NULL;
  pi1_act_mtx=NULL;
  Dist_group=NULL;
  Num_mtx=NULL;
  Gamma=NULL;
  PI0=NULL;
  
  #############
  ##observed points
  lambda_ga=seq(0,1,0.001);
  gamma_ga=sapply(lambda_ga,f1,rawp=rawp);
  Gamma=c(Gamma,gamma_ga);
  alpha_mtx=c(alpha_mtx,gamma_ga[which(lambda_ga==0.05)]);
  
  
  ###esimation
  pi0_mtx=NULL;
  x.axis=NULL;
  y.axis=NULL;
  itv=(1-STA)/Divi;
  for (i in 1:Divi)##10 for uniform data
  {
    cutoff=STA+(i/Divi)*(1-STA);##10 for uniform data
    lambda=seq(cutoff-itv,cutoff,itv/10);
    gamma_mtx=sapply(lambda,f1,rawp=rawp);
    LModel=lm(gamma_mtx~lambda);
    pi0_mtx=c(pi0_mtx,coefficients(LModel)[2]);
  }


  ##################################
  ########searching
  N_COM=NULL;
  N_rawp=NULL;
  maxFDR_mtx=NULL;
  quapoint_mtx=NULL; 
  if (B<=1) B=100;
  quapoint_mtx=seq(0.01,0.99,1/B);
  for (k in 1:length(quapoint_mtx))
  {
    qua_point=quapoint_mtx[k];
    
    pi0_combLR=min(quantile(pi0_mtx,qua_point),1);#mean(pi0_mtx);#median();# qua_point=0.78 for desreasing distribution;
    ##0.4 for uniform or normal distribution;
    pi0_est=pi0_combLR;
    
    ###########Calculate independent index of raw p vlaues
    PI0=rbind(PI0,pi0_mtx);
    
    pi0s_est_COM=c(pi0s_est_COM,pi0_est);
    ##Condition1
    P_pi1=sort(rawp)[max(length(rawp)*(1-pi0_est),1)];##
    P_pi1_mtx=c(P_pi1_mtx,P_pi1);
    
    pi0=pi0_est;
    if (is.null(Pz)) Pz=0.05;
    maxFDR=Pz*pi0/(1-(1-Pz)*pi0);
    maxFDR_mtx=c(maxFDR_mtx,maxFDR);
    
    qvalues_combLR=QValuesfun(rawp,pi0);
    qvalue_cf=maxFDR;
    selected=which(qvalues_combLR<qvalue_cf);
    Sel_qvalues_combLR=selected;
    
    pi1_act_mtx=c(pi1_act_mtx,length(Sel_qvalues_combLR)/length(rawp));
    N_COM=c(N_COM,list(Sel_qvalues_combLR));
    Num_mtx=c(Num_mtx,length(Sel_qvalues_combLR));
  }
  length(N_COM)
  length(quapoint_mtx)
  
  ####doing judging
  ##by max FDR
  pi1s_est_COM=1-pi0s_est_COM;
  Diff=sum(rawp<=Pz)/length(rawp)-pi1_act_mtx;
  
  ###
  loc=which.min(abs(Diff));
  Diff.loc=Diff[loc];
  selQuantile=quapoint_mtx[loc];
  pi0_Est=min(1,pi0s_est_COM[loc]);
  maxFDR.Pz=Pz*pi0_Est/(1-(1-Pz)*pi0_Est);
  
  if(Bplot)
  {
    #windows();
    par(mfrow=c(1,2));
    hist(rawp,main="Histogram of p-value");
    gamma_ga=sapply(lambda_ga,f1,rawp=rawp);
    plot(lambda_ga,gamma_ga,type="l",main="Relationship of p- and q value",xlab=expression(lambda),ylab=expression(gamma),cex.lab=1.45,cex.axis=1.42)
    #par(xaxp=c(0,1,10));
    #axis(1);
    ##qvalues
    qValues=QValuesfun(rawp,pi0=pi0_Est);
    gammaq_ga=sapply(lambda_ga,f1,rawp=qValues);
    lines(lambda_ga,gammaq_ga,col="blue",lwd=2,lty="dashed")
    abline(v=Pz,col="black",lwd=2,lty="dotdash")
    abline(v=maxFDR.Pz,col="blue",lwd=2,lty="dotdash")
    text(0.75,0.6,labels=paste("L=",round(abs(Diff.loc),4),sep=""));
    leg=list(bquote("CPD of p-value"),bquote("CPD of q-value"),bquote("Pmax"==.(Pz)),bquote("FDRmax"==.(round(maxFDR.Pz,2))));
    legend("bottomright",legend=as.expression(leg),lwd=2,lty=c("solid","dashed","dotdash","dotdash"),col=c("black","blue","black","blue"));
  }
  
  return(list(pi0_Est=pi0_Est,selQuantile=selQuantile));
}

####
f1<-function(cutoff,rawp){sum(rawp<cutoff)/length(rawp)};

###
QValuesfun<-function(rawp,pi0)
{
  order_rawp=sort(rawp);
  #order_rawp=sort(rawp,method="qu");
  qvalues=pi0*length(order_rawp)*order_rawp/c(1:length(order_rawp));
  temp=cummin(qvalues[seq(length(qvalues),1,-1)])
  qvalues=temp[seq(length(temp),1,-1)];
  qvalues=qvalues[order(order(rawp))]
}

glm.set<-function(set,numC1.ind,numC2.ind,numT1.ind,numT2.ind)
{

  Treat <-  c( rep(1,length(numC1.ind)),rep(0,length(numC2.ind)) ) # get the treatment vector

  Cs=as.matrix(set[,c(numC1.ind,numC2.ind)])
  Ts=as.matrix(set[,c(numT1.ind,numT2.ind)])

  sums=Cs+Ts # get number of Cs and Ts
  Ps  =Cs/sums # get probability of Cs
  
  Tmod=model.matrix(~Treat)
  glm.bare<-function(X,Y) # X probs > Ps, Y total number > sums
  {
    #cat(Tmod)
    obj=glm.fit(Tmod,X,weights=Y,family=binomial(link=logit))
    deviance <- obj$null.deviance - obj$deviance
    dispersion=1 #(if binomial or poisson)
    aliased <- is.na(coef(obj))
    p <- obj$rank
    if (p > 0) { # if clause and the rest to get the t-value or wald statistic
        p1 <- 1L:p
        Qr <- obj$qr
        coef.p <- obj$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- (coef.p/s.err)[2]    
     }else if (obj$df.residual<=0)
     { tvalue=NaN }
 
    wald <- tvalue
    beta1 <- obj$coefficients[2]
    p.value <- 1-pchisq(deviance,df=1)

    return( c(wald,beta1,p.value) )
  }

  Ps.list=split(Ps,1:nrow(Ps) ) # get Probs
  Ws.list=split(sums,1:nrow(sums) ) # get weights to feed into function
  res=t(mapply(glm.bare,Ps.list,Ws.list))
  colnames(res)=c("wald","beta1","pvalue")
  return(res)
}


glm.set.mc<-function(set,numC1.ind,numC2.ind,numT1.ind,numT2.ind,n.mc)
{

  Treat <-  c( rep(1,length(numC1.ind)),rep(0,length(numC2.ind)) ) # get the treatment vector

  Cs=as.matrix(set[,c(numC1.ind,numC2.ind)])
  Ts=as.matrix(set[,c(numT1.ind,numT2.ind)])

  sums=Cs+Ts # get number of Cs and Ts
  Ps  =Cs/sums # get probability of Cs
  
  Tmod=model.matrix(~Treat)
  glm.bare<-function(dat,indX,indY,Tmod) # X probs > Ps, Y total number > sums
  {
    #cat(Tmod)
    obj=glm.fit(Tmod,dat[indX],weights=dat[indY],family=binomial(link=logit))
    deviance <- obj$null.deviance - obj$deviance
    dispersion=1 #(if binomial or poisson)
    aliased <- is.na(coef(obj))
    p <- obj$rank
    if (p > 0) { # if clause and the rest to get the t-value or wald statistic
        p1 <- 1L:p
        Qr <- obj$qr
        coef.p <- obj$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- (coef.p/s.err)[2]    
     }else if (obj$df.residual<=0)
     { tvalue=NaN }
 
    wald <- tvalue
    beta1 <- obj$coefficients[2]
    p.value <- 1-pchisq(deviance,df=1)

    return( c(wald,beta1,p.value) )
  }

  PWs=cbind(Ps,sums)
  PWs.list=split(PWs,1:nrow(PWs) ) # get weights to feed into function
  res=simplify2array(mclapply(PWs.list,glm.bare,indX=1:ncol(Ps),indY=ncol(Ps)+(1:ncol(sums)),Tmod=Tmod,  mc.cores=n.mc))
  res=data.frame(t(res))
  colnames(res)=c("wald","beta1","pvalue")
  return(res)
}

### GLMs that can deal with NA values

glm.set.v1<-function(set,numC1.ind,numC2.ind,numT1.ind,numT2.ind)
{

  Treat <-  c( rep(1,length(numC1.ind)),rep(0,length(numC2.ind)) ) # get the treatment vector

  Cs=as.matrix(set[,c(numC1.ind,numC2.ind)])
  Ts=as.matrix(set[,c(numT1.ind,numT2.ind)])

  sums=Cs+Ts # get number of Cs and Ts
  Ps  =Cs/sums # get probability of Cs
  
  #Tmod=model.matrix(~Treat)
  #X=c(0.1,NA,0.6,0.8)
  #Y=c(30,NA,40,100)
  glm.bare<-function(X,Y,Treat) # X probs > Ps, Y total number > sums
  {
    #cat(Tmod)
    Treat2=Treat[!is.na(X)]
    Tmod=model.matrix(~Treat2)
    obj=glm.fit(Tmod,X[!is.na(X)],weights=Y[! is.na(Y)],family=binomial(link=logit))
    deviance <- obj$null.deviance - obj$deviance
    dispersion=1 #(if binomial or poisson)
    aliased <- is.na(coef(obj))
    p <- obj$rank
    if (p > 0) { # if clause and the rest to get the t-value or wald statistic
        p1 <- 1L:p
        Qr <- obj$qr
        coef.p <- obj$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- (coef.p/s.err)[2]    
     }else if (obj$df.residual<=0)
     { tvalue=NaN }
 
    wald <- tvalue
    beta1 <- obj$coefficients[2]
    p.value <- 1-pchisq(deviance,df=1)

    return( c(wald,beta1,p.value) )
  }

  Ps.list=split(Ps,1:nrow(Ps) ) # get Probs
  Ws.list=split(sums,1:nrow(sums) ) # get weights to feed into function
  res=t(mapply(glm.bare,Ps.list,Ws.list,MoreArgs = list(Treat = Treat)))
  colnames(res)=c("wald","beta1","pvalue")
  return(res)
}



glm.set.mc.v1<-function(set,numC1.ind,numC2.ind,numT1.ind,numT2.ind,n.mc)
{

  Treat <-  c( rep(1,length(numC1.ind)),rep(0,length(numC2.ind)) ) # get the treatment vector

  Cs=as.matrix(set[,c(numC1.ind,numC2.ind)])
  Ts=as.matrix(set[,c(numT1.ind,numT2.ind)])

  sums=Cs+Ts # get number of Cs and Ts
  Ps  =Cs/sums # get probability of Cs
  
  #Tmod=model.matrix(~Treat)
  glm.bare<-function(dat,indX,indY,Treat) # X probs > Ps, Y total number > sums
  {
    #cat(Tmod)

    X=dat[indX]
    Y=dat[indY]
    Treat2=Treat[!is.na(X)]
    Tmod=model.matrix(~Treat2)
    obj=glm.fit(Tmod,X[! is.na(X)],weights=Y[! is.na(Y)],family=binomial(link=logit))
    deviance <- obj$null.deviance - obj$deviance
    dispersion=1 #(if binomial or poisson)
    aliased <- is.na(coef(obj))
    p <- obj$rank
    if (p > 0) { # if clause and the rest to get the t-value or wald statistic
        p1 <- 1L:p
        Qr <- obj$qr
        coef.p <- obj$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- (coef.p/s.err)[2]    
     }else if (obj$df.residual<=0)
     { tvalue=NaN }
 
    wald <- tvalue
    beta1 <- obj$coefficients[2]
    p.value <- 1-pchisq(deviance,df=1)

    return( c(wald,beta1,p.value) )
  }

  PWs=cbind(Ps,sums)
  PWs.list=split(PWs,1:nrow(PWs) ) # get weights to feed into function
  res=simplify2array(mclapply(PWs.list,glm.bare,indX=1:ncol(Ps),indY=ncol(Ps)+(1:ncol(sums)),Treat=Treat,  mc.cores=n.mc))
  res=data.frame(t(res))
  colnames(res)=c("wald","beta1","pvalue")
  return(res)
}


### END OF GLMs that can deal with NA values


# function to fix logistic regression pvalues and make qvalues
fix.q.values.glm<-function(pvals,slim=FALSE)
{
   if(slim==FALSE){
    #qvals=qvalue::qvalue(pvals[,3])$qvalues # get qvalues
    qvals=p.adjust(pvals,"BH")
  }else{
    slimObj=SLIMfunc(pvals[,3]);qvals=QValuesfun(pvals[,3], slimObj$pi0_Est)
  }                

  pvals=cbind(pvals,qvalue=qvals) # merge pvals and qvals
  return(pvals)
}

# function to fix fisher.test pvalues and make qvalues
fix.q.values.fisher<-function(pvals,slim=FALSE)
{
  if(slim==FALSE){
    #qvals=qvalue::qvalue(pvals)$qvalues # get qvalues
    qvals=p.adjust(pvals,"BH")
    
  }else{
  slimObj=SLIMfunc(pvals);qvals=QValuesfun(pvals, slimObj$pi0_Est)
  }                

  pvals=data.frame(pvalue=pvals,qvalue=qvals) # merge pvals and qvals
  return(pvals)
}
 
# A FASTER VERSION OF FISHERs EXACT
fast.fisher<-function (x, y = NULL, workspace = 2e+05, hybrid = FALSE, control = list(), 
    or = 1, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95, 
    simulate.p.value = FALSE, B = 2000, cache=F) 
{
    if (nrow(x)!=2 | ncol(x)!=2) stop("Incorrect input format for fast.fisher")
    #if (cache) {
    #  key = paste(x,collapse="_")
    # cachedResult = hashTable[[key]]
    #  if (!is.null(cachedResult)) {
    #    return(cachedResult)
    #  }
    #}
    # ---- START: cut version of fisher.test ----
    DNAME <- deparse(substitute(x))
    METHOD <- "Fisher's Exact Test for Count Data"
    nr <- nrow(x)
    nc <- ncol(x)
    PVAL <- NULL
    if ((nr == 2) && (nc == 2)) {
        m <- sum(x[, 1])
        n <- sum(x[, 2])
        k <- sum(x[1, ])
        x <- x[1, 1]
        lo <- max(0, k - n)
        hi <- min(k, m)
        NVAL <- or
        names(NVAL) <- "odds ratio"
        support <- lo:hi
        logdc <- dhyper(support, m, n, k, log = TRUE)
        dnhyper <- function(ncp) {
            d <- logdc + log(ncp) * support
            d <- exp(d - max(d))
            d/sum(d)
        }
        mnhyper <- function(ncp) {
            if (ncp == 0) 
                return(lo)
            if (ncp == Inf) 
                return(hi)
            sum(support * dnhyper(ncp))
        }
        pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
            if (ncp == 1) {
                if (upper.tail) 
                  return(phyper(x - 1, m, n, k, lower.tail = FALSE))
                else return(phyper(x, m, n, k))
            }
            if (ncp == 0) {
                if (upper.tail) 
                  return(as.numeric(q <= lo))
                else return(as.numeric(q >= lo))
            }
            if (ncp == Inf) {
                if (upper.tail) 
                  return(as.numeric(q <= hi))
                else return(as.numeric(q >= hi))
            }
            d <- dnhyper(ncp)
            if (upper.tail) 
                sum(d[support >= q])
            else sum(d[support <= q])
        }
        if (is.null(PVAL)) {
            PVAL <- switch(alternative, less = pnhyper(x, or), 
                greater = pnhyper(x, or, upper.tail = TRUE), 
                two.sided = {
                  if (or == 0) 
                    as.numeric(x == lo)
                  else if (or == Inf) 
                    as.numeric(x == hi)
                  else {
                    relErr <- 1 + 10^(-7)
                    d <- dnhyper(or)
                    sum(d[d <= d[x - lo + 1] * relErr])
                  }
                })
            RVAL <- list(p.value = PVAL)
        }
        mle <- function(x) {
            if (x == lo) 
                return(0)
            if (x == hi) 
                return(Inf)
            mu <- mnhyper(1)
            if (mu > x) 
                uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
            else if (mu < x) 
                1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 
                  1))$root
            else 1
        }
        ESTIMATE <- mle(x)
        #names(ESTIMATE) <- "odds ratio"
        RVAL <- c(RVAL, estimate = ESTIMATE, null.value = NVAL)
    }
    RVAL <- c(RVAL, alternative = alternative, method = METHOD, data.name = DNAME)
    attr(RVAL, "class") <- "htest"
    # ---- END: cut version of fisher.test ----    
    #if (cache) hashTable[[key]] <<- RVAL # write to global variable
    return(RVAL)                                                                         
}

mc.fish<-function(my.list,num.cores)
{

unlist( parallel::mclapply( my.list,function(x) fast.fisher(matrix(as.numeric( x) ,ncol=2,byrow=T),conf.int = F)$p.value,
                                                         mc.cores=num.cores,mc.preschedule = TRUE) ) 
}

# end of S3 functions

##############################################################################
## S4 OBJECTS
##############################################################################


#' An S4 class that holds differential methylation information
#'
#' This class is designed to hold statistics and locations for differentially methylated regions/bases. It extends \code{\link{data.frame}} class.
#'  \code{\link[methylKit]{calculateDiffMeth}} function returns an object with \code{methylDiff} class.
#'          
#' @section Slots:\describe{
#'    \item{\code{sample.ids}}{ids/names of samples in a vector}
#'    \item{\code{assembly}}{a name of genome assembly, such as :hg18,mm9, etc}
#'    \item{\code{context}}{numeric vector identifying which samples are which
#'    group }
#'    \item{\code{treatment}}{numeric vector identifying which samples are which
#'     group }
#'    \item{\code{destranded}}{logical denoting if methylation inormation is
#'     destranded or not}
#'    \item{\code{resolution}}{string either 'base' or 'region' defining the 
#'    resolution of methylation information}
#'    \item{\code{.Data}}{data.frame holding the locations and statistics}
#'
#' }
#' 
#' @section Details:
#' \code{methylDiff} class extends \code{\link{data.frame}} class therefore
#'  providing novice and experienced R users with a data structure that is 
#'  well known and ubiquitous in many R packages.
#' 
#' @section Subsetting:
#'  In the following code snippets, \code{x} is a \code{methylDiff}.
#'  Subsetting by \code{x[i,]} will produce a new object if subsetting is done on
#'  rows. Column subsetting is not directly allowed to prevent errors in the 
#'  downstream analysis. see ?methylKit[ .
#' 
#' @section Coercion:
#'   \code{methylDiff} object can be coerced to \code{\link[GenomicRanges]{GRanges}} object via \code{\link{as}} function.
#' 
#' @section Accessors: 
#' The following functions provides access to data slots of methylDiff:
#' \code{\link[methylKit]{getData}},\code{\link[methylKit]{getAssembly}},\code{\link[methylKit]{getContext}}
#' 
#' @examples
#' data(methylKit)
#' library(GenomicRanges)
#' my.gr=as(methylDiff.obj,"GRanges")
#' 
#' @name methylDiff-class
#' @aliases methylDiff
#' @rdname methylDiff-class
#' @export
#' @docType class
setClass("methylDiff",representation(
  sample.ids = "character", assembly = "character",context = "character",treatment="numeric",destranded="logical",resolution="character"),contains="data.frame")


##############################################################################
## S4 FUNCTIONS
##############################################################################


#' Calculate differential methylation statistics
#' 
#' The function calculates differential methylation statistics between two groups 
#' of samples. The function uses either logistic regression test
#' or Fisher's Exact test to calculate differential methylation. 
#' See references for detailed explanation on statistics.
#' 
#' @param .Object a methylBase object to calculate differential methylation
#' @param slim If TRUE(default) SLIM method will be used for P-value adjustment.
#'             If FALSE, \code{\link{p.adjust}} with method="BH" option will be 
#'             used for P-value correction.
#' @param weigthed.mean calculate the mean methylation difference between groups 
#'                      using read coverage as weights
#' @param num.cores  integer for denoting how many cores should be used for 
#'                   differential methylation calculations (only can be used in
#'                    machines with multiple cores)
#' @usage calculateDiffMeth(.Object,slim=TRUE,weigthed.mean=TRUE,num.cores=1)
#' @examples
#' 
#' data(methylKit)
#' 
#' # Logistic regression test will be applied since there are multiple samples in each group
#' # in methylBase.obj object
#' my.diffMeth=calculateDiffMeth(methylBase.obj,slim=TRUE,weigthed.mean=TRUE,num.cores=1)
#' 
#' # pool samples in each group
#' pooled.methylBase=pool(methylBase.obj,sample.ids=c("test","control"))
#' 
#' # After applying pool() function, there is one sample in each group.
#' # Fisher's exact test will be applied for differential methylation
#' my.diffMeth2=calculateDiffMeth(pooled.methylBase,slim=TRUE,
#'                                weigthed.mean=TRUE,num.cores=1)
#' 
#' 
#' 
#' @return a methylDiff object containing the differential methylation 
#'                      statistics and locations
#' @section Details:
#'  The function either uses a logistic regression 
#'  (when there are multiple samples per group) or fisher's exact 
#'  when there is one sample per group.
#' @references Altuna Akalin, Matthias Kormaksson, Sheng Li,
#'             Francine E. Garrett-Bakelman, Maria E. Figueroa, Ari Melnick, 
#'             Christopher E. Mason. (2012). 
#'             "methylKit: A comprehensive R package for the analysis 
#'             of genome-wide DNA methylation profiles." Genome Biology. 
#' @seealso \code{\link[methylKit]{pool}}, \code{\link[methylKit]{reorganize}}
#' 
#' @export
#' @docType methods
#' @rdname calculateDiffMeth-methods
setGeneric("calculateDiffMeth", function(.Object,slim=TRUE,weigthed.mean=TRUE,
                              num.cores=1) standardGeneric("calculateDiffMeth"))

#' @aliases calculateDiffMeth,methylBase-method
#' @rdname calculateDiffMeth-methods
setMethod("calculateDiffMeth", "methylBase",
                    function(.Object,slim,weigthed.mean,num.cores){
    
    #get CpGs with the cutoff
    #inds=rowSums( S3Part(.Object)[,.Object@coverage.index]>=coverage.cutoff) == length(.Object@coverage.index,na.rm=TRUE)
    #subst=S3Part(.Object)[inds,]
    subst=S3Part(.Object,strictS3 = TRUE)
    
    if(length(.Object@treatment)<2 ){
      stop("can not do differential methylation calculation with less than two samples")
    }
    if(length(unique(.Object@treatment))<2 ){
      stop("can not do differential methylation calculation when there is no control\n
           treatment option should have 0 and 1 designating treatment and control samples")
    }
    
    if(length(unique(.Object@treatment))>2 ){
      stop("can not do differential methylation calculation when there are more than\n
           two groups, treatment vector indicates more than two groups")
    }
       
       
    # get the indices for numCs and numTs in each set
    set1.Cs=.Object@numCs.index[.Object@treatment==1]
    set2.Cs=.Object@numCs.index[.Object@treatment==0]
    set1.Ts=.Object@numTs.index[.Object@treatment==1]
    set2.Ts=.Object@numTs.index[.Object@treatment==0]
  
    # if one control, one treatment case to the fisher's exact test
    if(length(.Object@treatment)==2 )
    {
      if(num.cores>1){
        my.f.list=split(as.matrix(subst[,c(set1.Cs,set1.Ts,set2.Cs,set2.Ts)]),1:nrow(subst))
        f.fisher=function(x){ fast.fisher(matrix( x ,ncol=2,byrow=T),conf.int = F)$p.value }
        pvals    = unlist( parallel::mclapply( my.f.list ,f.fisher, mc.cores=num.cores) ) # apply fisher test
        #pvals    =mc.fish(my.f.list,num.cores)
  
      }else{
        pvals =apply( subst[,c(set1.Cs,set1.Ts,set2.Cs,set2.Ts)],1,function(x) fast.fisher(matrix(as.numeric(x),ncol=2,byrow=T),conf.int = F)$p.value ) # apply fisher test
      }
      pvals    = fix.q.values.fisher(pvals,slim=slim)   
      
      # calculate mean methylation change
      mom.meth1    = 100*(subst[,set1.Cs]/subst[,set1.Cs-1]) # get % methylation
      mom.meth2    = 100*(subst[,set2.Cs]/subst[,set2.Cs-1])
      mom.mean.diff=mom.meth1-mom.meth2 # get difference between percent methylations
      x=data.frame(subst[,1:4],pvals,meth.diff=mom.mean.diff,stringsAsFactors=F) # make a data frame and return it
      obj=new("methylDiff",x,sample.ids=.Object@sample.ids,assembly=.Object@assembly,context=.Object@context,
          treatment=.Object@treatment,destranded=.Object@destranded,resolution=.Object@resolution)
      obj
    }
    else # else do the GLM - logistic regression
    { 
                              
        # pvalues
        if(num.cores>1){
          #pvals  = glm.set.mc(subst,set1.Cs,set2.Cs,set1.Ts,set2.Ts,num.cores)
          pvals  = glm.set.mc.v1(subst,set1.Cs,set2.Cs,set1.Ts,set2.Ts,num.cores)
  
        }else{
          #pvals  = glm.set(subst,set1.Cs,set2.Cs,set1.Ts,set2.Ts) # get p-values
          pvals  = glm.set.v1(subst,set1.Cs,set2.Cs,set1.Ts,set2.Ts) # get p-values
        }
  
        # get qvalues
        pvals  = fix.q.values.glm(pvals,slim=slim)   
        
        # calculate mean methylation change
        if(length(set1.Cs) > 1){
          mom.meth1=100*rowMeans(subst[,set1.Cs]/subst[,set1.Cs-1],na.rm=TRUE) # get means of means
          pm.meth1=100*rowSums(subst[,set1.Cs],na.rm=TRUE)/rowSums(subst[,set1.Cs-1],na.rm=TRUE) # get weigthed means
        }else{
          mom.meth1    = 100*(subst[,set1.Cs]/subst[,set1.Cs-1]) # get % methylation
          pm.meth1     = mom.meth1
        }
  
        if(length(set2.Cs)>1){
          mom.meth2=100*rowMeans(subst[,set2.Cs]/subst[,set2.Cs-1],na.rm=TRUE)
          pm.meth2=100*rowSums(subst[,set2.Cs],na.rm=TRUE)/rowSums(subst[,set2.Cs-1],na.rm=TRUE) # get weigthed means
        }else{
          mom.meth2    = 100*(subst[,set2.Cs]/subst[,set2.Cs-1]) # get % methylation
          pm.meth2     = mom.meth2
        }
        pm.mean.diff=pm.meth1-pm.meth2
        mom.mean.diff=mom.meth1-mom.meth2
        
        if(weigthed.mean){
          x=data.frame(subst[,1:4],pvals[,3:4],meth.diff=pm.mean.diff,stringsAsFactors=F) # make a data frame and return it
          obj=new("methylDiff",x,sample.ids=.Object@sample.ids,assembly=.Object@assembly,context=.Object@context,
            treatment=.Object@treatment,destranded=.Object@destranded,resolution=.Object@resolution)
          obj
  
        }
        else{
          x=data.frame(subst[,1:4],pvals[,3:4],meth.diff=mom.mean.diff,stringsAsFactors=F) # make a data frame and return it
          obj=new("methylDiff",x,sample.ids=.Object@sample.ids,assembly=.Object@assembly,context=.Object@context,
            treatment=.Object@treatment,destranded=.Object@destranded,resolution=.Object@resolution)
          obj
        }
      
    }
  }    
)


# differential methylation summary
# outputs a summary.methylDiff object
# that shows:
# total number of DMCs
# total number of CpGs covered in given assays
# total number of DMCs per Chromosome
# total number CpGs covered per Chr
#setGeneric(name="synopsis", def=function(.Object,difference=25,qvalue=0.01) standardGeneric("synopsis"))
#setMethod(f="synopsis", signature="methylDiff", 
#          definition=function(.Object,difference,qvalue) {
#                    cat("hi")                         
#})

# a class that holds differential methylation information
# 
#setClass("synopsis.methylDiff",representation(
#qvalue="numeric",
#  difference="numeric",
#  sample.ids = "character", 
#  assembly = "character",
#  treatment="numeric",
#  destranded="logical"),contains="data.frame")



##############################################################################
## ACESSOR FUNCTIONS FOR methylDiff OBJECT
##############################################################################


#' @rdname show-methods
#' @aliases show,methylDiff
setMethod("show", "methylDiff", function(object) {
  
  cat("methylDiff object with",nrow(object),"rows\n--------------\n")
  print(head(object))
  cat("--------------\n")
  cat("sample.ids:",object@sample.ids,"\n")
  cat("destranded",object@destranded,"\n")
  cat("assembly:",object@assembly,"\n")
  cat("context:", object@context,"\n")
  cat("treament:", object@treatment,"\n")
  cat("resolution:", object@resolution,"\n")
})

#' @rdname getContext-methods
#' @aliases getContext,methylDiff-method
setMethod("getContext", signature="methylDiff", definition=function(x) {
                return(x@context)
        })


#' @rdname getAssembly-methods
#' @aliases getAssembly,methylDiff-method
setMethod("getAssembly", signature="methylDiff", definition=function(x) {
  return(x@assembly)
})


# a function for getting data part of methylDiff    
#' @rdname getData-methods
#' @aliases getData,methylDiff-method
setMethod(f="getData", signature="methylDiff", definition=function(x) {
                #return(as(x,"data.frame"))
                return(S3Part(x, strictS3 = TRUE))
        }) 



##############################################################################
## CONVERTOR FUNCTIONS FOR methylDiff OBJECT
##############################################################################

setAs("methylDiff", "GRanges", function(from)
{
  
  GRanges(seqnames=from$chr,ranges=IRanges(start=from$start, end=from$end),
          strand=from$strand, 
          qvalue=from$qvalue,
          meth.diff=from$meth.diff
  )
  
})



### subset methylDiff


#' @aliases select,methylDiff-method
#' @rdname select-methods
setMethod("select", "methylDiff",
          function(x, i)
          {
            
            new("methylDiff",getData(x)[i,],
                sample.ids = x@sample.ids,
                assembly = x@assembly,
                context = x@context,
                treatment=x@treatment,
                destranded=x@destranded,
                resolution=x@resolution)
          }
)


#' @rdname extract-methods
#' @aliases [,methylDiff-method
setMethod("[","methylDiff", 
          function(x,i,j){
            #cat(missing(i),"\n",missing(j),"\n",missing(drop))
            if(!missing(j)){
              stop(paste("subsetting on columns is not allowed for",class(x),
                         "object\nif you want to do extraction on the data part", 
                         "of the object use getData() first"),
                   call. = FALSE)
            }
            select(x,i)
          }
)





                      
#' get differentially methylated regions/bases based on cutoffs 
#' 
#' The function subsets a \code{\link{methylDiff}} object in order to get 
#' differentially methylated bases/regions
#' satisfying thresholds.
#' 
#' @param .Object  a \code{\link{methylDiff}} object
#' @param difference  cutoff for absolute value of % methylation change between test and control (default:25)
#' @param qvalue  cutoff for qvalue of differential methylation statistic (default:0.01) 
#' @param type  one of the "hyper","hypo" or "all" strings. Specifies what type of differentially menthylated bases/regions should be returned.
#'              For retrieving Hyper-methylated regions/bases type="hyper", for hypo-methylated type="hypo" (default:"all") 
#' 
#' @return a methylDiff object containing the differential methylated locations satisfying the criteria 
#' 
#' @usage get.methylDiff(.Object,difference=25,qvalue=0.01,type="all")
#' @examples
#' 
#' data(methylKit)
#' 
#' # get differentially methylated bases/regions with specific cutoffs
#' all.diff=get.methylDiff(methylDiff.obj,difference=25,qvalue=0.01,type="all")
#' 
#' # get hyper-methylated
#' hyper=get.methylDiff(methylDiff.obj,difference=25,qvalue=0.01,type="hyper")
#' 
#' # get hypo-methylated
#' hypo=get.methylDiff(methylDiff.obj,difference=25,qvalue=0.01,type="hypo")
#' 
#'
#' @export
#' @docType methods
#' @rdname get.methylDiff-methods
setGeneric(name="get.methylDiff", def=function(.Object,difference=25,qvalue=0.01,type="all") standardGeneric("get.methylDiff"))

#' @aliases get.methylDiff,methylDiff-method
#' @rdname get.methylDiff-methods
setMethod(f="get.methylDiff", signature="methylDiff", 
          definition=function(.Object,difference,qvalue,type) {
            
                    if(type=="all"){
                      new.obj=new("methylDiff",.Object[.Object$qvalue<qvalue & abs(.Object$meth.diff) > difference,],
                              sample.ids=.Object@sample.ids,assembly=.Object@assembly,context=.Object@context,
                              treatment=.Object@treatment,destranded=.Object@destranded,resolution=.Object@resolution)
                      return(new.obj)
                    }else if(type=="hyper"){
                      new.obj=new("methylDiff",.Object[.Object$qvalue<qvalue & (.Object$meth.diff) > difference,],
                              sample.ids=.Object@sample.ids,assembly=.Object@assembly,context=.Object@context,
                              treatment=.Object@treatment,destranded=.Object@destranded,resolution=.Object@resolution)
                      return(new.obj)
                    }else if(type=="hypo"){
                      new.obj=new("methylDiff",.Object[.Object$qvalue<qvalue & (.Object$meth.diff) < -1*difference,],
                              sample.ids=.Object@sample.ids,assembly=.Object@assembly,context=.Object@context,
                              treatment=.Object@treatment,destranded=.Object@destranded,resolution=.Object@resolution) 
                      return(new.obj)
                    }else{
                      stop("Wrong 'type' argument supplied for the function, it can be 'hypo', 'hyper' or 'all' ")
                    }
          }) 

##############################################################################
## PLOTTING FUNCTIONS FOR methylDiff OBJECT
##############################################################################

#' Get and plot the number of hyper/hypo methylated regions/bases per chromosome
#'
#' This function gets number of  hyper/hypo methylated regions/bases from \code{\link{methylDiff}} object. 
#' It can also plot percentages of differentially methylated bases per chromosome.
#'
#' @param x a \code{\link{methylDiff}}  object
#' @param plot TRUE|FALSE. If TRUE horizontal barplots for proportion of hypo/hyper methylated bases/regions
#' @param qvalue.cutoff  cutoff for q-value
#' @param meth.cutoff cutoff for percent methylation difference
#' @param exclude names of chromosomes to be excluded
#' @param ... extra graphical parameters to be passed to \code{\link{barplot}} function
#' 
#' @return plots a piechart or a barplot for percentage of the target features overlapping with annotation
#' 
#' @usage diffMethPerChr(x,plot=T,qvalue.cutoff=0.01, meth.cutoff=25,exclude=NULL,...)
#' @examples
#' 
#' data(methylKit)
#'  
#' # get a list of differentially methylated bases/regions per chromosome and overall
#' diffMethPerChr(methylDiff.obj, plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25,exclude=NULL)
#'
#' @export
#' @docType methods
#' @rdname diffMethPerChr-methods
setGeneric("diffMethPerChr", def=function(x,plot=T,qvalue.cutoff=0.01, meth.cutoff=25,exclude=NULL,...) standardGeneric("diffMethPerChr"))

#' @aliases diffMethPerChr,methylDiff-method
#' @rdname  diffMethPerChr-methods
setMethod("diffMethPerChr", signature(x = "methylDiff"),
                    function(x,plot,qvalue.cutoff, meth.cutoff,exclude,...){
                      x=getData(x)
                      temp.hyper=x[x$qvalue < qvalue.cutoff & x$meth.diff >= meth.cutoff,]
                      temp.hypo =x[x$qvalue < qvalue.cutoff & x$meth.diff <= -meth.cutoff,]
                      
                      dmc.hyper=100*nrow(temp.hyper)/nrow(x) # get percentages of hypo/ hyper
                      dmc.hypo =100*nrow(temp.hypo )/nrow(x)
                      
                      all.hyper.hypo=data.frame(percentage.of.hypermethylated=dmc.hyper,
                                                number.of.hypermethylated=nrow(temp.hyper),
                                                percentage.of.hypomethylated=dmc.hypo  ,
                                                number.of.hypomethylated=nrow(temp.hypo))
                      
                      # plot barplot for percentage of DMCs per chr
                      dmc.hyper.chr=merge(as.data.frame(table(temp.hyper$chr)), as.data.frame(table(x$chr)),by="Var1")
                      dmc.hyper.chr=cbind(dmc.hyper.chr,perc=100*dmc.hyper.chr[,2]/dmc.hyper.chr[,3])

                      dmc.hypo.chr=merge(as.data.frame(table(temp.hypo$chr)), as.data.frame(table(x$chr)),by="Var1")
                      dmc.hypo.chr=cbind(dmc.hypo.chr,perc=100*dmc.hypo.chr[,2]/dmc.hypo.chr[,3])

                      dmc.hypo.hyper=merge(dmc.hypo.chr[,c(1,2,4)],dmc.hyper.chr[,c(1,2,4)],by="Var1") # merge hyper hypo per chromosome
                      dmc.hypo.hyper=dmc.hypo.hyper[order(as.numeric(sub("chr","",dmc.hypo.hyper$Var1))),] # order the chromosomes
                      
                      names(dmc.hypo.hyper)=c("chr","number.of.hypomethylated","percentage.of.hypomethylated","number.of.hypermethylated","percentage.of.hypermethylated")
                      if(plot){
                        
                        if(!is.null(exclude)){dmc.hypo.hyper=dmc.hypo.hyper[! dmc.hypo.hyper$chr %in% exclude,]}
                        
                        barplot(
                          t(as.matrix(data.frame(hyper=dmc.hypo.hyper[,5],hypo=dmc.hypo.hyper[,3],row.names=dmc.hypo.hyper[,1]) ))
                          ,las=2,horiz=T,col=c("magenta","aquamarine4"),main=paste("% of hyper & hypo methylated regions per chromsome",sep=""),xlab="% (percentage)",...)
                        mtext(side=3,paste("qvalue<",qvalue.cutoff," & methylation diff. >=",meth.cutoff," %",sep="") )
                        legend("topright",legend=c("hyper","hypo"),fill=c("magenta","aquamarine4"))
                      }else{
                        
                        list(diffMeth.per.chr=dmc.hypo.hyper,diffMeth.all=all.hyper.hypo)
                        
                      }

})




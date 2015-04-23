# rewrite of diffMeth operations
require(parallel)

#' @param overdispersion MN McCullagh and Nelder (1989) recommendation for 
#'        overdispersion correction. 
#' @param effect
logRegDF<-function(df,treatment,covariates,overdispersion=c("none","MN","shrinkMN"),
                   adjust=c("SLIM","holm","hochberg","hommel","bonferroni","BH","BY","fdr","none","qvalue"),
                   effect=c("wmean","mean","predicted"),parShrinkNM=list(),
                   test=c("F","Chisq"),mc.cores=1)
{
  # get the model matrix from treatment and (optional) covariates
  vars <- as.data.frame(cbind(treatment,covariates))
  
  # get formula from model matrix
  formula <-as.formula(paste("~ ", paste(colnames(vars), collapse= "+")))
  
  # get C and T cols from methylBase data.frame-part
  Tcols=seq(7,ncol(df),by=3)
  Ccols=Tcols-1
  # get count matrix and make list
  cntlist=split(as.matrix(df[,c(Ccols,Tcols)]),1:nrow(df))
  
  # get the result of tests
  tmp=simplify2array(
    mclapply(cntlist,logReg,formula,vars,overdispersion=overdispersion,effect=effect,
             parShrinkNM=parShrinkNM,test=test,mc.cores=mc.cores))
  
  # return the data frame part of methylDiff
  tmp <- as.data.frame(t(tmp))
  res <- cbind(df[,1:4],tmp$p.value,p.adjusted(tmp$q.value,method=adjust),tmp$meth.diff.1)
  colnames(res)[5:7] <- c("pvalue","qvalue","meth.diff")
  return(res)
}

logReg<-function(counts, formula, vars, overdispersion=c("none","MN","shrinkMN"),
                 effect=c("wmean","mean","predicted"), parShrinkNM=list(), test=c("F","Chisq")){
  
  w=counts[1:(length(counts)/2)]+counts[((length(counts)/2)+1):length(counts)]
  y=counts[1:(length(counts)/2)]
  prop=y/w
  
  # if there are more variables than treatment
  if(ncol(vars)>1){
    fmlaCov<-as.formula(paste("~ ", paste(colnames(vars)[-1], collapse= "+")))
    modelCov<-model.matrix(as.formula(fmlaCov),as.data.frame(vars[,-1,drop=FALSE]) )
    # this is only with covariates
    objCov=glm.fit(modelCov,prop,weights=w,family=binomial(link=logit))
  }
  
  # full model with all variables
  modelMat<-model.matrix( formula ,as.data.frame(vars) )
  obj=glm.fit(modelMat,prop,weights=w,family=binomial(link=logit))
  
  mu=fitted(obj)
  nprm=length(obj$coef) # number of parameters fitted
  
  #get dispersion
  overdispersion <- match.arg(overdispersion)
  phi=switch(overdispersion,
             none=1,
             MN={
               mu=fitted(obj)
               uresids <- (y-w*mu)/sqrt(mu*(w-w*mu)) # pearson residuals
               phi=sum( uresids^2 )/(length(w)-nprm) # correction factor  
               ifelse(phi>1,phi,1)
             },
             shrinkMN=1)
  
  if(ncol(vars)>1){
    deviance <- objCov$deviance - obj$deviance
    ddf=objCov$df.residual-obj$df.residual
  }else{
    deviance <- obj$null.deviance - obj$deviance
    ddf=obj$df.null-obj$df.residual # difference in degrees of freedom
  }
  
  test=match.arg(test)
  
  # do F-test when overdispersion >1 given
  test=ifelse(test=="F" & phi>1,"F","Chisq")
  
  p.value=switch(test,
                 F={
                   pf(deviance/phi, ddf, (length(w)-nprm), lower.tail = FALSE)     
                 },
                 Chisq={
                   pchisq(deviance/phi, 1, lower.tail = FALSE)
                 })
  
  #calculate effect size
  effect <- match.arg(effect)
  meths=switch(effect,
               wmean={
                 ms=tapply(y, as.factor(vars$treatment), sum )
                 ws=tapply(w, as.factor(vars$treatment), sum )
                 ms/ws  
               },
               mean={
                 tapply(prop, as.factor(vars$treatment), FUN = mean)
               },
               predicted={
                 tapply(mu, as.factor(vars$treatment), FUN = mean)
               }
  )

  # calculate the difference
  # if more than two groups calculate is as abs(max difference between two groups)
  if(length(unique(treatment))>2){
    meth.diff=max(meths)-min(meths)
  }else{
    meth.diff=meths[2]-meths[1]
  }
  
  names(meths)=paste0("meth_",names(meths))
  c(meth.diff=100*meth.diff,p.value=p.value,q.value=p.value,100*meths)
}

# mutates methylDiff object
estimateShrinkageNM<-function(df,treatment,covariates,sample.size=100000){
  
}

# estimates phi and Pearson ChiSq statistic per CpG or region
# for the logistic regression model
estimateDispersionParam<-function(counts,modelMat){
  
}

fishersTestDF<-function(df,treatment,mc.cores=1){
  
  Tcols=seq(7,ncol(df),by=3) # get C and T cols
  Ccols=Tcols-1
  cntlist=split(as.matrix(df[,c(Ccols,Tcols)]),1:nrow(df)) # get count matrix and make list
  
  # get the result of tests
  res=simplify2array( mclappy(cntlst,fishersTest,treatment,mc.cores))
  
  # return the data frame part of methylDiff
  res
}

fishersTest<-function(counts,treatment){
  
  p.value=fast.fisher(matrix(counts,ncol=2) ,conf.int = F)$p.value
  m1=counts[,1]/(counts[,1]+counts[,3])
  m2=counts[,2]/(counts[,2]+counts[,4])
  meths=c(m1,m2)
  names(meths)=paste0("meth_",treatment)
  c(meth.diff=m2-m1,p.value=p.value,q.value=q.value,meths)
  
}


# A Faster version of fishers exact test
# it basically chops of some parts of the default base function, invoked same
# way
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

#p.adjust2<-function(pvalues,method=c("SLIM","BH","BY","qvalue")){}
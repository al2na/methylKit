## PART THAT DEALS WITH DIFFERENTIAL METHYLATION CALCULATIONS

##############################################################################
## S3 functions to be used in S4 stuff
##############################################################################

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


# New calulateDiffMeth-part

logReg<-function(counts, formula, vars, treatment, overdispersion=c("none","MN","shrinkMN"),
                 effect=c("wmean","mean","predicted"), parShrinkNM=list(), test=c("F","Chisq")){
  
  # correct counts and treatment factor for NAs in counts
  treatment<-ifelse(is.na(counts),NA,treatment)[1:length(treatment)]
  treatment<-treatment[!is.na(treatment)]
  counts<-counts[!is.na(counts)]
  
  w=counts[1:(length(counts)/2)]+counts[((length(counts)/2)+1):length(counts)]
  y=counts[1:(length(counts)/2)]
  prop=y/w
  
  # get the model matrix from treatment and (optional) covariates
  vars <- as.data.frame(cbind(treatment,vars))
  
  # get formula from model matrix
  formula <-as.formula(paste("~ ", paste(colnames(vars), collapse= "+")))
  
  # if there are covariates other than treatment
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
             shrinkMN={
               mu=fitted(obj)
               uresids <- (y-w*mu)/sqrt(mu*(w-w*mu)) # pearson residuals
               phi=sum( uresids^2 )/(length(w)-nprm) # correction factor
               
               # change phi with parameters from parShrinkMN
               df.prior=parShrinkMN$df.prior
               var.prior=parShrinkMN$var.prior
               df.total=(length(w)-nprm)+df.prior
               phi=((length(w)-nprm)*phi + df.prior*var.prior)/df.total
               ifelse(phi>1,phi,1)
             })
  
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
                 #ms=tapply(y, as.factor(vars$treatment), sum )
                 ms=tapply(y,treatment,sum)
                 #ws=tapply(w, as.factor(vars$treatment), sum )
                 ws=tapply(w,treatment,sum)
                 ms/ws  
               },
               mean={
                 tapply(prop, treatment, FUN = mean)
               },
               predicted={
                 tapply(mu, treatment, FUN = mean)
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

# estimation of shrinkage parameters
estimateShrinkageMN<-function(cntlist,treatment,covariates,
                              sample.size=100000,mc.cores=1){
  
  message("Estimating shrinkage for scaling factor phi...")
  
  # get formula and construct model matrix
  vars <- as.data.frame(cbind(treatment,covariates))
  formula <-as.formula(paste("~ ", paste(colnames(vars), collapse= "+")))
  modelMat<-model.matrix( formula ,as.data.frame(vars) )
  
  # check number of samples
  sample.size <- ifelse(length(cntlist)<=100000,length(cntlist),sample.size)
  
  # calculate phis up to the first 100000 sites (stabilizing point)
  estimation=simplify2array(
    mclapply(cntlist[1:sample.size],estimatePhi,modelMat=modelMat,treatment=treatment,
             mc.cores=mc.cores))
  
  # for each phi, take the correct df (depending on number of model parameters)
  phis<-estimation[1,]
  df<-estimation[2,]
  
  # squeeze sample variances 
  shr=squeezeVar(phis,df)
  
  # output prior df and variances (to be used later as input for parShrinkMN)
  list(df.prior=shr$df.prior,var.prior=shr$var.prior)
}

estimatePhi<-function(counts,modelMat,treatment){
  
  # correct counts and treatment factor for NAs in counts
  treatment<-ifelse(is.na(counts),NA,treatment)[1:length(treatment)]
  treatment<-treatment[!is.na(treatment)]
  counts<-counts[!is.na(counts)]
  
  n=counts[1:(length(counts)/2)]+counts[((length(counts)/2)+1):length(counts)]
  y=counts[1:(length(counts)/2)]
  prop=y/n
  
  glmfit=glm.fit(modelMat,prop,weights=n,family=binomial(link=logit))
  
  # fit glm
  mu <- fitted(glmfit)
  
  # calculate and record results
  resids <- (y-n*mu)/sqrt(mu*(n-n*mu))

  # get phi correction coefficients 
  phi <- sum( resids^2 )/(length(n)-2)
  c(phi,(length(n)-2))
}

# end of S3 functions

##############################################################################
## S4 OBJECTS
##############################################################################


#' An S4 class that holds differential methylation information
#'
#' This class is designed to hold statistics and locations for differentially 
#' methylated regions/bases. It extends \code{\link{data.frame}} class.
#' \code{\link[methylKit]{calculateDiffMeth}} function returns an object 
#' with \code{methylDiff} class.
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
#'  Subsetting by \code{x[i,]} will produce a new object if subsetting is done 
#'  on rows. Column subsetting is not directly allowed to prevent errors in the 
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
#' @param .Object a methylBase or methylBaseDB object to calculate differential methylation                    
#' @param covariates a data.frame containing covariates, which should be included in the test.                   
#' @param overdispersion If set to "none"(default), no overdispersion correction will be attempted.
#'              If set to "MN", basic overdispersion correction will be applied.
#'              (EXPERIMENTAL: If set to "shrinkMN", overdisperison correction with squeezeVar() 
#'              from the limma-package will be applied (not thoroughly tested as of yet).
#' @param adjust different methods to correct the p-values for multiple testing. 
#'              Default is "SLIM" from methylKit. For "qvalue" please see \code{\link[qvalue]{qvalue}} 
#'              and for all other methods see \code{\link[stats]{p.adjust}}.
#' @param effect method to calculate the mean methylation different between groups 
#'              using read coverage as weights (default). When set to "mean", the generic mean is applied
#'              and when set to "predicted", a logistic model is used instead.
#' @param parShrinkMN a list for squeezeVar(). (NOT IMPLEMENTED)
#' @param test the statistical test used to determine the methylation differences. 
#'              The Chisq-test is used by default, while the F-test can be chosen 
#'              if overdispersion control ist applied.
#' @param mc.cores integer denoting how many cores should be used for parallel
#'              differential methylation calculations (can only be used in
#'              machines with multiple cores).
#' @param slim If set to FALSE, \code{adjust} will be set to "BH" (default behaviour of earlier versions)
#' @param weighted.mean If set to FALSE, \code{effect} will be set to "mean" (default behaviour of earlier versions)  
#' @param chunk.size Number of rows to be taken as a chunk for processing the \code{methylBaseDB} objects (default: 1e6)
#' @param save.db A Logical to decide whether the resulting object should be saved as flat file database or not, default: explained in Details sections  
#' @param ... optional Arguments used when save.db is TRUE
#'            
#'            \code{suffix}
#'                  A character string to append to the name of the output flat file database, 
#'                  only used if save.db is true, default actions: append \dQuote{_filtered} to current filename 
#'                  if database already exists or generate new file with filename \dQuote{sampleID_filtered}
#'                  
#'            \code{dbdir} 
#'                  The directory where flat file database(s) should be stored, defaults
#'                  to getwd(), working directory for newly stored databases
#'                  and to same directory for already existing database
#'                  
#            \code{dbtype}
#                  The type of the flat file database, currently only option is "tabix"
#                  (only used for newly stored databases)             
#'                    
#' @usage calculateDiffMeth(.Object,covariates,overdispersion=c("none","MN","shrinkMN"),
#'                          adjust=c("SLIM","holm","hochberg","hommel","bonferroni","BH",
#'                          "BY","fdr","none","qvalue"),effect=c("wmean","mean","predicted"), 
#'                          parShrinkNM=list(),test=c("F","Chisq"),mc.cores=1,slim=TRUE, 
#'                          weighted.mean=TRUE, chunk.size, save.db,...)
#' 
#' @examples
#' 
#' data(methylKit)
#' 
#' # The Chisq-test will be applied when no overdispersion control is chosen.
#' my.diffMeth=calculateDiffMeth(methylBase.obj,covariates=NULL,overdispersion=c("none"),
#'                               adjust=c("SLIM"),effect=c("wmean"),parShrinkMN=list(),
#'                               test=c("Chisq"),mc.cores=1)
#' 
#' # pool samples in each group
#' pooled.methylBase=pool(methylBase.obj,sample.ids=c("test","control"))
#'  
#' # After applying the pool() function, there is one sample in each group.
#' # The Fisher's exact test will be applied for differential methylation.
#' my.diffMeth2=calculateDiffMeth(pooled.methylBase,covariates=NULL,overdispersion=c("none"),
#'                                adjust=c("SLIM"),effect=c("wmean"),test=c("F"))
#'                                
#' # Covariates and overdispersion control:
#' # generate a methylBase object with age as a covariate
#' covariates=data.frame(age=c(30,80,30,80))
#' sim.methylBase<-dataSim(replicates=4,sites=1000,treatment=c(1,1,0,0),
#'                         covariates=covariates,
#'                         sample.ids=c("test1","test2","ctrl1","ctrl2"))
#' 
#' # Apply overdispersion correction and include covariates 
#' # in differential methylation calculations.
#' my.diffMeth3<-calculateDiffMeth(sim.methylBase,
#'                                 covariates=covariates,
#'                                 overdispersion="MN",test="Chisq",mc.cores=1)
#'                                
#' @return a methylDiff object containing the differential methylation 
#'                      statistics and locations
#' @section Details:
#' Covariates can be included in the analysis. The function will then try to separate the 
#' influence of the covariates from the treatment effect via a linear model.\cr
#' The Chisq-test is used per default only when no overdispersion correction is applied.
#' If overdispersion correction is applied, the function automatically switches to the 
#' F-test. The Chisq-test can be manually chosen in this case as well, but the F-test only 
#' works with overdispersion correction switched on.
#' 
#' The parameter \code{chunk.size} is only used when working with \code{methylBaseDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects as \code{methylBaseDB}, 
#' while being per default FALSE for \code{methylBase}. If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the database allows the calculation in-memory, 
#' then you might want to change the value of this parameter.
#' 
#' @references Altuna Akalin, Matthias Kormaksson, Sheng Li,
#'             Francine E. Garrett-Bakelman, Maria E. Figueroa, Ari Melnick, 
#'             Christopher E. Mason. (2012). 
#'             "methylKit: A comprehensive R package for the analysis 
#'             of genome-wide DNA methylation profiles." Genome Biology. 
#' @seealso \code{\link[methylKit]{pool}}, \code{\link[methylKit]{reorganize}}
#'          \code{\link[methylKit]{dataSim}}
#' 
#' @export
#' @docType methods
#' @aliases calculateDiffMeth,methylBase-method
#' @rdname calculateDiffMeth-methods

setGeneric("calculateDiffMeth", function(.Object,covariates=NULL,
                                         overdispersion=c("none","MN","shrinkMN"),
                                         adjust=c("SLIM","holm","hochberg","hommel",
                                                  "bonferroni","BH","BY","fdr",
                                                  "none","qvalue"),
                                         effect=c("wmean","mean","predicted"),
                                         parShrinkNM=list(),
                                         test=c("F","Chisq"),mc.cores=1,slim=TRUE,
                                         weighted.mean=TRUE,
                                         chunk.size=1e6,save.db=FALSE,...) 
  standardGeneric("calculateDiffMeth"))

setMethod("calculateDiffMeth", "methylBase",
          function(.Object,covariates,overdispersion,
                   adjust,effect,parShrinkNM,
                   test,mc.cores,slim,weighted.mean,save.db=FALSE,...){
            
            # extract data.frame from methylBase
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
            
            # add backwards compatibility with old parameters
            if(slim==FALSE) adjust="BH" else adjust=adjust
            if(weighted.mean==FALSE) effect="mean" else effect=effect
            
            vars <- covariates
            
            # get C and T cols from methylBase data.frame-part
            Tcols=seq(7,ncol(subst),by=3)
            Ccols=Tcols-1
            
            #### check if covariates+intercept+treatment more than replicates ####
            if(!is.null(covariates)){if(ncol(covariates)+2 >= length(Tcols)){stop("Too many covariates/too few replicates.")}}
            
            # get count matrix and make list
            cntlist=split(as.matrix(subst[,c(Ccols,Tcols)]),1:nrow(subst))
            
            # call estimate shrinkage before logReg
            if(overdispersion=="shrinkMN"){
              parShrinkMN<-estimateShrinkageMN(cntlist,treatment=.Object@treatment,covariates=vars,
                                               sample.size=100000,mc.cores=mc.cores)
            }
            
            # get the result of tests
            tmp=simplify2array(
              mclapply(cntlist,logReg,formula,vars,treatment=.Object@treatment,overdispersion=overdispersion,effect=effect,
                       parShrinkMN=parShrinkMN,test=test,mc.cores=mc.cores))
            
            # return the data frame part of methylDiff
            tmp <- as.data.frame(t(tmp))
            x=data.frame(subst[,1:4],tmp$p.value,p.adjusted(tmp$q.value,method=adjust),meth.diff=tmp$meth.diff.1,stringsAsFactors=FALSE)
            colnames(x)[5:7] <- c("pvalue","qvalue","meth.diff")
            
            if(!save.db) {
              obj=new("methylDiff",x,sample.ids=.Object@sample.ids,assembly=.Object@assembly,context=.Object@context,
                      destranded=.Object@destranded,treatment=.Object@treatment,resolution=.Object@resolution)
              obj
            } else {
              
              # catch additional args 
              args <- list(...)
              
              if( !( "dbdir" %in% names(args)) ){
                dbdir <- .check.dbdir(getwd())
              } else { dbdir <- .check.dbdir(args$dbdir) }
              #                         if(!( "dbtype" %in% names(args) ) ){
              #                           dbtype <- "tabix"
              #                         } else { dbtype <- args$dbtype }
              if(!( "suffix" %in% names(args) ) ){
                suffix <- "_diffMeth"
              } else { 
                suffix <- paste0("_",args$suffix)
              }
              
              # create methylDiffDB
              makeMethylDiffDB(df=x,dbpath=dbdir,dbtype="tabix",sample.ids=.Object@sample.ids,
                               assembly=.Object@assembly,context=.Object@context,
                               destranded=.Object@destranded,treatment=.Object@treatment,
                               resolution=.Object@resolution,suffix=suffix )
            }

        }
)



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


#' @rdname getTreatment-methods
#' @aliases getTreatment,methylDiff-method
setMethod("getTreatment", signature = "methylDiff", function(x) {
  return(x@treatment)
})

#' @rdname getTreatment-methods
#' @aliases getTreatment,methylDiff-method
setReplaceMethod("getTreatment", signature = "methylDiff", function(x, value) {
  
  if(! ( length(x@treatment) == length(value) ) ){
    stop("The new treatment vector is not valid, check the length of input")
  } else {
    x@treatment <- value
    return(x)
  }
  
})


#' @rdname getSampleID-methods
#' @aliases getSampleID,methylDiff-method
setMethod("getSampleID", signature = "methylDiff", function(x) {
  return(x@sample.id)
})

#' @rdname getSampleID-methods
#' @aliases getSampleID,methylDiff-method
setReplaceMethod("getSampleID", signature = "methylDiff", function(x, value) {
  
  if(! ( length(value) == length(x@sample.ids) ) ){
    stop("The vector of new sample ids is not valid, check the length of input")
  } else {
    x@sample.ids <- value
    return(x)
  }
  
})


##############################################################################
## CONVERTOR FUNCTIONS FOR methylDiff OBJECT
##############################################################################

setAs("methylDiff", "GRanges", function(from)
{
  
  GRanges(seqnames=as.character(from$chr),ranges=IRanges(start=from$start, end=from$end),
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
# @aliases [,methylDiff-method
#' @aliases extract,methylDiff-method
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

#' @aliases selectByOverlap,methylDiff-method
#' @rdname selectByOverlap-methods
setMethod("selectByOverlap", "methylDiff",
          function(object, ranges){
            
            if(missing(ranges) | class(ranges)!="GRanges") {
              stop("No ranges specified or given ranges object not of class GRanges, please check your input!")
            }
            hits <- findOverlaps(ranges,as(object,"GRanges"))@subjectHits
            
            return(object[hits])
          }
)


#' get differentially methylated regions/bases based on cutoffs 
#' 
#' The function subsets a \code{\link{methylDiff}} or \code{\link{methylDiffDB}} object in order to get 
#' differentially methylated bases/regions
#' satisfying thresholds.
#' 
#' @param .Object  a \code{\link{methylDiff}} or \code{\link{methylDiffDB}} object
#' @param difference  cutoff for absolute value of methylation percentage change between test and control (default:25)
#' @param qvalue  cutoff for qvalue of differential methylation statistic (default:0.01) 
#' @param type  one of the "hyper","hypo" or "all" strings. Specifies what type of differentially menthylated bases/regions should be returned.
#'              For retrieving Hyper-methylated regions/bases type="hyper", for hypo-methylated type="hypo" (default:"all") 
#' @param chunk.size Number of rows to be taken as a chunk for processing the \code{methylDiffDB} objects (default: 1e6)
#' @param save.db A Logical to decide whether the resulting object should be saved as flat file database or not, default: see Details  
#' @param ... optional Arguments used when save.db is TRUE
#'            
#'            \code{suffix}
#'                  A character string to append to the name of the output flat file database, 
#'                  only used if save.db is true, default actions: append \dQuote{_filtered} to current filename 
#'                  if database already exists or generate new file with filename \dQuote{sampleID_filtered}
#'                  
#'            \code{dbdir} 
#'                  The directory where flat file database(s) should be stored, defaults
#'                  to getwd(), working directory for newly stored databases
#'                  and to same directory for already existing database
#'                  
#            \code{dbtype}
#                  The type of the flat file database, currently only option is "tabix"
#                  (only used for newly stored databases)
#' 
#' @return a methylDiff or methylDiffDB object containing the differential methylated locations satisfying the criteria 
#' 
#' @usage get.methylDiff(.Object,difference=25,qvalue=0.01,type="all",chunk.size,save.db,...)
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
#' @section Details:
#' The parameter \code{chunk.size} is only used when working with \code{methylDiffDB} objects, 
#' as they are read in chunk by chunk to enable processing large-sized objects which are stored as flat file database.
#' Per default the chunk.size is set to 1M rows, which should work for most systems. If you encounter memory problems or 
#' have a high amount of memory available feel free to adjust the \code{chunk.size}.
#' 
#' The parameter \code{save.db} is per default TRUE for methylDB objects as \code{methylDiffDB}, 
#' while being per default FALSE for \code{methylDiff}. If you wish to save the result of an 
#' in-memory-calculation as flat file database or if the size of the database allows the calculation in-memory, 
#' then you might want to change the value of this parameter.
#'
#' @export
#' @docType methods
#' @rdname get.methylDiff-methods
setGeneric(name="get.methylDiff", def=function(.Object,difference=25,qvalue=0.01,type="all",chunk.size=1e6,save.db=FALSE,...) standardGeneric("get.methylDiff"))

#' @aliases get.methylDiff,methylDiff-method
#' @rdname get.methylDiff-methods
setMethod(f="get.methylDiff", signature="methylDiff", 
          definition=function(.Object,difference,qvalue,type,save.db=FALSE,...) {
            
            if(!save.db) {
            
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
            
            } else {
              
              # catch additional args 
              args <- list(...)
              
              if( !( "dbdir" %in% names(args)) ){
                dbdir <- .check.dbdir(getwd())
              } else { dbdir <- .check.dbdir(args$dbdir) }
              #                         if(!( "dbtype" %in% names(args) ) ){
              #                           dbtype <- "tabix"
              #                         } else { dbtype <- args$dbtype }
              if(!( "suffix" %in% names(args) ) ){
                suffix <- paste0("_",type)
              } else { 
                suffix <- paste0("_",args$suffix)
              }
              
              # create methylBaseDB
              if(type=="all"){
                new.obj= makeMethylDiffDB(df=.Object[.Object$qvalue<qvalue & abs(.Object$meth.diff) > difference,],
                                   dbpath=dbdir,dbtype="tabix",sample.ids=.Object@sample.ids,
                                   assembly=.Object@assembly,context=.Object@context,
                                   destranded=.Object@destranded,treatment=.Object@treatment,
                                   resolution=.Object@resolution,suffix=suffix )
              }else if(type=="hyper"){
              new.obj= makeMethylDiffDB(df=.Object[.Object$qvalue<qvalue & (.Object$meth.diff) > difference,],
                                   dbpath=dbdir,dbtype="tabix",sample.ids=.Object@sample.ids,
                                   assembly=.Object@assembly,context=.Object@context,
                                   destranded=.Object@destranded,treatment=.Object@treatment,
                                   resolution=.Object@resolution,suffix=suffix )
              }else if(type=="hypo"){
                new.obj= makeMethylDiffDB(df=.Object[.Object$qvalue<qvalue & (.Object$meth.diff) < -1*difference,],
                                   dbpath=dbdir,dbtype="tabix",sample.ids=.Object@sample.ids,
                                   assembly=.Object@assembly,context=.Object@context,
                                   destranded=.Object@destranded,treatment=.Object@treatment,
                                   resolution=.Object@resolution,suffix=suffix )
              }else{
                stop("Wrong 'type' argument supplied for the function, it can be 'hypo', 'hyper' or 'all' ")
              }
              return(new.obj) 
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
            
            all.hyper.hypo=data.frame(number.of.hypermethylated=nrow(temp.hyper),
                                      percentage.of.hypermethylated=dmc.hyper,
                                      number.of.hypomethylated=nrow(temp.hypo),
                                      percentage.of.hypomethylated=dmc.hypo)
            
            # plot barplot for percentage of DMCs per chr
            dmc.hyper.chr=merge(as.data.frame(table(temp.hyper$chr)), as.data.frame(table(x$chr)),by="Var1")
            dmc.hyper.chr=cbind(dmc.hyper.chr,perc=100*dmc.hyper.chr[,2]/dmc.hyper.chr[,3])
            
            dmc.hypo.chr=merge(as.data.frame(table(temp.hypo$chr)), as.data.frame(table(x$chr)),by="Var1")
            dmc.hypo.chr=cbind(dmc.hypo.chr,perc=100*dmc.hypo.chr[,2]/dmc.hypo.chr[,3])
            
            dmc.hyper.hypo=merge(dmc.hyper.chr[,c(1,2,4)],dmc.hypo.chr[,c(1,2,4)],by="Var1") # merge hyper hypo per chromosome
            dmc.hyper.hypo=dmc.hyper.hypo[order(as.numeric(sub("chr","",dmc.hyper.hypo$Var1))),] # order the chromosomes
            
            names(dmc.hyper.hypo)=c("chr","number.of.hypermethylated","percentage.of.hypermethylated","number.of.hypomethylated","percentage.of.hypomethylated")
            if(plot){
              
              if(!is.null(exclude)){dmc.hyper.hypo=dmc.hyper.hypo[! dmc.hyper.hypo$chr %in% exclude,]}
              
              barplot(
                t(as.matrix(data.frame(hyper=dmc.hyper.hypo[,3],hypo=dmc.hyper.hypo[,5],row.names=dmc.hyper.hypo[,1]) ))
                ,las=2,horiz=T,col=c("magenta","aquamarine4"),main=paste("% of hyper & hypo methylated regions per chromsome",sep=""),xlab="% (percentage)",...)
              mtext(side=3,paste("qvalue<",qvalue.cutoff," & methylation diff. >=",meth.cutoff," %",sep="") )
              legend("topright",legend=c("hyper","hypo"),fill=c("magenta","aquamarine4"))
            }else{
              
              list(diffMeth.per.chr=dmc.hyper.hypo,diffMeth.all=all.hyper.hypo)
              
            }
            
          })


# functions used to simulate methylation data

#' Simulate DNA methylation data
#' 
#' The function simulates DNA methylation data from multiple samples.
#' See references for detailed explanation on statistics.
#'
#' @param replicates  the number of samples that should be simulated.
#' @param sites       the number of CpG sites per sample.
#' @param treatment   a vector containing treatment information.
#' @param percentage  the proportion of sites which should be affected by the 
#' treatment.
#' @param effect      a number between 0 and 100 specifying the effect size of 
#'                    the treatment. This is essentially describing the average
#'                    percent methylation difference between differentially
#'                    methylated bases.See 'Examples' and 'Details'.
#' @param alpha       shape1 parameter for beta distribution (used for 
#'                    initial sampling of methylation  proportions)
#' @param beta        shape2 parameter for beta distribution (used for 
#'                    initial sampling of methylation  proportions) 
#' @param theta       dispersion parameter for beta distribution (initial 
#'                    sampling of methylation  proportions)
#' @param covariates  a data.frame containing covariates (optional)
#' @param sample.ids  will be generated automatically from \code{treatment}, 
#'                    but can be 
#'                    overwritten by a character vector containing sample names.
#' @param assembly    the assembly description (e.g. "hg18").Only
#'                    needed for book keeping. 
#' @param context     the experimanteal context of the data (e.g. "CpG"). Only
#'                    needed for book keeping. 
#' @param add.info if set to TRUE, the output will be a list with the first 
#'                    element being 
#'                    the methylbase object and a vector of indices that 
#'                    indicate which CpGs should be differentially
#'                    methylated. This vector can be used to subset
#'                    simulated methylBase or methylDiff object with
#'                    differentially methylated bases.
#'
#' @examples
#' 
#' data(methylKit)
#' 
#' # Simulate data for 4 samples with 20000 sites each.
#' # The methylation in 10% of the sites are elevated by 25%.
#' my.methylBase=dataSim(replicates=4,sites=2000,treatment=c(1,1,0,0),
#' percentage=10,effect=25)
#' 
#'
#'
#' @return a methylBase object containing simulated methylation data, 
#' or if add.info=TRUE a list containing the methylbase object and 
#' the indices of all treated sites (differentially methylated bases or regions) 
#' as the second element.
#' 
#' @section Details:
#' The function uses 
#' a Beta distribution to simulate the methylation proportion 
#' background across all samples. 
#' The parameters \code{alpha}, \code{beta} used in a beta distribution to draw
#' methylation proportions,\eqn{\mu}, from a typical bimodal distribution. 
#' For each initial methylation proportion drawn using the parameters above, 
#' a range of 
#' methylation proportions is distributed around the original \eqn{\mu} with
#' overdispersion parameter \eqn{\theta}, this is using an alternative 
#' parameterization of Beta distribution:   \eqn{Beta(\mu,\theta)}. 
#' The parameters \code{percentage} and \code{effect} determine the proportion 
#' of sites that are 
#' affected by the treatment (meaning differential sites) and the strength of
#' this influence, respectively. \code{effect} is added on top of \eqn{\mu} for
#' the CpGs that are affected by the treament. The affected group of samples 
#' for that
#' particular CpG will now be distributed by \eqn{Beta(\mu+effect,\theta)}.
#' The coverage is modeled with a negative binomial distribution, using
#' \code{rnbinom} function with \code{size=1} and \code{prob=0.01}. 
#' The additional information needed for a valid methylBase object, such as 
#' CpG start, end and strand, is generated as "dummy values", 
#' but can be overwritten as needed.
#' 
#' @export
#' @docType methods
#' @aliases dataSim
#' @rdname dataSim-methods

dataSim <- function(replicates,sites,treatment,percentage=10,effect=25,
                     alpha=0.4,beta=0.5,theta=10,
                     covariates=NULL,sample.ids=NULL,assembly="hg18",
                     context="CpG",add.info=FALSE){
  
  # check if length(treatment) == # replicates
  if(length(treatment) != replicates){
    stop("treatment and replicates must be of same length")} 
  # check if # covariates == # replicates
  if(!is.null(covariates)){
    if(nrow(covariates)!=replicates){
      stop("nrow(covariates) must be equal to replicates")}
  }
  
  # create sample.ids (if not given by user)
  if(is.null(sample.ids)){
    sample.ids<-ifelse(treatment==1,paste0("test",cumsum(treatment)),
                       paste0("ctrl",cumsum(!treatment)))
  }
  
  # create data.frame
  raw <- matrix(ncol=replicates*3,nrow=sites)
  index<-seq(1,replicates*3,by=3) # for easier access of TCols, CCols, coverage
  
  # draw substitution probabilities from beta distribution 
  # (same for all samples)getData(a[[1]])[,6]
  x <- rbeta(sites,alpha,beta)
  #x<- rbeta(sites,alpha/4,beta/4)
  
  
  # get treatment and covariate indices for all samples
  treatment_indices<-sample((1:sites)[x < (1-(effect/100))],
                            size=sites*(percentage/100))
  covariate_indices<-sample(1:sites,size=sites*0.05)
  
  # get treatment effect indices
  addinfo<-rep(0,sites)
  
  # if more than one effect size is supplied, randomize effect sizes
  effects <- if(length(effect)==1) rep(effect,length(x)) else sample(effect,
                                                                     length(x),replace=TRUE)
  
  # fill data.frame with raw counts for each sample
  for(i in 1:replicates){
    
    # draw coverage from neg. binomial distribution
    coverage <- rnbinom(sites,1,0.01)
    # fix sites w/o reads 
    coverage <- ifelse(coverage<10,10,coverage)
    
    # fill in coverage:
    raw[,index[i]] <- coverage
    
    # fill in CCols: coverage * percentage of Cs
    raw[,index[i]+1] <- ceiling(coverage * rbetabinom(n=sites,prob=x,
                                                      size=50,theta=theta)/50)
    
    # add treatment information
    if(treatment[i]==1){
      # increase base probabilities
      y<-x+(effects/100)
      
      ifelse(y>1,1-x,effect)
      
      # update treatment indices (if effect result is > 1, take 1-x as real effect size)
      addinfo<-ifelse(y>1,1-x,effect)
      # effect size for all non-treated indices are set to 0
      addinfo[-treatment_indices]<-0
      
      # correct base probabilities > 1
      y<-ifelse(y>1,1,y)
      
      # TCols: coverage * percentage of Cs with modified probability y
      raw[treatment_indices,index[i]+1] <- 
        ceiling(coverage[treatment_indices] * rbetabinom(
          n=length(treatment_indices),prob=y[treatment_indices],
          size=50,theta=theta)/50)
    }
    
    # add covariate information
    if(!is.null(covariates[i,,drop=FALSE])){
      # CCols: coverage * percentage of Cs with modified probability
      raw[covariate_indices,index[i]+1] <- 
        ceiling(coverage[covariate_indices] * rbetabinom(
          n=length(covariate_indices),
          prob=influence(p=x[covariate_indices],
                         x=rep(covariates[i,],
                               times=length(covariate_indices))),
          size=50,theta=theta)/50)
    }
    
    # fill in CCols: coverage - Cs
    raw[,index[i]+2] <- coverage - raw[,index[i]+1]
  }
  
  # name data.frame and return it
  df<-raw
  df=as.data.frame(df) 
  
  #add dummy chromosome, start, end and strand info 
  info<-data.frame(chr=rep("chr1",times=sites),start=1:sites,
                   end=2:(sites+1),strand=rep("+",times=sites))
  df<-cbind(info,df)
  
  # get indices of coverage,numCs and numTs in the data frame 
  coverage.ind=seq(5,by=3,length.out=replicates)
  numCs.ind   =coverage.ind+1
  numTs.ind   =coverage.ind+2
  
  # change column names
  names(df)[coverage.ind]=paste(c("coverage"),1:replicates,sep="" )
  names(df)[numCs.ind]   =paste(c("numCs"),1:replicates,sep="" )
  names(df)[numTs.ind]   =paste(c("numTs"),1:replicates,sep="" )
  
  #make methylbase object and return the object
  obj=new("methylBase",(df),sample.ids=sample.ids,
          assembly=assembly,context=context,treatment=treatment,
          coverage.index=coverage.ind,numCs.index=numCs.ind,numTs.index=numTs.ind,
          destranded=FALSE,resolution="base")
  if(add.info==FALSE){
    obj
  }
  else{
    addinfo=treatment_indices
    list(obj,addinfo)
  }
}



# hardcoded function to transform prob p via logistic function for covariate
influence <- function(p,x=NULL,b0=0.1,b1=0.1,k=100){
  
  if(is.null(x)){
    p<-p
  }
  else{
    # logistc function for covariate x 
    y <- exp(b0 + b1*x) / (k + exp(b0 + b1*x))
    # take mean of both probabilies
    p <- (p+y)/2 
  }
}





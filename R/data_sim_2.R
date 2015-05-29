require(emdbook)

#' @param replicates  number of replicates
#' @param sites       number of CpG sites per replicate
#' @param percentage  percentage of sites which are affected by treatment
#' @param effect      effect size of treatment
#' @param alpha       shape1 parameter for beta distribution (used for substitution probabilites)
#' @param beta        shape2 parameter for beta distribution (used for substitution probabilites)
#' @param theta       dispersion parameter for beta distribution (used for substitution probabilites)
#' @param treatment   vector containing treatment information
#' @param covariates  data.frame containing covariates
#' @param experiment  ... not implemented
#' 
#' @return returns a data.frame containing simulated methylation data
data_sim <- function(replicates,sites=20000,percentage=20,effect=20,alpha=0.4,beta=0.5,theta=5,
                     treatment,covariates=NULL,experiment=NULL){
  
  # check if length(treatment) == # replicates
  if(length(treatment) != replicates){stop("treatment and replicates must be of same length")} 
  # check if # covariates == # replicates
  if(!is.null(covariates)){
    if(nrow(covariates)!=replicates){stop("nrow(covariates) must be equal to replicates")}
  }
  
  # create data.frame
  raw <- matrix(ncol=replicates*3,nrow=sites)
  index<-seq(1,replicates*3,by=3) # for easier access of TCols, CCols, coverage
  
  # draw substitution probabilities from beta distribution (same for all samples)
  x <- rbeta(sites,alpha,beta)
  
  # fill data.frame with raw counts for each sample
  for(i in 1:replicates){
    
    # draw coverage from neg. binomial distribution
    coverage <- rnbinom(sites,1,0.01)
    # fix sites w/o reads 
    coverage <- ifelse(coverage<10,10,coverage)
    
    # fill in coverage:
    raw[,index[i]] <- coverage
    
    # fill in TCols: coverage * percentage of Ts
    raw[,index[i]+2] <- ceiling(coverage * rbetabinom(n=sites,prob=x,size=50,theta=theta)/50)
    # add treatment information
    if(treatment[i]==1){
      # randomly choose sites
      treatment_indices<-sample(sites,size=sites*(percentage/100))
      # increase base probabilities
      y<-x+(effect/100);y<-ifelse(y>1,1,y)
      # TCols: coverage * percentage of Ts with modified probability y
      raw[treatment_indices,index[i]+2] <- 
        ceiling(coverage[treatment_indices] * rbetabinom(n=1,prob=y[treatment_indices],size=50,theta=theta)/50)
    } 
    # add covariate information
    if(!is.null(covariates[i,,drop=FALSE])){
      # randomly choose 5% of all sites
      covariate_indices<-sample(sites,size=sites*0.05)
      # TCols: coverage * percentage of Ts with modified probability
      raw[covariate_indices,index[i]+2] <- 
        ceiling(coverage[covariate_indices] * rbetabinom(n=1,prob=influence(p=x[covariate_indices],x=rep(covariates[i,],times=length(covariate_indices))),size=50,theta=theta)/50)
    }
    
    # fill in CCols: coverage - Cs
    raw[,index[i]+1] <- coverage - raw[,index[i]+2]
  }
  
  # name data.frame and return it
  res<-raw
  colnames(res) <- paste0(c("coverage","numCs","numTs"),rep(1:replicates,each=3))
  as.data.frame(res)
}

### hardcoded function to transform prob p via logistic function for covariate
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

#### Testing ####
set.seed(100)
replicates=4
sites=20000
percentage=0.2
effect=0.2
alpha=0.4
beta=0.5
theta=5
treatment=c(0,0,1,1)
covariates=data.frame(age=c(20,60,30,50))
experiment=NULL

# test-run
result <- data_sim(replicates=4,sites=20000,percentage=20,effect=20,alpha=0.4,beta=0.5,theta=5,
                   treatment=c(0,0,1,1),covariates=covariates,experiment=NULL)

# checks validity of the data.frame supposedly returned by segMeth or segDiffMeth
check.seg.df<-function(seg.df,diff=FALSE,meth=FALSE)
{

  if( diff==FALSE & meth==FALSE){

    if( ! all(names(seg.df) %in% c("chr","start","end","num.cov.base","avg.meth.diff","med.meth.diff","hmm.state","avg.meth","med.meth" ))  ){
      stop("Provided data frame doesn't look like the output of segMeth() or segDiffMeth() functions", .call=FALSE)}
    if( ncol(seg.df) != 7 ){
      stop("Provided data frame doesn't look like the output of segMeth() or segDiffMeth() functions", .call=FALSE)}
  }else if(diff){
    if( ! all(names(seg.df) %in% c("chr","start","end","num.cov.base","avg.meth.diff","med.meth.diff","hmm.state" ))  ){
      stop("Provided data frame doesn't look like the output of segDiffMeth() function", .call=FALSE)}
    if( ncol(seg.df) != 7 ){
      stop("Provided data frame doesn't look like the output of segDiffMeth() function", .call=FALSE)}
  }else if(meth){
    if( ! all(names(seg.df) %in% c("chr","start","end","num.cov.base","hmm.state","avg.meth","med.meth" ))  ){
      stop("Provided data frame doesn't look like the output of segMeth() functions", .call=FALSE)}
    if( ncol(seg.df) != 7 ){
      stop("Provided data frame doesn't look like the output of segMeth() functions", .call=FALSE)}

  }
  
}

# function initializes parameters for semi-supervised hmm
initializeParam<-function(obs,hi.th=5,lo.th=-5){

  my.states=rep(1,length(obs))
  my.states[obs>hi.th]=2
  my.states[obs<lo.th]=3

  t.df=cbind(my.states[-length(my.states)],my.states[-1])

  # make transition matrix
  #mat=matrix(ncol=3,nrow=3)
  #for(k in 1:3){
  #  for(l in 1:3){
  #    mat[k,l]=sum( t.df[,1]==k & t.df[,2]==l)     
  #  }
  #}
  mat=matrix(table(t.df[,1],t.df[,2]),ncol=3) #as.matrix(table(t.df[,1],t.df[,2]))
  
  mat=mat/rowSums(mat)

  # get the variation and mean in each state
  means=c(mean(obs[my.states==1]),mean(obs[my.states==2]),mean(obs[my.states==3]))
  vars= c(var(obs[my.states==1]),var(obs[my.states==2]),var(obs[my.states==3]))
  initp=c(1/3,1/3,1/3) # initial starting probabilities
  #d.set=distributionSet(dis="NORMAL", mean=means,var=vars)
  HMMSet(initProb=initp, transMat=mat, dis="NORMAL", means,vars)
}

# function chops large segments into smaller ones, effectively removing regions without base coverage on the segments
# only data.table and tricks via GenomicRanges package should be sufficient to do this
chop.shop<-function(seg.df,df,chop.th){

  # get gaps between covered bases
  cp=GRanges(seqnames=df$chr,ranges=IRanges(start=df$start,end=df$end),meth.diff=df$meth.diff)
  #start(cp)=start(cp)-chop.th
  #end(cp)  =end(cp)  +chop.th

  gap=IRanges::gaps(IRanges::reduce(cp))
  gap=gap[width(gap)>chop.th,]  # only gaps bigger than chop.th should be considered for removal

  #put segment to GRanges
  seg=GRanges(seqnames=seg.df$chr,ranges=IRanges(start=seg.df$start,end=seg.df$end),strand="*",DataFrame(seg.df[,4:7]) )

  #get segments that contain gaps within
  mat=IRanges::findOverlaps(gap,seg,type = "within",ignore.strand =TRUE)
  seg2=seg[ unique(as.matrix(mat)[,2]),] # the segments overlap with no coverage regions
  seg.df1=seg.df[-unique(as.matrix(mat)[,2]),] # the segments that don't overlap with no coverage regions, will use this at the end

  gap1=gap[unique(as.matrix(mat)[,1]),] # get gaps contained within a segment
  #start(gap1)=start(gap1)-chop.th       # arrange gap length using chop.th 
  #end(gap1)  =end(gap1)  +chop.th
  
  # remove gaps from segments ( have to do the for loop because setdiff merges adjacent bases to one region)
  seg3=IRanges::setdiff(seg2[elementMetadata(seg2)$hmm.state== 1,],gap1)
  for(i in 2:length(unique(seg.df$hmm.state)))
  {
      seg3=c(seg3,IRanges::setdiff(seg2[elementMetadata(seg2)$hmm.state==i,],gap1))
  }
  
  # get hmm state of chopped windows
  mat=IRanges::findOverlaps(seg3,seg,type = "within",ignore.strand = TRUE)
  #mat=IRanges::findOverlaps(seg3,seg,ignore.strand = TRUE)

  mat=as.matrix(mat)
  #nrow(mat)
  #seg3[(1:length(seg3) )[! 1:length(seg3)  %in% mat[,1] ],]
    
  my.states=elementMetadata(seg)[mat[,2],4] # get the states per new segment

  # get average.meth diff etc for chopped windows
  #start(cp)=start(cp)+chop.th
  #end(cp)  =end(cp)  -chop.th
  mat=IRanges::findOverlaps(seg3,cp,ignore.strand = TRUE) # overlap CpGs with segments
  mat=as.matrix(mat)

  # get average meth values
  dt=data.table(row.id=mat[,1], meth.diff=elementMetadata(cp)[mat[,2],] )
  sum.dt=dt[,list( num.cov.base=length(meth.diff),
                   avg.meth.diff=mean(meth.diff),
                   med.meth.diff=median(meth.diff)),by=list(row.id)]

  # create the new data frame
  new.seg.df=cbind( as.data.frame(seg3)[unique(mat[,1]), -c(4,5)],as.data.frame(sum.dt)[,2:4],hmm.state=my.states )
  names(new.seg.df)=names(seg.df)

  res.df=rbind(seg.df1,new.seg.df)
  res.df=res.df[order(res.df$chr,res.df$start),] #order by start

  return(res.df)
}




.segDiffMeth<-function(df,semi.supervised=TRUE,lo.th=-5,hi.th=5,tol=1e-6,nStates=3, chop.th=NULL,per.base=FALSE){
  require(data.table)
  require(RHmm)
  require(methylKit)
  
  # get observations from data.frame
  ####df =methylKit::getData(methylDiff.obj)
  df =df[order(df$chr,df$start),]
  obs=df$meth.diff
  
  if(semi.supervised){
    my.param=initializeParam(obs) # initialize the search parameters by reasonable assumptions
    x=HMMFit(obs,nStates=3,control=list(verbose=0,initPoint=my.param,tol=tol) ) # fit the HMM
    #vit= RHmm::viterbi(x, obs)
    vit.states=c()
    for(my.chr in unique(as.character(df$chr)) )
    {
      vit.states=c(vit.states, RHmm::viterbi(x, obs[df$chr==my.chr])$states ) # get state predictions per chr and concat them 
    }


  }else{
    x=HMMFit(obs,nStates=nStates,control=list(verbose=0,init="KMEANS",tol=tol) );
    vit.states=c()
    for(my.chr in unique(as.character(df$chr)) )
    {
      vit.states=c(vit.states, RHmm::viterbi(x, obs[df$chr==my.chr])$states )
    }
  }

  # print summary
  cat("\nLog-likelihood of the model: ",x$LLH,"\n")
  cat("\nDoes it converge ? : ",x$convergence,"\n\n")
  for(i in sort(unique(vit.states)) ){
    cat("Summary of methylation info for state", i,": \n")
    print( summary( obs[vit.states==i]) )
  }

  if(per.base){return(data.frame(id=df[,1],meth.diff=obs,hmm.state=vit.states))}
  
  # create segment ids
  seg.num=1:length(rle(vit.states)$lengths) 
  seg.len=rle(vit.states)$lengths
  seg.id=rep(seg.num,seg.len)
  
  dfs=data.table( cbind(df,vit.states,seg.id) ) # make a data table to be summarized

  # summarize data on predicted segments
  sum.df=dfs[,list(chr=unique(chr),
                   start=min(start),
                   end  =max(end),
                   num.cov.base=length(start),
                   avg.meth.diff=mean(meth.diff),
                   med.meth.diff=median(meth.diff),
                   hmm.state=unique(vit.states) ),by=list(chr,seg.id)]
  res=as.data.frame(sum.df)[,c(-1,-2)] # return data as a data.frame
  
  if(!is.null(chop.th) & is.numeric(chop.th)){

    #stop("chop.th option not supported yet\n")
    res=chop.shop(seg.df=res,df=df,chop.th=chop.th)
  }
  

  res

}

.segDiffMeth.train<-function(df,semi.supervised=TRUE,tol=1e-6,nStates=3 ){
  require(RHmm)
  require(methylKit)

  # get observations from data.frame
  ####df =methylKit::getData(methylDiff.obj)
  df =df[order(df$chr,df$start),]
  obs=df$meth.diff

  if(semi.supervised){
    my.param=initializeParam(obs) # initialize the search parameters by reasonable assumptions
    x=HMMFit(obs,nStates=3,control=list(verbose=0,initPoint=my.param,tol=tol) ) # fit the HMM
  }else{
    x=HMMFit(obs,nStates=nStates,control=list(verbose=0,init="KMEANS",tol=tol) )
  }
  # print summary
  cat("\nLog-likelihood of the model: ",x$LLH,"\n")
  cat("\nDoes it converge ? : ",x$convergence,"\n\n")

  return(x)
}

.segDiffMeth.predict<-function(df, hmm.fit, chop.th=NULL, per.base=FALSE){
  require(data.table)
  require(RHmm)
  require(methylKit)

  # get observations from data.frame
  ####df =methylKit::getData(methylDiff.obj)
  df =df[order(df$chr,df$start),]
  obs=df$meth.diff

  vit.states=c()
  for(my.chr in unique(as.character(df$chr)) )
  {
    vit.states=c(vit.states, RHmm::viterbi(x, obs[df$chr==my.chr])$states ) # get state predictions per chr and concat them 
  }

  # print summary
  for(i in sort(unique(vit.states)) ){
    cat("Summary of methylation info for state", i,": \n")
    print( summary( obs[vit.states==i]) )
  }

  if(per.base){return(data.frame(id=df[,1],meth.diff=obs,hmm.state=vit.states))}

  # create segment ids
  seg.num=1:length(rle(vit.states)$lengths)
  seg.len=rle(vit.states)$lengths
  seg.id=rep(seg.num,seg.len)

  dfs=data.table( cbind(df,vit.states,seg.id) ) # make a data table to be summarized
  # summarize data on predicted segments
  sum.df=dfs[,list(chr=unique(chr),
                   start=min(start),
                   end  =max(end),
                   num.cov.base=length(start),
                   avg.meth.diff=mean(meth.diff),
                   med.meth.diff=median(meth.diff),
                   hmm.state=unique(vit.states) ),by=list(chr,seg.id)]
  res=as.data.frame(sum.df)[,c(-1,-2)] # return data as a data.frame

  if(!is.null(chop.th) & is.numeric(chop.th)){

    #stop("chop.th option not supported yet\n")
    res=chop.shop(seg.df=res,df=df,chop.th=chop.th)
  }

  res
}

####################### S4 CLASS  ########################################################################################################################################## 
#' An S4 class for holding RHmm fit data from the \code{segDiffMeth.train} 
#'
#' This object stores the HMMFit output data that is generated by the \code{segDiffMeth.train} function utilizing the HMMFit function from the \code{RHmm} package.
#'
#' @name HMMFitClass-class
#' @rdname HMMFitClass-class
#' @export
setClass("HMMFitClass", contains="list")

####################### S4 FUNCTIONS  ##########################################################################################################################################  

#' function that segments the genome based on differential methylation scores per base
#'
#' The functions uses a 3 state HMM to segment the genome into hypo,hyper and not differentially methylated regions in the semi-supervised setting. Otherwise there is no limit on states. 
#' @param obj a \code{methylDiff} or \code{methylBase} object with two groups from methylKit package
#' @param semi.supervised if TRUE HMM parameters are initialized with reasonable assumptions (default:TRUE)
#' @param lo.th Needed for paramater initialization of semi-supervised HMM (ignored when semi.supervised=FALSE). Every base below this % methylation difference threshold are assumed to be hypo-methylated (default: -5)
#' @param hi.th Needed for paramater initialization of semi-supervised HMM (ignored when semi.supervised=FALSE). Every base above this % methylation difference threshold are assumed to be hyper-methylated (default: 5)
#' @param tol tolerance threshold for HMMfit() function (Baum-Welch algorithm) (default:1e-6)
#' @param nStates number of states to be used in the unsupervised algorithm, only taken into account when semi.supervised=FALSE
#' @param per.base If TRUE state predictions for each base will be returned instead of concatanated segment coordinates (default:FALSE)
#' @param chop.th a base-pair threshold used to chop large segments with low read coverage to smaller segments. If NULL, this operation is not invoked (default:NULL)
#' @return a data.frame with the region boundaries, HMM state for the regions, average diff.meth value ofthe region and number of covered bases. If per.base=TRUE, the function will return hmm.state per base. This option is provided just for convenience.
#'
#' @usage segDiffMeth(obj,semi.supervised=TRUE,lo.th=-5,hi.th=5,tol=1e-6,nStates=3,chop.th=NULL,per.base=FALSE)
#'
#' @note if you get "Error: protect(): protection stack overflow" error, save your data and start R with " --max-ppsize=400000" option or higher
#' @author Altuna Akalin
#' @examples
#' library(methylKit)
#' data(methylKit)
#' seg.df=segDiffMeth(methylDiff.obj,semi.supervised=TRUE,lo.th=-5,hi.th=5,tol=1e-6 )
#' seg.df2=segDiffMeth(methylDiff.obj,semi.supervised=FALSE,lo.th=-5,hi.th=5,tol=1e-6,nStates=4 )
#' @rdname segDiffMeth
setGeneric("segDiffMeth", function(obj,semi.supervised=TRUE,lo.th=-5,hi.th=5,tol=1e-6,nStates=3,chop.th=NULL,per.base=FALSE) standardGeneric("segDiffMeth"))

#' @rdname segDiffMeth
#' @aliases segDiffMeth,methylDiff-method
setMethod("segDiffMeth", "methylDiff",
                    function(obj, semi.supervised, lo.th, hi.th, tol, nStates, chop.th, per.base){
                      df=methylKit::getData(obj)
                      .segDiffMeth(df,semi.supervised=semi.supervised,lo.th=lo.th,hi.th=hi.th,tol=tol, nStates=nStates, chop.th=chop.th, per.base=per.base)
})

#' @rdname segDiffMeth
#' @aliases segDiffMeth,methylBase-method
setMethod("segDiffMeth", "methylBase",
                    function(obj, semi.supervised, lo.th, hi.th, tol, nStates, chop.th, per.base){
                      if(unique(obj@treatment) != 2){stop("There must be two groups in the treatment vector, segDiffMeth will only work with two groups")}

                      df=methylKit::getData(obj)

                      # get the indices for numCs and numTs in each set
                      set1.Cs=obj@numCs.index[obj@treatment==1]
                      set2.Cs=obj@numCs.index[obj@treatment==0]
                      set1.Ts=obj@numTs.index[obj@treatment==1]
                      set2.Ts=obj@numTs.index[obj@treatment==0]

                      if(length(set1.Cs)>1){
                        pm.meth1 = 100*rowSums(df[,set1.Cs])/rowSums(df[,set1.Cs-1],na.rm=TRUE) # get weigthed means
                      }else{
                        pm.meth1 = 100*(df[,set1.Cs]/df[,set1.Cs-1]) # get % methylation
                      }

                      if(length(set2.Cs)>1){
                        pm.meth2 = 100*rowSums(df[,set2.Cs])/rowSums(df[,set2.Cs-1],na.rm=TRUE) # get weigthed means
                      }else{
                        pm.meth2 = 100*(df[,set2.Cs]/df[,set2.Cs-1]) 
                      }
                      
                      pm.mean.diff=pm.meth1-pm.meth2
                      df=cbind(df[,1:7],meth.diff=pm.mean.diff)
                      .segDiffMeth(df,semi.supervised=semi.supervised,lo.th=lo.th,hi.th=hi.th,tol=tol,nStates=nStates,chop.th=chop.th, per.base=per.base)

})

#' Differential methylatin segmentation HMM fitting 
#'
#' The functions uses a 3 state HMM to fit the genome into hypo,hyper and not differentially methylated regions in the semi-supervised setting. Otherwise there is no limit on states. 
#' @param obj a \code{methylDiff} or \code{methylBase} object with two groups from methylKit package
#' @param semi.supervised if TRUE HMM parameters are initialized with reasonable assumptions (default:TRUE)
#' @param tol tolerance threshold for HMMfit() function (Baum-Welch algorithm) (default:1e-6)
#' @param nStates number of states to be used in the unsupervised algorithm, only taken into account when semi.supervised=FALSE
#' @return a \code{HMMFitClass} object
#'
#' @usage segDiffMeth.train(obj,semi.supervised=TRUE,tol=1e-6,nStates=3)
#'
#' @note if you get "Error: protect(): protection stack overflow" error, save your data and start R with " --max-ppsize=400000" option or higher
#' @author Altuna Akalin and Sheng Li
#' @examples
#' library(methylKit)
#' data(methylKit)
#' seg.fit=segDiffMeth.train(methylDiff.obj,semi.supervised=TRUE,tol=1e-6 )
#' seg.fit2=segDiffMeth.train(methylDiff.obj,semi.supervised=TRUE,tol=1e-6, nStates=4)
#' @rdname segDiffMeth.train
setGeneric("segDiffMeth.train", function(obj,semi.supervised=TRUE,tol=1e-6,nStates=3) standardGeneric("segDiffMeth.train"))

#' @rdname segDiffMeth.train
#' @aliases segDiffMeth.train,methylDiff-method
setMethod("segDiffMeth.train", "methylDiff",
                    function(obj, semi.supervised,tol, nStates){
                      df=methylKit::getData(obj)
                      .segDiffMeth.train(df,semi.supervised=semi.supervised,tol=tol, nStates=nStates)
})

#' @rdname segDiffMeth.train
#' @aliases segDiffMeth.train,methylBase-method
setMethod("segDiffMeth.train", "methylBase",
                    function(obj, semi.supervised, tol, nStates){
                      if(unique(obj@treatment) != 2){stop("There must be two groups in the treatment vector, segDiffMeth will only work with two groups")}

                      df=methylKit::getData(obj)

                      # get the indices for numCs and numTs in each 
                      set1.Cs=obj@numCs.index[obj@treatment==1]
                      set2.Cs=obj@numCs.index[obj@treatment==0]
                      set1.Ts=obj@numTs.index[obj@treatment==1]
                      set2.Ts=obj@numTs.index[obj@treatment==0]

                      if(length(set1.Cs)>1){
                        pm.meth1 = 100*rowSums(df[,set1.Cs])/rowSums(df[,set1.Cs-1],na.rm=TRUE) # get weigthed means
                      }else{
                        pm.meth1 = 100*(df[,set1.Cs]/df[,set1.Cs-1]) # get % methylation
                      }

                      if(length(set2.Cs)>1){
                        pm.meth2 = 100*rowSums(df[,set2.Cs])/rowSums(df[,set2.Cs-1],na.rm=TRUE) # get weigthed means
                      }else{
                        pm.meth2 = 100*(df[,set2.Cs]/df[,set2.Cs-1])
                      }

                      pm.mean.diff=pm.meth1-pm.meth2
                      df=cbind(df[,1:7],meth.diff=pm.mean.diff)
                      .segDiffMeth.train(df,semi.supervised=semi.supervised,tol=tol,nStates=nStates,)

})

#' Predict differetial methylation states using \code{HMMFitClass} object and segments the genome 
#'
#' The functions uses a \code{HMMFitClass} object to segment the genome into hypo,hyper and not differentially methylated regions in the semi-supervised setting. Otherwise there is no limit on states. 
#' @param obj a \code{methylDiff} or \code{methylBase} object with two groups from methylKit package
#' @param hmm.fit a \code{HMMFitClass} object from function \code{segDiffMeth.train}
#' @param per.base If TRUE state predictions for each base will be returned instead of concatanated segment coordinates (default:FALSE)
#' @param chop.th a base-pair threshold used to chop large segments with low read coverage to smaller segments. If NULL, this operation is not invoked (default:NULL)
#' @return a data.frame with the region boundaries, HMM state for the regions, average diff.meth value ofthe region and number of covered bases. If per.base=TRUE, the function will return hmm.state per base. This option is provided just for convenience.

#'
#' @usage segDiffMeth.predict(obj,hmm.fit,chop.th=600, per.base=FALSE)
#'
#' @note if you get "Error: protect(): protection stack overflow" error, save your data and start R with " --max-ppsize=400000" option or higher
#' @author Altuna Akalin and Sheng Li
#' @examples
#' library(methylKit)
#' data(methylKit)
#' seg.fit=segDiffMeth.train(methylDiff.obj,semi.supervised=TRUE,tol=1e-6, nStates=4)
#' seg.df=segDiffMeth.predict(methylDiff.obj, seg.fit, chop.th=600, per.base=FALSE)
#' @rdname segDiffMeth.predict
setGeneric("segDiffMeth.predict", function(obj,hmm.fit,chop.th=NULL, per.base=FALSE) standardGeneric("segDiffMeth.predict"))

#' @rdname segDiffMeth.predict 
#' @aliases segDiffMeth.predict,methylDiff-method,HMMFitClass-method
setMethod("segDiffMeth.predict", signature(obj="methylDiff", hmm.fit="HMMFitClass"),
                    function(obj, hmm.fit, chop.th, per.base){
                      df=methylKit::getData(obj)
                      .segDiffMeth.predict(df,hmm.fit, chop.th=chop.th, per.base=per.base)
})

#' @rdname segDiffMeth.predict
#' @aliases segDiffMeth.predict,methylBase-method,HMMFitClass-method
setMethod("segDiffMeth.predict",  signature(obj="methylBase", hmm.fit="HMMFitClass"),
                    function(obj, hmm.fit, chop.th, per.base){
                      if(unique(obj@treatment) != 2){stop("There must be two groups in the treatment vector, segDiffMeth will only work with two groups")}

                      df=methylKit::getData(obj)

                      # get the indices for numCs and numTs in each set
                      set1.Cs=obj@numCs.index[obj@treatment==1]
                      set2.Cs=obj@numCs.index[obj@treatment==0]
                      set1.Ts=obj@numTs.index[obj@treatment==1]
                      set2.Ts=obj@numTs.index[obj@treatment==0]

                      if(length(set1.Cs)>1){
                        pm.meth1 = 100*rowSums(df[,set1.Cs])/rowSums(df[,set1.Cs-1],na.rm=TRUE) # get weigthed means
                      }else{
                        pm.meth1 = 100*(df[,set1.Cs]/df[,set1.Cs-1]) # get % methylation
                      }

                      if(length(set2.Cs)>1){
                        pm.meth2 = 100*rowSums(df[,set2.Cs])/rowSums(df[,set2.Cs-1],na.rm=TRUE) # get weigthed means
                      }else{
                        pm.meth2 = 100*(df[,set2.Cs]/df[,set2.Cs-1])
                      }

                      pm.mean.diff=pm.meth1-pm.meth2
                      df=cbind(df[,1:7],meth.diff=pm.mean.diff)
                      .segDiffMeth.predict(df,hmm.fit,chop.th=chop.th, per.base=per.base)

})


#' prints out BED track file for HMM segments
#'
#' The functions prints out color coded BED track file from the data frame returned by \code{segDiffMeth} or \code{segMeth}
#' @param seg.df data.frame resulting from \code{segDiffMeth} or \code{segMeth} function
#' @param state.names a character vector for the names of the states, should follow the numerical order of the states
#' @param cols color names for states
#' @param file.name name of the bed file to be written out. If NULL, the BED track will be returned as a data frame, nothing will be written out
#' @param description description string for the track line of the BED file
#' @return bed file written out
#'
#' @usage seg2bed(seg.df,state.names=c("none","hyper","hypo"),cols=c("gray","magenta","green"),file.name=NULL,description )
#'
#' @author Altuna Akalin
#' @examples
#' library(methylKit)
#' data(methylKit)
#' seg.df=segDiffMeth(methylDiff.obj,semi.supervised=TRUE,lo.th=-5,hi.th=5,tol=1e-6 )
#' seg2bed(seg.df,file.name="example.seg.bed")
#'
#' @rdname seg2bed
setGeneric("seg2bed", function(seg.df,state.names=c("none","hyper","hypo"),cols=c("gray","magenta","green"),file.name=NULL,description="my bed track") standardGeneric("seg2bed"))

#' @rdname seg2bed
#' @aliases seg2bed,data.frame-method
setMethod("seg2bed", "data.frame", function(seg.df, state.names, cols ,file.name, description){

  # check if it is a proper data.frame from segMeth or segDiffMeth
  if( ! all(names(seg.df) %in% c("chr","start","end","num.cov.base","avg.meth.diff","med.meth.diff","hmm.state","avg.meth","med.meth" ))  ){
    stop("Provided data frame doesn't look like the output of segMeth() or segDiffMeth() functions")}
  if( ncol(seg.df) != 7 ){
    stop("Provided data frame doesn't look like the output of segMeth() or segDiffMeth() functions")}
  if( length(state.names) != length(unique(seg.df$hmm.state)) ){
    stop("length of state.names should be equal to number of HMM states")}
  if( length(cols) != length(unique(seg.df$hmm.state)) ){
    stop("length of cols vector should be equal to number of HMM states")}

  
  # initialize rgb color vector rcols
  rcols=rep(paste(col2rgb(cols[1]),collapse=","),nrow(seg.df))
  for(i in 2:length(unique(seg.df$hmm.state)) ){
    rcols[seg.df$hmm.state==i]=paste(col2rgb(cols[i]),collapse=",")
  }

  # initialize state names vector state.vec
  state.vec=rep(state.names[1],nrow(seg.df))
  for(i in 2:length(unique(seg.df$hmm.state)) ){
    state.vec[seg.df$hmm.state==i]= state.names[i]
  }
  
  df=cbind( seg.df[,1:3],name=paste("num.cov.base",seg.df[,4],state.vec,sep="_"),
           avg.meth.diff=round(seg.df[,5],2),strand=rep(".",nrow(seg.df)),thickStart=rep(0,nrow(seg.df)),thickEnd=rep(0,nrow(seg.df)),cols=rcols)
  if(is.null(file.name)){return(df)}
  track.line=paste("track type=bed name='",file.name,"' description='",description,"' itemRgb='On' ",sep="")
  #unlink(file.name)                      
  cat(track.line,"\n",file=file.name)
  df[,2]=df[,2]-1 # convert them to 0-base coordinates before writing out
  write.table(df,file=file.name,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t",append=TRUE)
}) 



#' Filter segments based on DMC content, HMM state and number of covered bases
#'
#' The functions filters segments based on their DMC content and HMM state. The segment data frame should be obtained by segDiffMeth() function
#' @param seg.df data.frame resulting from \code{segDiffMeth}  function
#' @param methylDiff.obj a \code{methylDiff} object used to create seg.df via \code{segDiffMeth}  function
#' @param DMC.fraction threshold for fraction of DMCs over number of covered bases. Only segments having DMC fraction over this threshold will be retained
#' @param difference cutoff for absolute value of methylation change between test and control when extracting DMCs from \code{methylDiff} object (default:25)
#' @param qvalue cutoff for qvalue of differential methylation when extracting DMCs from \code{methylDiff} object (default:0.01)
#' @param covered.bases threshold for number of covered bases per segment. only segments with number of bases covered above this threshold will be retained
#' @param retain.state vector containining HMM state numbers. Only these states in this vector will be retained. If NULL, every state will be retained.
#' @return a data frame which is a subset of seg.df input data.frame
#'
#' @usage segFilter(seg.df,methylDiff.obj,DMC.fraction=0.25,difference=25,qvalue=0.01, covered.bases=NULL, retain.state=NULL)
#'
#' @author Altuna Akalin
#'
#' @rdname segFilter
setGeneric("segFilter", function(seg.df,methylDiff.obj,DMC.fraction=0.25,difference=25,qvalue=0.01, covered.bases=NULL, retain.state=NULL) standardGeneric("segFilter"))

#' @rdname segFilter
#' @aliases segFilter,data.frame,methylDiff-method
setMethod("segFilter", c("data.frame","methylDiff"), function(seg.df,methylDiff.obj,DMC.fraction ,difference ,qvalue , covered.bases , retain.state){

  # check if it is a proper data.frame from segMeth or segDiffMeth
  if( ! all(names(seg.df) %in% c("chr","start","end","num.cov.base","avg.meth.diff","med.meth.diff","hmm.state" ))  ){
    stop("Provided data frame doesn't look like the output of segDiffMeth() function, segFilter() only works with the output of segDiffMeth() not with segMeth()")}
  if( ncol(seg.df) != 7 ){
    stop("Provided data frame doesn't look like the output of segDiffMeth() function, it has to have 7 columns")}

  
  dmc=get.methylDiff(methylDiff.obj,difference=difference,qvalue=qvalue,type="all")
  g.dmc=as(dmc,"GRanges")

  #put segment to GRanges
  seg=GRanges(seqnames=seg.df$chr,ranges=IRanges(start=seg.df$start,end=seg.df$end),strand="*" )

  #get segments that contain DMCs within
  mat=IRanges::findOverlaps(g.dmc,seg,type = "within",ignore.strand =TRUE)
  mat=as.matrix(mat)

  # get the number of DMCs in each segment
  tab=table(mat[,2])
  num.dmc1=data.frame(row.id=as.numeric(names(tab)),Freq=as.vector(tab)) ; nrow(num.dmc1); #num.dmc1$Var1=as.numeric( num.dmc1$Var1 ) # get number of DMCs per segments 
  num.dmc2=data.frame(row.id=(1:nrow(seg.df))[!  (1:nrow(seg.df)) %in% num.dmc1$row.id ], Freq=0   ) ; nrow(num.dmc2) # give 0 DMCs for segments that don't have any DMCs
  num.dmc =rbind(num.dmc1, num.dmc2)
  num.dmc =num.dmc[order(num.dmc$row.id),]

  if(nrow(num.dmc) != nrow(seg.df)){stop("Error when getting DMC frequencies of segments, something is wrong")}
  
  freq.DMC=num.dmc$Freq/seg.df$num.cov.base

  res.df=seg.df[freq.DMC > DMC.fraction,]

  if( (! is.null(covered.bases) ) & is.numeric(covered.bases) ){ res.df=res.df[res.df$num.cov.base > covered.bases ,] }

  if( (! is.null(retain.state) ) & is.numeric(retain.state)){ res.df=res.df[res.df$hmm.state %in% retain.state ,]  }

  
  res.df

})



#' Convert data frame from segDiffMeth or segMeth to GRanges
#'
#' The function converts segment data frame from segDiffMeth or segMeth to GRanges object.
#' @param seg.df data.frame resulting from \code{segDiffMeth} or \code{segMeth} function
#' @return a GRanges object for the segments
#'
#' @usage seg2GRanges(seg.df)
#'
#' @author Altuna Akalin
#'
#' @rdname seg2GRanges
setGeneric("seg2GRanges", function(seg.df) standardGeneric("seg2GRanges"))

#' @aliases seg2GRanges,data.frame-method
setMethod("seg2GRanges",  "data.frame",  function(seg.df) {

  # check if it is a proper data.frame from segMeth or segDiffMeth
  check.seg.df(seg.df)
  
  seg=GRanges(seqnames=seg.df$chr,ranges=IRanges(start=seg.df$start,end=seg.df$end),strand="*",DataFrame(seg.df[,4:7]) )
  seg
})

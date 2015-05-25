require(qvalue)

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
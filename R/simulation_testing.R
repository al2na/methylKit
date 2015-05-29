set.seed(100)

# 1. test run and histograms w/o covariates
result <- data_sim(replicates=2,treatment=c(0,1),increase=0.3)
head(result)

h11 <- hist(result$TCols1/50,breaks=50)
h12 <- hist(result$TCols2/50,breaks=50)
plot(h11, col=rgb(0,0,1,1/4),main=NULL)
plot(h12, col=rgb(1,0,0,1/4),add=TRUE,main=NULL)
title("simulated methylation counts w/o covariates")
legend("top",legend=c("control","treated"),fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

# correlation of CCols and TCols
smoothScatter(result$TCols1,result$CCols1, main=paste("Correlation of TCols and CCols for untreated sample 1:",
                                                      round(cor(result$TCols1,result$CCols1),2)))
smoothScatter(result$TCols2,result$CCols2,main=paste("Correlation of TCols and CCols for treated sample 2:",
                                                     round(cor(result$TCols2,result$CCols2),2)))

# 2. test run and histograms with 1 covariate
covariates=data.frame(age=c(20,60,30,50))
result2 <- data_sim(replicates=4,treatment=c(0,0,1,1),increase=0.3,covariates=covariates)
head(result2)
h21 <- hist(result2$TCols1/50,breaks=50)
h22 <- hist(result2$TCols2/50,breaks=50)
h23 <- hist(result2$TCols3/50,breaks=50)
h24 <- hist(result2$TCols4/50,breaks=50)

plot(h21, col=rgb(0,0,1,1/4),main=NULL)
plot(h22, col=rgb(1,0,0,1/4),add=TRUE,main=NULL)
title("simulated methylation counts with covariates: samples 1 and 2")
legend("top",legend=c("control,age=20","control,age=60"),fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

plot(h21, col=rgb(0,0,1,1/4),main=NULL)
plot(h23, col=rgb(1,0,0,1/4),add=TRUE,main=NULL)
title("simulated methylation counts with covariates: samples 1 and 3")
legend("top",legend=c("control,age=20","treated,age=30"),fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))

plot(h21, col=rgb(0,0,1,1/4),main=NULL)
plot(h24, col=rgb(1,0,0,1/4),add=TRUE,main=NULL)
title("simulated methylation counts with covariates: samples 1 and 4")
legend("top",legend=c("control,age=20","treated,age=50"),fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))


# correlation of CCols and TCols
smoothScatter(result2$TCols1,result2$CCols1, main=paste("Correlation of TCols and CCols for untreated sample 1:",
                                                      round(cor(result2$TCols1,result2$CCols1),2)))
smoothScatter(result2$TCols2,result2$CCols2,main=paste("Correlation of TCols and CCols for untreated sample 2:",
                                                     round(cor(result2$TCols2,result2$CCols2),2)))
smoothScatter(result2$TCols3,result2$CCols3,main=paste("Correlation of TCols and CCols for treated sample 3:",
                                                       round(cor(result2$TCols3,result2$CCols3),2)))
smoothScatter(result2$TCols4,result2$CCols4,main=paste("Correlation of TCols and CCols for treated sample 4:",
                                                       round(cor(result2$TCols4,result2$CCols4),2)))

# load test data(-head), extract treatment and data.frame-part
data(methylKit)
df <- getData(methylBase.obj)
treatment <- methylBase.obj@treatment
covariates=NULL

res <- logRegDF(df,treatment,covariates,overdispersion="none",effect="wmean",test="Chisq",mc.cores=4)
head(res)
head(methylDiff.obj)

covariates=data.frame(batch=c("a","b","a","a"))
res2 <- logRegDF(df,treatment,covariates,overdispersion="none",effect="wmean",test="Chisq",mc.cores=4)
res3 <- logRegDF(df,treatment,covariates,overdispersion="none",effect="wmean",test="F",mc.cores=4)
res4 <- logRegDF(df,treatment,covariates,overdispersion="MN",effect="wmean",test="F",mc.cores=4)
res5 <- logRegDF(df,treatment,covariates,overdispersion="MN",effect="wmean",test="F",adjust="bonferroni",mc.cores=4)

# => Issue 1: when length(w) == nprm, phi is Inf and the tests produce NA (division by 0)

covariates=data.frame(batch=c("a","b","a","a"),age=c(20,40,30,40))
res6 <- logRegDF(df,treatment,covariates,overdispersion="none",effect="wmean",test="Chisq",mc.cores=4)
res7 <- logRegDF(df,treatment,covariates,overdispersion="MN",effect="wmean",test="Chisq",mc.cores=4)
res8 <- logRegDF(df,treatment,covariates,overdispersion="MN",effect="wmean",test="F",mc.cores=4)

# check Fisher-test

#counts <- c(20,19,22,12,19,17,10,37)
#covariates=data.frame(batch=c("a","b","a","a"))
#covariates=data.frame(batch=c("a","b","a","a"),age=c(20,40,30,40))
#vars <- as.data.frame(cbind(treatment,covariates))
#formula <-as.formula(paste("~ ", paste(colnames(vars), collapse= "+")))
#overdispersion="MN"
#test="F"
#effect="wmean"
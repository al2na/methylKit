#library("methylKit")
context("calculateDiffMeth checks")

data("methylKit")


dbdir <- "methylDB"
methylBaseDB <- makeMethylDB(methylBase.obj,dbdir = dbdir)
mc.cores <- 1


test_that("check if default calculateDiffMeth works", {
  expect_is(
    calculateDiffMeth(
      methylBase.obj,
      save.db = TRUE,
      dbdir = dbdir,
      suffix = "mdiff",
      mc.cores = mc.cores
    ),
    'methylDiffDB'
  )
  expect_is(calculateDiffMeth(methylBaseDB, suffix = "mDiff",mc.cores = mc.cores),
            'methylDiffDB')
  expect_equal(
    calculateDiffMeth(methylBase.obj, save.db = FALSE, mc.cores = mc.cores),
    calculateDiffMeth(methylBaseDB, save.db = FALSE, mc.cores = mc.cores)
  )
  expect_equal(
    getData(calculateDiffMeth(
      methylBase.obj,
      save.db = TRUE,
      dbdir = dbdir,
      mc.cores = mc.cores
    )),
    getData(calculateDiffMeth(
      methylBaseDB,
      save.db = TRUE,
      dbdir = dbdir,
      mc.cores = mc.cores
    )))
})


## create artificial single samples methylBase 
suppressWarnings(single.methylBase<- dataSim(replicates=1,sites=1000,treatment=c(1),
                         sample.ids=c("test1")))
## create artificial single group methylBase 
singleGroup.methylBase <- dataSim(replicates=4,sites=1000,treatment=c(0,0,0,0),
                                  sample.ids=c("test1","test2","ctrl1","ctrl2"))

## check if covariates+intercept+treatment more than replicates 
## 2 + 1 + 1 >= 4   
covariates.failed=data.frame(age=c(30,80,30,80),
                             gender=c("male","female","male","female"))

test_that("check if calculateDiffMeth errors", {
  expect_error(calculateDiffMeth(single.methylBase, mc.cores = mc.cores))
  expect_error(calculateDiffMeth(singleGroup.methylBase, mc.cores = mc.cores))
  expect_error(
    calculateDiffMeth(methylBase.obj,
                      mc.cores = mc.cores,
                      covariates = covariates.failed)
  )
  
})

# pool samples in each group
pooled.methylBase=pool(methylBase.obj,sample.ids=c("test","control"))

test_that("check if calculateDiffMeth with pooled methylBase works", {
  expect_is(calculateDiffMeth(pooled.methylBase), 'methylDiff')
})

# Covariates and overdispersion control:
# generate a methylBase object with age as a covariate
set.seed(123)
covariates=data.frame(age=c(30,80,30,80))
sim.methylBase<-dataSim(replicates=4,sites=1000,treatment=c(1,1,0,0),
                        covariates=covariates,
                        sample.ids=c("test1","test2","ctrl1","ctrl2"))

test_that("check if calculateDiffMeth with covariates works", {
  expect_is(suppressWarnings(
    calculateDiffMeth(
      sim.methylBase,
      covariates = covariates,
      overdispersion = "MN",
      test = "Chisq",
      mc.cores = 1
    )
  ),
  'methylDiff')
})

## set more than two groups
methylBase.threeGroup <- methylBase.obj 
getTreatment(methylBase.threeGroup) <- c(1,1,2,3)

test_that("check if calculateDiffMeth with three treatment groups works", {
  expect_is(suppressWarnings(calculateDiffMeth(methylBase.threeGroup)), 'methylDiff')
})


test_that("check different arguments for calculateDiffMeth", {
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "SLIM",test = "fast.fisher"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "holm",test = "F"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "hochberg",test = "Chisq"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "hommel",test = "midPval"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "bonferroni",overdispersion = "none"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "BH",overdispersion = "MN"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "BY",overdispersion = "shrinkMN"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "fdr",effect = "wmean"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "none",effect = "mean"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,adjust = "qvalue",effect = "predicted"), 'methylDiff')
  expect_is(calculateDiffMeth(methylBase.obj,slim = FALSE ,weighted.mean = FALSE), 'methylDiff')
  
})


unlink("tests/testthat/methylDB",recursive = TRUE)
unlink(dbdir,recursive = TRUE)

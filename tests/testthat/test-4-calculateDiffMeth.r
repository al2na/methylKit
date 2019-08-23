#library("methylKit")
context("calculateDiffMeth checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

myobj=methRead( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

mydblist = suppressMessages(methRead( file.list,
                                      sample.id=list("t1","t2","c1","c2"),assembly="hg18",
                                      pipeline="amp",treatment=c(1,1,0,0),dbtype = "tabix",dbdir="methylDB"))

# unite function
methidh=unite(myobj)
methidh2=unite(myobj,min.per.group=1L)
methidh3=pool(methidh,sample.ids = c("test","control"))
methidh4 <- methidh; getTreatment(methidh4) <- c(1,1,2,3)

suppressMessages(methidhDB <- unite(mydblist))
suppressMessages(methidh2DB <- unite(mydblist,min.per.group=1L))
suppressMessages(methidh3DB <- pool(methidhDB,sample.ids = c("test","control")))
suppressMessages(methidh4DB <- makeMethylDB(obj = methidh4,dbdir="methylDB"))

# differential methylation
myDiff =calculateDiffMeth(methidh)
myDiff2=calculateDiffMeth(methidh2)
myDiff3=calculateDiffMeth(methidh3)
myDiff4=calculateDiffMeth(methidh4)

suppressMessages(myDiffDB <- calculateDiffMeth(methidhDB))
suppressMessages(myDiff2DB <- calculateDiffMeth(methidh2DB))
suppressMessages(myDiff3DB <- calculateDiffMeth(methidh3DB))
suppressMessages(myDiff4DB <- calculateDiffMeth(methidh4DB))

hypo <- getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
hyper <- getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")

suppressMessages(hypoDB <- getMethylDiff(myDiffDB,difference=25,qvalue=0.01,type="hypo"))
suppressMessages(hyperDB <- getMethylDiff(myDiffDB,difference=25,qvalue=0.01,type="hyper"))


myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
suppressMessages(mydblist2 <- reorganize(mydblist,sample.ids=c("t1","c2"),treatment=c(1,0) ))

test_that("check if reorganize works with methylRawList", {
  expect_is(myobj2, 'methylRawList')
})

test_that("check if reorganize works with methylRawListDB", {
  expect_is(mydblist2, 'methylRawListDB')
})


test_that("check if calculateDiffMeth output is a methylDiff object", {
  expect_is(myDiff, 'methylDiff')
})

test_that("check if calculateDiffMeth output is a methylDiffDB object", {
  expect_is(myDiffDB, 'methylDiffDB')
})


test_that("check if calculateDiffMeth output from unite(...,min.per.group=1) is a methylDiff object", {
  expect_is(myDiff2, 'methylDiff')
})

test_that("check if calculateDiffMeth output from unite(...,min.per.group=1) is a methylDiffDB object", {
  expect_is(myDiff2DB, 'methylDiffDB')
})


test_that("check if calculateDiffMeth output from pooled methylBase is a methylDiff object", {
  expect_is(myDiff3, 'methylDiff')
})

test_that("check if calculateDiffMeth output from pooled methylBaseDB is a methylDiffDB object", {
  expect_is(myDiff3DB, 'methylDiffDB')
})


test_that("check if calculateDiffMeth output with three treatment groups is a methylDiff object", {
  expect_is(myDiff4, 'methylDiff')
})

test_that("check if calculateDiffMeth output from three treatment groups is a methylDiffDB object", {
  expect_is(myDiff4DB, 'methylDiffDB')
})




test_that("check getting hypo/hyper meth works with methylDiff", {
  expect_is(hypo, 'methylDiff')
  expect_is(hyper, 'methylDiff')
})

test_that("check getting hypo/hyper meth works with methylDiffDB", {
  expect_is(hypoDB, 'methylDiffDB')
  expect_is(hyperDB, 'methylDiffDB')
})

test_that("check different arguments for calculateDiffMeth", {
  expect_is(calculateDiffMeth(methidh,adjust = "SLIM",mc.cores = 1,test = "fast.fisher"), 'methylDiff')
  expect_is(calculateDiffMeth(methidh,adjust = "SLIM",mc.cores = 1,test = "F"), 'methylDiff')
  expect_is(calculateDiffMeth(methidh,adjust = "SLIM",mc.cores = 1,test = "Chisq"), 'methylDiff')
  expect_is(calculateDiffMeth(methidh,adjust = "SLIM",mc.cores = 1,test = "midPval"), 'methylDiff')
  expect_is(calculateDiffMeth(methidh3,adjust = "SLIM",mc.cores = 1,test = "fast.fisher"), 'methylDiff')
  expect_is(calculateDiffMeth(methidh3,adjust = "SLIM",mc.cores = 1,test = "F"), 'methylDiff')
  expect_is(calculateDiffMeth(methidh3,adjust = "SLIM",mc.cores = 1,test = "Chisq"), 'methylDiff')
  expect_is(calculateDiffMeth(methidh3,adjust = "SLIM",mc.cores = 1,test = "midPval"), 'methylDiff')
  
})


unlink("tests/testthat/methylDB",recursive = TRUE)
library("methylKit")
context("calculateDiffMeth and annotate.WithGenicPart checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
file.list

myobj=read( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

mydblist = suppressMessages(read( file.list,
                                  sample.id=list("t1","t2","c1","c2"),assembly="hg18",
                                  pipeline="amp",treatment=c(1,1,0,0),dbtype = "tabix",dbdir="methylDB"))

# unite function
methidh=unite(myobj)
methidh2=unite(myobj,min.per.group=1L)

suppressMessages(methidhDB <- unite(mydblist))
suppressMessages(methidh2DB <- unite(mydblist,min.per.group=1L))

# differential methylation
myDiff =calculateDiffMeth(methidh)
myDiff2=calculateDiffMeth(methidh2)

suppressMessages(myDiffDB <- calculateDiffMeth(methidhDB))
suppressMessages(myDiff2DB <- calculateDiffMeth(methidh2DB))

hypo <- get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
hyper <- get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")

suppressMessages(hypoDB <- get.methylDiff(myDiffDB,difference=25,qvalue=0.01,type="hypo"))
suppressMessages(hyperDB <- get.methylDiff(myDiffDB,difference=25,qvalue=0.01,type="hyper"))


# load annotation
gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))
cpg.obj=read.feature.flank(system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit"),feature.flank.name=c("CpGi","shores"))

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


test_that("check getting hypo/hyper meth works with methylDiff", {
    expect_is(hypo, 'methylDiff')
    expect_is(hyper, 'methylDiff')
})

test_that("check getting hypo/hyper meth works with methylDiffDB", {
  expect_is(hypoDB, 'methylDiffDB')
  expect_is(hyperDB, 'methylDiffDB')
})


test_that("annotate.WithGenicParts with methylDiff", {
    expect_is(annotate.WithGenicParts(myDiff,gene.obj),
        'annotationByGenicParts')
})

test_that("annotate.WithGenicParts with methylDiffDB", {
  expect_is(annotate.WithGenicParts(myDiffDB,gene.obj),
            'annotationByGenicParts')
})


test_that("annotate.WithGenicParts with methylDiff", {
    expect_is(annotate.WithFeature.Flank(myDiff,cpg.obj$CpGi,cpg.obj$shores),
        'annotationByFeature')
})

test_that("annotate.WithGenicParts with methylDiffDB", {
  expect_is(annotate.WithFeature.Flank(myDiffDB,cpg.obj$CpGi,cpg.obj$shores),
            'annotationByFeature')
})
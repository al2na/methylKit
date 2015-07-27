library("methylKit")
context("calculateDiffMethDSS checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
file.list

myobj=read( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

# unite function
methidh=unite(myobj)
methidh2=unite(myobj,min.per.group=1L)

# differential methylation
myDiff =calculateDiffMethDSS(methidh)
myDiff2=calculateDiffMethDSS(methidh2)

# load annotation
gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))
cpg.obj=read.feature.flank(system.file("extdata", "cpgi.hg18.bed.txt", package = "methylKit"),feature.flank.name=c("CpGi","shores"))

myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )

test_that("check if reorganize works", {
    expect_that(myobj2, 
        is_a('methylRawList'))
})

test_that("check if calculateDiffMeth output is a methylDiff object", {
    expect_that(myDiff, 
        is_a('methylDiff'))
})



test_that("check if calculateDiffMeth output from unite(...,min.per.group=1) is a methylDiff object", {
    expect_that(myDiff2, 
        is_a('methylDiff'))
})

test_that("check getting hypo/hyper meth works", {
    expect_that(get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hypo"),
        is_a('methylDiff'))
    expect_that(get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hyper"),
        is_a('methylDiff'))
})

test_that("annotate.WithGenicParts", {
    expect_that(annotate.WithGenicParts(myDiff,gene.obj),
        is_a('annotationByGenicParts'))
})


test_that("annotate.WithGenicParts", {
    expect_that(annotate.WithFeature.Flank(myDiff,cpg.obj$CpGi,cpg.obj$shores),
        is_a('annotationByFeature'))
})
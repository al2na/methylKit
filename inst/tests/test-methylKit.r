library(methylKit)

print("check library")
file.list=list( system.file("tests", "test1.myCpG.txt", package = "methylKit"),
                system.file("tests", "test2.myCpG.txt", package = "methylKit"),
                system.file("tests", "control1.myCpG.txt", package = "methylKit"),
                system.file("tests", "control2.myCpG.txt", package = "methylKit") )

myobj=read( file.list,
           sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

# unite function
t=unite(myobj,destrand=T)
methidh=unite(myobj)

# differential methylation
myDiff=calculateDiffMeth(methidh)

# load annotation 
gene.obj=read.transcript.features(system.file("tests", "refseq.hg18.bed.txt", package = "methylKit"))


# testthat
context("methylKit functions test")

test_that("If all the test files are in the file.list",{
    expect_that(length(file.list), equals(4))
})

test_that("read function works", {
    expect_that(myobj, is_a('methylRawList'))
})
test_that("getMethylationStats on myobj[[2]] works", {
    expect_that (getMethylationStats(myobj[[2]],plot=F,both.strands=F), prints_text('methylation statistics per basesummary:   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    0.00   20.00   82.79   63.17   94.74  100.00 percentiles:       0%       10%       20%       30%       40%       50%       60%       70%       80%       90%       95%       99%     99.5%     99.9%      100%   0.00000   0.00000   0.00000  48.38710  70.00000  82.78556  90.00000  93.33333  96.42857 100.00000 100.00000 100.00000 100.00000 100.00000 100.00000'))
})

test_that("unite function works", {
    expect_that(t, is_a('methylBase'))
    expect_that(methidh, is_a('methylBase'))
})

test_that("clusterSamples function works", {
    expect_that(clusterSamples(methidh, dist="correlation", method="ward", plot=FALSE), is_a('hclust'))
})

test_that("PCASamples function works", {
    expect_that(PCASamples(methidh), is_a('summary.princomp'))
})

test_that("calculateDiffMeth function works", {
    expect_that(myDiff, is_a('methylDiff'))
})

test_that("annotate.WithGenicParts", {
    expect_that(annotate.WithGenicParts(myDiff,gene.obj), is_a('annotationByGenicParts'))
})

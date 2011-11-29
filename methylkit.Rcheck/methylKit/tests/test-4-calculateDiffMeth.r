library("methylKit")
context("calculateDiffMeth and annotate.WithGenicPart checks")

file.list=list( system.file("data", "test1.myCpG.txt", package = "methylKit"),
                system.file("data", "test2.myCpG.txt", package = "methylKit"),
                system.file("data", "control1.myCpG.txt", package = "methylKit"),
                system.file("data", "control2.myCpG.txt", package = "methylKit") )
file.list

myobj=read( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

# unite function
t=unite(myobj,destrand=T)
methidh=unite(myobj)

# differential methylation
myDiff=calculateDiffMeth(methidh)

# load annotation
gene.obj=read.transcript.features(system.file("data", "refseq.hg18.bed.txt", package = "methylKit"))

test_that("check if there are 4 test files in the file.list",{
    expect_that(length(file.list),
        equals(4))
})

test_that("check if calculateDiffMeth output is a methylDiff object", {
    expect_that(myDiff, 
        is_a('methylDiff'))
})

test_that("annotate.WithGenicParts", {
    expect_that(annotate.WithGenicParts(myDiff,gene.obj),
        is_a('annotationByGenicParts'))
})

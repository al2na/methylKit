library("methylKit")
context("calculateDiffMeth and annotate.WithGenicPart checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
file.list

myobj=read( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

# unite function
t=unite(myobj,destrand=T)
methidh=unite(myobj)

# differential methylation
myDiff=calculateDiffMeth(methidh)

# load annotation
gene.obj=read.transcript.features(system.file("extdata", "refseq.hg18.bed.txt", package = "methylKit"))



test_that("check if calculateDiffMeth output is a methylDiff object", {
    expect_that(myDiff, 
        is_a('methylDiff'))
})

test_that("annotate.WithGenicParts", {
    expect_that(annotate.WithGenicParts(myDiff,gene.obj),
        is_a('annotationByGenicParts'))
})

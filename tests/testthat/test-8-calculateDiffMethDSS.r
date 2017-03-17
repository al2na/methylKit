library("methylKit")
context("calculateDiffMethDSS checks")

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )
file.list

myobj=methRead( file.list,
                sample.id=list("test1","test2","ctrl1","ctrl2"),assembly="hg18",
            pipeline="amp",treatment=c(1,1,0,0))

# unite function
methidh=unite(myobj)
methidh2=unite(myobj,min.per.group=1L)

# differential methylation
myDiff =calculateDiffMethDSS(methidh)
myDiff2=calculateDiffMethDSS(methidh2)


myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )

test_that("check if reorganize works", {
    expect_that(myobj2, 
        is_a('methylRawList'))
})

test_that("check if calculateDiffMeth output is a methylDiff object", {
    expect_that(myDiff, 
        is_a('methylDiff'))
})



test_that(paste("check if calculateDiffMeth output from", 
          "unite(...,min.per.group=1) is a methylDiff object"), {
    expect_that(myDiff2, 
        is_a('methylDiff'))
})

test_that("check getting hypo/hyper meth works", {
    expect_that(getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo"),
        is_a('methylDiff'))
    expect_that(getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper"),
        is_a('methylDiff'))
})


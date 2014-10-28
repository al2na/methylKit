context("regionCounts checks")



myobj=read( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                sample.id="test1",assembly="hg18",pipeline="amp")


library(GenomicRanges)

promoters=GRanges(seqnames=c("chr21","chr21"),
                             ranges=IRanges(start=c(17806094,1011979),end=c(17808094,10121796)),strand=c("+","-") )

test_that("check that regionCounts on methylRaw works", {
    expect_that( regionCounts(myobj,promoters),
        is_a('methylRaw'))
})

data(methylKit)
test_that("check that tileMethylCounts on methylBase works", {
  expect_that( tileMethylCounts(object=methylBase.obj,win.size=1000,step.size=1000,cov.bases=0) ,
  is_a('methylBase'))
})

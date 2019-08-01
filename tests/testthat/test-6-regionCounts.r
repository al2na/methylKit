context("regionCounts checks")



myobj=methRead( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                sample.id="test1",assembly="hg18",pipeline="amp")
myobjDB=suppressMessages(methRead( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                         sample.id="test1",assembly="hg18",pipeline="amp",
                         dbtype = "tabix",dbdir="methylDB"))


suppressMessages(library(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE))

promoters=GRanges(seqnames=c("chr21","chr21"),
                             ranges=IRanges(start=c(17806094,1011979),end=c(17808094,10121796)),strand=c("+","-") )

test_that("check that regionCounts on methylRaw works", {
    expect_is( regionCounts(myobj,promoters),
        'methylRaw')
    expect_is( suppressMessages(regionCounts(myobj,promoters,save.db = T,dbdir="methylDB",suffix="promoters")),
        'methylRawDB')
})

test_that("check that regionCounts on methylRawDB works", {
  expect_is( regionCounts(myobjDB,promoters),
               'methylRawDB')
  expect_is( regionCounts(myobjDB,promoters,save.db = F),
             'methylRaw')
})

data(methylKit)
test_that("check that tileMethylCounts on methylRaw works", {
  expect_is( tileMethylCounts(object=myobj,win.size=1000,step.size=1000,cov.bases=0) ,
  'methylRaw')
})

test_that("check that tileMethylCounts on methylRawDB works", {
  expect_is( tileMethylCounts(object=myobjDB,win.size=1000,step.size=1000,cov.bases=0) ,
               'methylRawDB')
})

unlink("tests/testthat/methylDB",recursive = TRUE)
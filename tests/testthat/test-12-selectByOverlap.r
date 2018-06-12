context("selectByOverlap checks")
data(methylKit)

file.list=list( system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )

methylRawListDB.obj=methRead(file.list,
                         sample.id=list("test1","test2","ctrl1","ctrl2"),
                         assembly="hg18",treatment=c(1,1,0,0),
                         dbtype = "tabix",dbdir = "methylDB")

methylBaseDB.obj=unite(methylRawListDB.obj)

methylDiffDB.obj = calculateDiffMeth(methylBaseDB.obj)

# define the windows of interest as a GRanges object, this can be any set
# of genomic locations
library(GenomicRanges)
my.win=GRanges(seqnames="chr21",
ranges=IRanges(start=seq(from=9764513,by=10000,length.out=20),width=5000) )

test_that("check if selectByOverlap works for methylDB objects." ,{
  expect_is(selectByOverlap(methylRawListDB.obj[[1]],my.win),"methylRaw")
  expect_is(selectByOverlap(methylRawListDB.obj,my.win),"methylRawList")
  expect_is(selectByOverlap(methylBaseDB.obj,my.win),"methylBase")
  expect_is(selectByOverlap(methylDiffDB.obj,my.win),"methylDiff")
})

test_that("check if selectByOverlap works for normal methyl objects." ,{
  expect_is(selectByOverlap(methylRawList.obj[[1]],my.win),"methylRaw")
  expect_is(selectByOverlap(methylRawList.obj,my.win),"methylRawList")
  expect_is(selectByOverlap(methylBase.obj,my.win),"methylBase")
  expect_is(selectByOverlap(methylDiff.obj,my.win),"methylDiff")
})

context("methSeg and methSeg2bed checks")
#library("methylKit")
data("methylKit")

tileRaw <- tileMethylCounts(methylRawList.obj[[1]])
tileBase <- tileMethylCounts(methylBase.obj)
tileDiff <- calculateDiffMeth(tileBase)

test_that("check if methSeg works for methylRaw or methylDiff with resolution region" ,{
  expect_is(methSeg(tileRaw,diagnostic.plot = FALSE),"GRanges")
  expect_is(methSeg(tileDiff,diagnostic.plot = FALSE),"GRanges")
})

test_that("check if methSeg works for methylRaw or methylDiff with resolution base " ,{
  expect_is(methSeg(methylRawList.obj[[1]],diagnostic.plot = FALSE),"GRanges")
  expect_is(methSeg(methylDiff.obj,diagnostic.plot = FALSE),"GRanges")
})

test_that("check if methSeg errors for a single ranged methylRaw" ,{
  expect_error(methSeg(methylRawList.obj[[1]][1],diagnostic.plot = FALSE))
})

test_that("check if methSeg returns error for methylBase objects" ,{
  expect_error(methSeg(methylBase.obj,diagnostic.plot = FALSE))
})

seg <- methSeg(tileRaw,diagnostic.plot = FALSE)
methSeg2bed(segments = seg, filename = "test.bed")

test_that("check if methSeg2bed returns bed file", {
  expect_true(file.exists("test.bed"))
})

unlink("test.bed")
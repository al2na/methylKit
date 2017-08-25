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
  expect_warning(methraw.gr <- methSeg(methylRawList.obj[[1]],diagnostic.plot = FALSE))
  expect_is(methraw.gr,"GRanges")
  expect_is(methSeg(methylDiff.obj,diagnostic.plot = FALSE),"GRanges")
})

test_that("check if methSeg errors for a single ranged methylRaw" ,{
  expect_error(methSeg(methylRawList.obj[[1]][1],diagnostic.plot = FALSE))
})

test_that("check if methSeg returns error for methylBase objects" ,{
  expect_error(methSeg(methylBase.obj,diagnostic.plot = FALSE))
})


gr = as(methylRawList.obj[[1]],"GRanges")
mcols(gr)$meth=100*gr$numCs/gr$coverage

test_that("check if methSeg works for GRanges object" ,{
  expect_is(methSeg(gr,diagnostic.plot = FALSE),"GRanges")
})

gr_2 <- gr
mcols(gr_2) <- NULL
test_that("check if methSeg errors for GRanges object without meta columns" ,{
  expect_error(methSeg(gr_2,diagnostic.plot = FALSE))
})


gr.list <- lapply(methylRawList.obj, FUN = function(obj){ 
  obj= as(obj,"GRanges")
  mcols(obj)$meth=100*obj$numCs/obj$coverage
  obj = obj[,"meth"]
})
large.gr <- do.call("c",gr.list)

res1 <- methSeg(large.gr)
res2 <- methSeg(large.gr,join.neighbours = TRUE)

test_that("check if joining neighbours works" ,{
  expect_false(all(rle(res1$seg.group)$lengths == 1))
  expect_true(all(rle(res2$seg.group)$lengths == 1))
})


seg <- methSeg(tileRaw,diagnostic.plot = FALSE)
methSeg2bed(segments = seg, filename = "test.bed")

test_that("check if methSeg2bed returns bed file", {
  expect_true(file.exists("test.bed"))
})

unlink("test.bed")
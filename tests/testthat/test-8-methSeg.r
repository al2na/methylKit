context("methSeg and methSeg2bed checks")
#library("methylKit")
data("methylKit")

tileRaw <- tileMethylCounts(methylRawList.obj[[1]])
tileBase <- tileMethylCounts(methylBase.obj)
tileDiff <- calculateDiffMeth(tileBase)

MethDiff_multiChrom <- methylDiff.obj
MethDiff_multiChrom$chr <- rep(paste0("chr",1:5), c(nrow(MethDiff_multiChrom)-10,1,3,1,5))

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

test_that("check if methSeg warns for on single ranged chrom" ,{
  expect_warning(MethDiff_multiChrom.gr <- methSeg(MethDiff_multiChrom, diagnostic.plot = FALSE))
  expect_is(MethDiff_multiChrom.gr,"GRanges")
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

test_that("check if initialization works" ,{
  expect_is(methSeg(methylDiff.obj,diagnostic.plot = FALSE,verbose=FALSE),"GRanges")
  expect_is(methSeg(methylDiff.obj,diagnostic.plot = FALSE,verbose=FALSE,initialization=list(subset = sample(1:13,9,replace = FALSE))),"GRanges")
  expect_error(methSeg(methylDiff.obj,diagnostic.plot = FALSE,verbose=FALSE,initialize.on.subset = 0.6))
  expect_is(methSeg(methylDiff.obj,diagnostic.plot = FALSE,verbose=FALSE,initialize.on.subset = 10),"GRanges")
})

seg <- methSeg(tileRaw,diagnostic.plot = FALSE)
methSeg2bed(segments = seg, filename = "test.bed")

test_that("check if methSeg2bed returns bed file", {
  expect_true(file.exists("test.bed"))
})

segments <- GRanges(seqnames = c("chr2","chr2","chr2","chrM","chrX"),
                    ranges = IRanges(
                      start = c(52445739,58194119,58582091,4848,65761515),
                      end = c(52453651,58211798,58601283,15960,65769140)),
                    strand = "*", ID = "meth", num.mark = c(4,6,8,11,5),
                    seg.mean = 0, startRow = c(918976,988985,992438,1,561162),
                    endRow = c(918980,988991,992446,12,561167),
                    seg.group = 1, name = 1)
methSeg2bed(segments = segments,filename = "test2.bed")


test_that("check if methSeg2bed returns bed file", {
  expect_true(file.exists("test2.bed"))
})

unlink("test.bed")
context("test getMethylationStats")

data("methylKit")

mydblist <- makeMethylDB(methylRawList.obj)

outDir = "figDir"
dir.create(outDir,showWarnings = FALSE)



pngFile <- tempfile(tmpdir = outDir,fileext = ".png")
png(pngFile)
getMethylationStats(methylRawList.obj[[1]],plot=T,both.strands=F)
dev.off()
pngFile2 <- tempfile(tmpdir = outDir,fileext = ".png")
png(pngFile2)
getMethylationStats(methylRawList.obj[[1]],plot=T,both.strands=T)
dev.off()

test_that("getMethylationStats on methylRaw works", {
  expect_output(getMethylationStats(methylRawList.obj[[1]],plot=F,both.strands=F),
                'methylation statistics per base')
  expect_output(getMethylationStats(methylRawList.obj[[1]],plot=F,both.strands=T),
                'methylation statistics per base')
  expect_true(file.size(pngFile) > 0)
  expect_true(file.size(pngFile2) > 0)
})


pngFile <- tempfile(tmpdir = outDir,fileext = ".png")
png(pngFile)
getMethylationStats(mydblist[[2]],plot=T,both.strands=F)
dev.off()
pngFile2 <- tempfile(tmpdir = outDir,fileext = ".png")
png(pngFile2)
getMethylationStats(mydblist[[2]],plot=T,both.strands=T)
dev.off()

test_that("getMethylationStats on methylRawDB works", {
  expect_output(getMethylationStats(mydblist[[2]],plot=F,both.strands=F),
                'methylation statistics per base')
  expect_output(getMethylationStats(mydblist[[2]],plot=F,both.strands=T),
                'methylation statistics per base')
  expect_true(file.size(pngFile) > 0)
  expect_true(file.size(pngFile2) > 0)
})

unlink(dirname(getDBPath(mydblist)),recursive = TRUE)
unlink(outDir,recursive = TRUE)

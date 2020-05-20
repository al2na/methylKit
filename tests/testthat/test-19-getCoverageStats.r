context("test getCoverageStats")

data("methylKit")

mydblist <- makeMethylDB(methylRawList.obj)

outDir = "figDir"
dir.create(outDir,showWarnings = FALSE)



pngFile <- tempfile(tmpdir = outDir,fileext = ".png")
png(pngFile)
getCoverageStats(methylRawList.obj[[1]],plot=T,both.strands=F)
dev.off()
pngFile2 <- tempfile(tmpdir = outDir,fileext = ".png")
png(pngFile2)
getCoverageStats(methylRawList.obj[[1]],plot=T,both.strands=T)
dev.off()

test_that("getCoverageStats on methylRaw works", {
  expect_output(getCoverageStats(methylRawList.obj[[1]],plot=F,both.strands=F),
                'read coverage statistics per base')
  expect_output(getCoverageStats(methylRawList.obj[[1]],plot=F,both.strands=T),
                'read coverage statistics per base')
  expect_true(file.size(pngFile) > 0)
  expect_true(file.size(pngFile2) > 0)
})


pngFile <- tempfile(tmpdir = outDir,fileext = ".png")
png(pngFile)
getCoverageStats(mydblist[[2]],plot=T,both.strands=F)
dev.off()
pngFile2 <- tempfile(tmpdir = outDir,fileext = ".png")
png(pngFile2)
getCoverageStats(mydblist[[2]],plot=T,both.strands=T)
dev.off()

test_that("getCoverageStats on methylRawDB works", {
  expect_output(getCoverageStats(mydblist[[2]],plot=F,both.strands=F),
                'read coverage statistics per base')
  expect_output(getCoverageStats(mydblist[[2]],plot=F,both.strands=T),
                'read coverage statistics per base')
  expect_true(file.size(pngFile) > 0)
  expect_true(file.size(pngFile2) > 0)
})

unlink(dirname(getDBPath(mydblist)),recursive = TRUE)
unlink(outDir,recursive = TRUE)

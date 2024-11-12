context("test export to bedgraph files")

data(methylKit)

# remove the file
outDir <- "bedDir"
dir.create(outDir,showWarnings = FALSE)
dbdir <- "methylDB"

# getting a bedgraph file from a methylDiff object containing differential
# methylation percentages
outFile <- tempfile(pattern = "mdiff.",tmpdir = outDir,fileext = ".bed")
bedgraph(methylDiff.obj, file.name=outFile, col.name="meth.diff",
        unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")
tt <- rtracklayer::import.bedGraph(outFile)
outFile2 <- tempfile(pattern = "mdiff.",tmpdir = outDir,fileext = ".bed")


test_that("export methylDiff worked", {
  expect_true(file.exists(outFile))
  expect_equal(length(tt),nrow(methylDiff.obj))
  expect_error(bedgraph(methylDiff.obj,col.name = "meth.diff",
         log.transform = TRUE,negative = TRUE,file.name = outFile2))
  expect_is(bedgraph(methylDiff.obj,col.name = "meth.diff"),
            class = "data.frame")
  expect_error(bedgraph(methylDiff.obj,col.name = "coverage"))
})

# getting a bedgraph file from a methylRaw object containing raw
#methylation values
outFile <- tempfile(pattern = "mRaw",tmpdir = outDir,fileext = ".bed")
bedgraph(methylRawList.obj[[1]], file.name=outFile, col.name="perc.meth",
        unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")
outFile2 <- tempfile(pattern = "mRaw",tmpdir = outDir,fileext = ".bed")
bedgraph(methylRawList.obj[[1]], file.name=outFile2, col.name="perc.meth",
        unmeth=TRUE,log.transform=FALSE,negative=FALSE,add.on="")
tt <- rtracklayer::import.bedGraph(outFile)
tt2 <- rtracklayer::import.bedGraph(outFile2)

test_that("export methylRaw worked", {
  expect_true(file.exists(outFile))
  expect_equal(length(tt),nrow(methylRawList.obj[[1]]))
  expect_equal(length(tt2),2)
  expect_is(bedgraph(methylRawList.obj[[1]],col.name = "coverage",
                     log.transform = TRUE,negative = TRUE),
            class = "data.frame")
  expect_is(bedgraph(methylRawList.obj[[1]],col.name = "perc.meth"),
            class = "data.frame")
  expect_error(bedgraph(methylRawList.obj[[1]],col.name = "meth.diff"))
})

# getting a bedgraph file from a methylRawList object containing raw
#methylation values
outFile <- tempfile(pattern = "mRawL",tmpdir = outDir,fileext = ".bed")
bedgraph(methylRawList.obj, file.name=outFile, col.name="perc.meth",
        unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")
outFile2 <- tempfile(pattern = "mRawL",tmpdir = outDir,fileext = ".bed")
bedgraph(methylRawList.obj, file.name=outFile2, col.name="perc.meth",
        unmeth=TRUE,log.transform=FALSE,negative=FALSE,add.on="")
tt <- rtracklayer::import.bedGraph(outFile)
tt2 <- rtracklayer::import.bedGraph(outFile2)

test_that("export methylRaw worked", {
  expect_true(file.exists(outFile))
  expect_equal(length(tt),length(methylRawList.obj))
  expect_equal(length(tt2),2*length(methylRawList.obj))
  expect_is(bedgraph(methylRawList.obj,col.name = "coverage",
                     log.transform = TRUE,negative = TRUE),
            class = "list")
  expect_error(bedgraph(methylRawList.obj,col.name = "meth.diff"))
})

# getting a bedgraph file from a methylDiffDB object containing differential
# methylation percentages
outFile <- tempfile(pattern = "mdiffDB.",tmpdir = outDir,fileext = ".bed")
methylDiffDB.obj <- makeMethylDB(methylDiff.obj, dbdir)
bedgraph(methylDiffDB.obj, file.name=outFile, col.name="meth.diff",
        unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")
tt <- rtracklayer::import.bedGraph(outFile)
outFile2 <- tempfile(pattern = "mdiffDB.",tmpdir = outDir,fileext = ".bed")

test_that("export methylDiffDB worked", {
  expect_true(file.exists(outFile))
  expect_equal(length(tt),methylDiffDB.obj@num.records)
  expect_error(bedgraph(methylDiffDB.obj,col.name = "meth.diff",
         log.transform = TRUE,negative = TRUE,file.name = outFile2))
  expect_is(bedgraph(methylDiffDB.obj,col.name = "meth.diff"),
            class = "data.frame")
  expect_error(bedgraph(methylDiffDB.obj,col.name = "coverage"))
})

# remove the file
unlink(outDir,recursive = TRUE)
unlink(dbdir,recursive = TRUE)

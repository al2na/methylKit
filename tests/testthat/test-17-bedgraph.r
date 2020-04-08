context("test export to bedgraph files")

data(methylKit)

# remove the file
unlink("test.bed")
unlink("test2.bed")

# getting a bedgraph file from a methylDiff object containing differential
# methylation percentages
bedgraph(methylDiff.obj, file.name="test.bed", col.name="meth.diff",
        unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

test_that("export worked", {
  expect_true(file.exists("test.bed"))
})

tt <- rtracklayer::import.bedGraph("test.bed")

test_that("import worked", {
  expect_equal(length(tt),nrow(methylDiff.obj))
})


# getting a bedgraph file from a methylBase object containing percent
#methylation values
bedgraph(methylRawList.obj[[1]], file.name="test2.bed", col.name="perc.meth",
        unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

test_that("export worked", {
  expect_true(file.exists("test2.bed"))
})

tt2 <- rtracklayer::import.bedGraph("test2.bed")

test_that("import worked", {
  expect_equal(length(tt2),nrow(methylRawList.obj[[1]]))
})


# remove the file
unlink("test.bed")
unlink("test2.bed")

